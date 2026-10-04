//! The encoder: PCM in, Vorbis I packets out (and Ogg pages, with
//! [`OggWriter`]).
//!
//! Two block sizes (256 and 2048) chosen per block by a transient
//! detector; a psychoacoustic mask per block and channel; a
//! floor 1 curve fitted to that mask; residues quantised against the floor
//! and coded with residue 2 over all channels, the front and rear pairs
//! square-polar coupled; codebooks designed from modelled symbol
//! statistics and sent in the setup header. Variable bit rate, steered by
//! one quality value.

mod books;
mod psy;

use std::io::Write;

use crate::bits::BitWriter;
use crate::codebook::Codebook;
use crate::error::{Error, Result};
use crate::floor1::Floor1;
use crate::header::{Comments, Floor, Identification, Setup, xiph_lacing};
use crate::mdct::Mdct;
use crate::ogg::PacketWriter;
use crate::residue::Residue;
use crate::tables::FLOOR1_INVERSE_DB;
use crate::window::Windows;

use books::{CLASSES, FLOOR_SUB_LIMITS, MAX_RESIDUE, PARTITION, RESIDUE_BOOKS};

/// The block sizes the encoder uses.
pub const BLOCKSIZE: [usize; 2] = [256, 2048];

/// Square polar coupling (the inverse of [`decouple`](crate::decouple)):
/// a left/right pair of integer residues to magnitude and angle. Exact:
/// `decouple(couple(l, r)) == (l, r)`.
pub fn couple(l: i32, r: i32) -> (i32, i32) {
    if l.abs() > r.abs() {
        if l > 0 { (l, l - r) } else { (l, r - l) }
    } else if r > 0 {
        (r, l - r)
    } else {
        (r, r - l)
    }
}

/// What to encode.
#[derive(Clone, Debug, PartialEq)]
pub struct EncoderConfig {
    /// Sample rate in Hz, 8000 to 192000.
    pub sample_rate: u32,
    /// Channels, 1 to 8, in Vorbis order (4.3.9): mono; L R; L C R;
    /// FL FR RL RR; FL C FR RL RR; FL C FR RL RR LFE; FL C FR SL SR RC LFE;
    /// FL C FR SL SR RL RR LFE.
    pub channels: u8,
    /// Quality, -1.0 (smallest) to 10.0 (best).
    pub quality: f32,
    /// User comments, as (name, value) pairs.
    pub comments: Vec<(String, String)>,
}

impl Default for EncoderConfig {
    fn default() -> Self {
        EncoderConfig { sample_rate: 44100, channels: 2, quality: 5.0, comments: Vec::new() }
    }
}

/// One encoded audio packet.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EncodedPacket {
    /// The packet.
    pub data: Vec<u8>,
    /// The granule position after it: samples a decoder has returned once
    /// it has decoded this packet (the last packet's is the input length).
    pub granule: i64,
    /// Its block size.
    pub blocksize: usize,
}

/// Per-block encoder statistics.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EncoderStats {
    /// Short and long blocks emitted.
    pub blocks: [u64; 2],
    /// Audio bytes emitted.
    pub bytes: u64,
    /// Bits spent on floors and on residues.
    pub floor_bits: u64,
    /// Bits spent on residues (classifications included).
    pub residue_bits: u64,
    /// Residue partitions coded in each classification.
    pub classes: [u64; 8],
}

/// A Vorbis I encoder.
pub struct Encoder {
    config: EncoderConfig,
    headers: [Vec<u8>; 3],
    setup: Setup,
    floors: [Floor1; 2],
    residues: [Residue; 2],
    psy: [psy::Psy; 2],
    mdct: [Mdct; 2],
    windows: Windows,
    /// Coded bins per block size (the low-pass), and the LFE's.
    lowpass: [usize; 2],
    lfe_lowpass: [usize; 2],
    adjust_db: f64,
    /// Quantiser dead zone, in steps (subtracted from |r| before rounding).
    deadzone: f32,
    channels: usize,
    /// Input per channel, from absolute sample `buf_start`.
    buf: Vec<Vec<f32>>,
    buf_start: i64,
    total: i64,
    finished: bool,
    done: bool,
    /// The next block: its centre, size (decided once known) and the
    /// previous block's size.
    center: i64,
    long: Option<bool>,
    prev_long: bool,
    stats: EncoderStats,
    /// Threads for a batch of blocks; 0 is the machine's count.
    threads: usize,
}

/// `20 log10` of one floor 1 table step.
fn step_db() -> f64 {
    20.0 * (FLOOR1_INVERSE_DB[255] as f64 / FLOOR1_INVERSE_DB[0] as f64).log10() / 255.0
}

impl Encoder {
    /// A new encoder.
    pub fn new(config: EncoderConfig) -> Result<Self> {
        if !(1..=8).contains(&config.channels) {
            return Err(Error::Config(format!("{} channels (1 to 8 are supported)", config.channels)));
        }
        if !(8000..=192_000).contains(&config.sample_rate) {
            return Err(Error::Config(format!("sample rate {} (8000 to 192000 Hz)", config.sample_rate)));
        }
        if !config.quality.is_finite() || !(-1.0..=10.0).contains(&config.quality) {
            return Err(Error::Config(format!("quality {} (-1 to 10)", config.quality)));
        }
        let q = config.quality as f64;
        let ch = config.channels as usize;
        let nyquist = config.sample_rate as f64 / 2.0;
        let lowpass_hz = (13000.0 + 1250.0 * q).clamp(11000.0, 24000.0).min(nyquist);
        let bins = |half: usize, hz: f64| ((hz / nyquist * half as f64).ceil() as usize).clamp(1, half);
        let lowpass = [bins(BLOCKSIZE[0] / 2, lowpass_hz), bins(BLOCKSIZE[1] / 2, lowpass_hz)];
        let lfe_lowpass = [bins(BLOCKSIZE[0] / 2, 150.0), bins(BLOCKSIZE[1] / 2, 150.0)];
        let end = |half: usize, lp: usize| ((lp * ch).div_ceil(PARTITION) * PARTITION).min(half * ch);
        let residue_end = [end(BLOCKSIZE[0] / 2, lowpass[0]), end(BLOCKSIZE[1] / 2, lowpass[1])];
        let setup = books::setup(config.channels, BLOCKSIZE, residue_end);
        let ident = Identification {
            channels: config.channels,
            sample_rate: config.sample_rate,
            bitrate_maximum: 0,
            bitrate_nominal: 0,
            bitrate_minimum: 0,
            blocksize: [BLOCKSIZE[0] as u16, BLOCKSIZE[1] as u16],
        };
        let comments = Comments {
            vendor: format!("rivet-vorbis {}", env!("CARGO_PKG_VERSION")),
            comments: config.comments.iter().map(|(k, v)| format!("{k}={v}")).collect(),
        };
        let headers = [ident.write(), comments.write(), setup.write(config.channels)];
        let floor = |i: usize| match &setup.floors[i] {
            Floor::One(f) => f.clone(),
            Floor::Zero(_) => unreachable!("the encoder uses floor 1"),
        };
        Ok(Encoder {
            floors: [floor(0), floor(1)],
            residues: [setup.residues[0].clone(), setup.residues[1].clone()],
            psy: [
                psy::Psy::new(BLOCKSIZE[0] / 2, config.sample_rate),
                psy::Psy::new(BLOCKSIZE[1] / 2, config.sample_rate),
            ],
            mdct: [Mdct::new(BLOCKSIZE[0]), Mdct::new(BLOCKSIZE[1])],
            windows: Windows::new(BLOCKSIZE),
            lowpass,
            lfe_lowpass,
            adjust_db: 4.0 - 1.6 * q,
            deadzone: (0.2 - 0.01 * q) as f32,
            channels: ch,
            buf: vec![Vec::new(); ch],
            buf_start: 0,
            total: 0,
            finished: false,
            done: false,
            center: 0,
            long: None,
            prev_long: true,
            stats: EncoderStats::default(),
            threads: 0,
            headers,
            setup,
            config,
        })
    }

    /// The three header packets: identification, comment, setup.
    pub fn headers(&self) -> &[Vec<u8>; 3] {
        &self.headers
    }

    /// The headers with Xiph lacing, as Matroska's `CodecPrivate`.
    pub fn codec_private(&self) -> Vec<u8> {
        xiph_lacing([&self.headers[0], &self.headers[1], &self.headers[2]])
    }

    /// The configuration.
    pub fn config(&self) -> &EncoderConfig {
        &self.config
    }

    /// Counts so far.
    pub fn stats(&self) -> &EncoderStats {
        &self.stats
    }

    /// How many threads code a batch of blocks (the blocks one
    /// [`encode`](Self::encode) call completes): 0, the default, is one per
    /// CPU; 1 codes everything on the caller's thread. The packets are the
    /// same byte for byte whatever the count.
    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    /// Encode interleaved samples (full scale ±1.0); returns the packets
    /// they complete.
    pub fn encode(&mut self, interleaved: &[f32]) -> Result<Vec<EncodedPacket>> {
        if self.finished {
            return Err(Error::Config("encode called after finish".into()));
        }
        if !interleaved.len().is_multiple_of(self.channels) {
            return Err(Error::Config("sample count is not a multiple of the channel count".into()));
        }
        for frame in interleaved.chunks(self.channels) {
            for (c, &s) in frame.iter().enumerate() {
                self.buf[c].push(if s.is_finite() { s } else { 0.0 });
            }
        }
        self.total += (interleaved.len() / self.channels) as i64;
        Ok(self.pump())
    }

    /// Encode planar samples, one slice per channel, of equal length.
    pub fn encode_planar(&mut self, planar: &[&[f32]]) -> Result<Vec<EncodedPacket>> {
        if planar.len() != self.channels || planar.iter().any(|c| c.len() != planar[0].len()) {
            return Err(Error::Config("planar input needs one equal-length slice per channel".into()));
        }
        if self.finished {
            return Err(Error::Config("encode called after finish".into()));
        }
        for (b, c) in self.buf.iter_mut().zip(planar) {
            b.extend(c.iter().map(|&s| if s.is_finite() { s } else { 0.0 }));
        }
        self.total += planar[0].len() as i64;
        Ok(self.pump())
    }

    /// End the stream: the remaining packets, the last carrying the input
    /// length as its granule position so decoders trim the padding.
    pub fn finish(&mut self) -> Result<Vec<EncodedPacket>> {
        self.finished = true;
        Ok(self.pump())
    }

    fn sample(&self, c: usize, i: i64) -> f32 {
        if i < self.buf_start || i >= self.total { 0.0 } else { self.buf[c][(i - self.buf_start) as usize] }
    }

    /// Whether `[from, to)` holds an attack: a 64-sample segment whose
    /// high-passed energy jumps well above the four before it.
    fn transient(&self, from: i64, to: i64) -> bool {
        const SEG: i64 = 64;
        let energy = |s: i64| -> f64 {
            let mut e = 0.0;
            for c in 0..self.channels {
                for i in s..s + SEG {
                    let d = (self.sample(c, i) - self.sample(c, i - 1)) as f64;
                    e += d * d;
                }
            }
            e
        };
        let start = from.div_euclid(SEG) * SEG;
        let mut history =
            [energy(start - 4 * SEG), energy(start - 3 * SEG), energy(start - 2 * SEG), energy(start - SEG)];
        let floor = 1e-6 * SEG as f64 * self.channels as f64;
        let mut s = start;
        while s < to {
            let e = energy(s);
            let before = history.iter().sum::<f64>() / 4.0;
            if e > floor && e > 12.0 * before.max(floor * 0.1) {
                return true;
            }
            history.rotate_left(1);
            history[3] = e;
            s += SEG;
        }
        false
    }

    fn pump(&mut self) -> Vec<EncodedPacket> {
        // First the block sequence: sizes follow from the transient
        // detector on the input alone, so it is decided block by block
        // without coding anything.
        struct Planned {
            c: i64,
            long: bool,
            prev_long: bool,
            next_long: bool,
            last: bool,
        }
        let mut plan: Vec<Planned> = Vec::new();
        while !self.done {
            if self.finished && self.total == 0 {
                self.done = true;
                break;
            }
            let c = self.center;
            let horizon = c + BLOCKSIZE[1] as i64 / 4 + 512 + 1024 + 128;
            if !self.finished && self.total < horizon {
                break;
            }
            let long = *self.long.get_or_insert(true);
            let first = c == 0 && self.stats.blocks == [0, 0] && plan.is_empty();
            let long = if first { !self.transient(c - 512, c + 1024) && long } else { long };
            let n = BLOCKSIZE[long as usize] as i64;
            let next_center_long = c + n / 4 + BLOCKSIZE[1] as i64 / 4;
            let next_long = !self.transient(next_center_long - 512, next_center_long + 1024);
            let last = self.finished && c >= self.total;
            plan.push(Planned { c, long, prev_long: self.prev_long, next_long, last });
            if last {
                self.done = true;
                break;
            }
            self.prev_long = long;
            self.long = Some(next_long);
            self.center = c + n / 4 + BLOCKSIZE[next_long as usize] as i64 / 4;
        }
        // Then the blocks, which are independent given the sequence, on
        // up to `threads` threads; packets and statistics in order.
        let threads =
            if self.threads == 0 { std::thread::available_parallelism().map_or(1, usize::from) } else { self.threads };
        let this = &*self;
        let coded = parallel_map(&plan, threads, |p| {
            let mut stats = EncoderStats::default();
            let data = this.encode_block(p.c, p.long, p.prev_long, p.next_long, &mut stats);
            (data, stats)
        });
        let mut out = Vec::with_capacity(plan.len());
        for (p, (data, stats)) in plan.iter().zip(coded) {
            self.stats.floor_bits += stats.floor_bits;
            self.stats.residue_bits += stats.residue_bits;
            for (a, b) in self.stats.classes.iter_mut().zip(stats.classes) {
                *a += b;
            }
            self.stats.blocks[p.long as usize] += 1;
            self.stats.bytes += data.len() as u64;
            let n = BLOCKSIZE[p.long as usize];
            out.push(EncodedPacket { data, granule: if p.last { self.total } else { p.c }, blocksize: n });
        }
        // Keep what the next block and the transient look-back need.
        let keep_from = self.center - BLOCKSIZE[1] as i64 - 1024;
        if keep_from - self.buf_start > 1 << 16 {
            let drop = (keep_from - self.buf_start) as usize;
            for b in self.buf.iter_mut() {
                b.drain(..drop.min(b.len()));
            }
            self.buf_start += drop as i64;
        }
        out
    }

    /// Fit floor 1 to per-bin targets (in floor Y units): each post takes
    /// the mean target of the bins nearest it, held within 2 units of their
    /// minimum; a post predicted within one unit is left uncoded.
    fn fit_floor(floor: &Floor1, target: &[i32]) -> Vec<i32> {
        let half = target.len();
        let xs = &floor.x_list;
        let sorted = floor.sorted_posts();
        let mut post_target = vec![0i32; xs.len()];
        for (idx, &i) in sorted.iter().enumerate() {
            let x = xs[i] as usize;
            let left = if idx > 0 { (xs[sorted[idx - 1]] as usize + x) / 2 } else { 0 };
            let right = if idx + 1 < sorted.len() { (x + xs[sorted[idx + 1]] as usize).div_ceil(2) } else { x + 1 };
            let (lo, hi) = (left.min(half - 1), right.clamp(left + 1, half).max(left.min(half - 1) + 1));
            let region = &target[lo..hi];
            let mean = region.iter().map(|&v| v as f64).sum::<f64>() / region.len() as f64;
            let min = *region.iter().min().expect("non-empty region");
            post_target[i] = (mean.round() as i32).min(min + 2).clamp(0, floor.range() - 1);
        }
        let mut fin = vec![0i32; xs.len()];
        let mut y = vec![0i32; xs.len()];
        fin[0] = post_target[0];
        fin[1] = post_target[1];
        y[0] = fin[0];
        y[1] = fin[1];
        for i in 2..xs.len() {
            let predicted = floor.predict(&fin, i);
            if (post_target[i] - predicted).abs() <= 1 {
                fin[i] = predicted;
                y[i] = 0;
            } else {
                y[i] = floor.encode_value(&fin, i, post_target[i]);
                fin[i] = post_target[i];
            }
        }
        y
    }

    fn write_floor(&self, w: &mut BitWriter, floor: &Floor1, y: &[i32]) {
        let books = &self.setup.codebooks;
        w.write_flag(true);
        let bits = crate::bits::ilog(floor.range() as i64 - 1);
        w.write(y[0] as u32, bits);
        w.write(y[1] as u32, bits);
        let sub = |v: i32| -> usize {
            if v == 0 { 0 } else { 1 + FLOOR_SUB_LIMITS.iter().position(|&l| v < l).expect("value within range") }
        };
        for part in y[2..].chunks(3) {
            let subs: Vec<usize> = part.iter().map(|&v| sub(v)).collect();
            let cval = subs[0] | subs[1] << 2 | subs[2] << 4;
            books[books::BOOK_FLOOR_MASTER].write_entry(w, cval as u32);
            for (&v, &s) in part.iter().zip(&subs) {
                if s > 0 {
                    books[books::BOOK_FLOOR_SUB[s - 1]].write_entry(w, v as u32);
                }
            }
        }
    }

    /// Code the interleaved residue vector `v` as residue 2 does (format 1
    /// over one vector).
    fn write_residue(&self, w: &mut BitWriter, residue: &Residue, v: &[i32], stats: &mut EncoderStats) {
        let books = &self.setup.codebooks;
        let end = (residue.end as usize).min(v.len());
        let parts = end / PARTITION;
        if parts == 0 {
            return;
        }
        let classes: Vec<usize> = (0..parts)
            .map(|p| {
                let m = v[p * PARTITION..(p + 1) * PARTITION].iter().map(|x| x.abs()).max().unwrap_or(0);
                CLASSES.iter().position(|&(max, _)| m <= max).unwrap_or(CLASSES.len() - 1)
            })
            .collect();
        for &c in &classes {
            stats.classes[c] += 1;
        }
        let classbook: &Codebook = &books[books::BOOK_CLASS];
        let per_word = classbook.dimensions as usize;
        let nclass = CLASSES.len();
        for pass in 0..8 {
            let mut p = 0;
            while p < parts {
                if pass == 0 {
                    let mut word = 0;
                    for i in 0..per_word {
                        word = word * nclass + classes.get(p + i).copied().unwrap_or(0);
                    }
                    classbook.write_entry(w, word as u32);
                }
                for _ in 0..per_word {
                    if p >= parts {
                        break;
                    }
                    let (_, passes) = CLASSES[classes[p]];
                    if let Some(&b) = passes.get(pass) {
                        let lattice = &RESIDUE_BOOKS[b];
                        let book = &books[books::BOOK_VQ0 + b];
                        let values: Vec<i32> = v[p * PARTITION..(p + 1) * PARTITION]
                            .iter()
                            .map(|&q| split(q, passes.len())[pass])
                            .collect();
                        for chunk in values.chunks(lattice.dim) {
                            let mut entry = 0usize;
                            for &x in chunk.iter().rev() {
                                let digit = ((x - lattice.min) / lattice.step) as usize;
                                debug_assert!(digit < lattice.values && (x - lattice.min) % lattice.step == 0);
                                entry = entry * lattice.values + digit;
                            }
                            book.write_entry(w, entry as u32);
                        }
                    }
                    p += 1;
                }
            }
        }
    }

    fn encode_block(&self, c: i64, long: bool, prev_long: bool, next_long: bool, stats: &mut EncoderStats) -> Vec<u8> {
        let li = long as usize;
        let n = BLOCKSIZE[li];
        let half = n / 2;
        let ch = self.channels;
        let floor = &self.floors[li];
        let lfe = books::lfe(self.config.channels);
        let step = step_db();
        let mut q = vec![vec![0i32; half]; ch];
        let mut floor_y: Vec<Option<Vec<i32>>> = vec![None; ch];
        for chan in 0..ch {
            let mut block: Vec<f32> = (0..n as i64).map(|i| self.sample(chan, c - n as i64 / 2 + i)).collect();
            self.windows.apply(&mut block, long, prev_long, next_long);
            let mut x = vec![0f32; half];
            self.mdct[li].forward(&block, &mut x);
            let energy: f64 = x.iter().map(|&v| (v as f64) * (v as f64)).sum();
            if energy < 1e-20 {
                continue;
            }
            let lp = if lfe == Some(chan) { self.lfe_lowpass[li] } else { self.lowpass[li] };
            let mask = self.psy[li].mask(&x, self.adjust_db);
            let mut target: Vec<i32> = mask
                .iter()
                .map(|&m| {
                    let amp = (12.0 * m).sqrt();
                    let t = 255.0 + 20.0 * amp.log10() / step;
                    ((t / 2.0).floor() as i32).clamp(0, floor.range() - 1)
                })
                .collect();
            let edge = target[lp - 1];
            for t in target[lp..].iter_mut() {
                *t = edge;
            }
            let y = Self::fit_floor(floor, &target);
            let curve = floor.synthesize(&y, half);
            let mut any = false;
            for k in 0..lp {
                let r = x[k] / curve[k];
                let v = ((r.abs() - self.deadzone).max(0.0).round().copysign(r) as i32)
                    .clamp(-MAX_RESIDUE / 2, MAX_RESIDUE / 2);
                q[chan][k] = v;
                any |= v != 0;
            }
            if any {
                floor_y[chan] = Some(y);
            } else {
                q[chan].fill(0);
            }
        }
        // Coupling (4.3.5 inverted), and the nonzero propagation of 4.3.3.
        let mut used: Vec<bool> = floor_y.iter().map(|f| f.is_some()).collect();
        for &(m, a) in self.setup.mappings[li].coupling.iter() {
            let (m, a) = (m as usize, a as usize);
            if !used[m] && !used[a] {
                continue;
            }
            used[m] = true;
            used[a] = true;
            let (lo, hi) = (m.min(a), m.max(a));
            let (left, right) = q.split_at_mut(hi);
            let (qm, qa) = if m < a { (&mut left[lo], &mut right[0]) } else { (&mut right[0], &mut left[lo]) };
            for (x, y) in qm.iter_mut().zip(qa.iter_mut()) {
                (*x, *y) = couple(*x, *y);
            }
        }
        let mut w = BitWriter::new();
        w.write(0, 1);
        w.write(long as u32, 1);
        if long {
            w.write_flag(prev_long);
            w.write_flag(next_long);
        }
        let before = w.bit_len();
        for y in &floor_y {
            match y {
                Some(y) => self.write_floor(&mut w, floor, y),
                None => w.write_flag(false),
            }
        }
        let floors_end = w.bit_len();
        stats.floor_bits += (floors_end - before) as u64;
        if used.iter().any(|&u| u) {
            let mut v = vec![0i32; half * ch];
            for (chan, qc) in q.iter().enumerate() {
                for (k, &x) in qc.iter().enumerate() {
                    v[k * ch + chan] = x;
                }
            }
            self.write_residue(&mut w, &self.residues[li], &v, stats);
        }
        stats.residue_bits += (w.bit_len() - floors_end) as u64;
        w.into_bytes()
    }
}

/// `f` over `items`, in order, on up to `threads` threads (the calling
/// thread among them), each taking the next item as it finishes one.
fn parallel_map<T: Sync, R: Send>(items: &[T], threads: usize, f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let n = items.len();
    if threads <= 1 || n <= 1 {
        return items.iter().map(f).collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let work = || {
        let mut done = Vec::new();
        loop {
            let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let Some(item) = items.get(i) else { return done };
            done.push((i, f(item)));
        }
    };
    let mut all: Vec<(usize, R)> = std::thread::scope(|s| {
        let helpers: Vec<_> = (1..threads.min(n)).map(|_| s.spawn(work)).collect();
        let mut all = work();
        for h in helpers {
            all.extend(h.join().expect("an encoder thread panicked"));
        }
        all
    });
    all.sort_unstable_by_key(|(i, _)| *i);
    all.into_iter().map(|(_, r)| r).collect()
}

/// A residue value as the passes of its class code it: one value; or a
/// coarse multiple of 33 and the rest; or a multiple of 495, of 33, and
/// the rest.
fn split(q: i32, passes: usize) -> [i32; 3] {
    match passes {
        1 => [q, 0, 0],
        2 => {
            let f = ((q as f64 / 33.0).round() as i32).clamp(-7, 7);
            [33 * f, q - 33 * f, 0]
        }
        _ => {
            let g = ((q as f64 / 495.0).round() as i32).clamp(-31, 31);
            let r = q - 495 * g;
            let f = ((r as f64 / 33.0).round() as i32).clamp(-7, 7);
            [495 * g, 33 * f, r - 33 * f]
        }
    }
}

/// Writes an encoder's stream as an Ogg Vorbis file: the identification
/// header alone on the first page, the comment and setup headers ending
/// the next, then the audio, the last page marked end-of-stream.
pub struct OggWriter<W: Write> {
    pages: PacketWriter<W>,
    held: Option<EncodedPacket>,
}

impl<W: Write> OggWriter<W> {
    /// Start a stream with serial number `serial`, writing `headers`.
    pub fn new(inner: W, headers: &[Vec<u8>; 3], serial: u32) -> Result<Self> {
        let mut pages = PacketWriter::new(inner, serial);
        pages.write_packet(&headers[0], 0, true, false)?;
        pages.write_packet(&headers[1], 0, false, false)?;
        pages.write_packet(&headers[2], 0, true, false)?;
        Ok(OggWriter { pages, held: None })
    }

    /// Queue packets (one is held back, to mark the last end-of-stream).
    pub fn write_packets(&mut self, packets: impl IntoIterator<Item = EncodedPacket>) -> Result<()> {
        for p in packets {
            if let Some(h) = self.held.replace(p) {
                self.pages.write_packet(&h.data, h.granule, false, false)?;
            }
        }
        Ok(())
    }

    /// Write the last packet on an end-of-stream page; returns the writer.
    pub fn finish(mut self) -> Result<W> {
        match self.held.take() {
            Some(h) => self.pages.write_packet(&h.data, h.granule, false, true)?,
            None => self.pages.finish(0)?,
        }
        Ok(self.pages.into_inner())
    }
}

/// Encode interleaved samples to a whole Ogg Vorbis file.
pub fn encode_ogg(config: &EncoderConfig, interleaved: &[f32]) -> Result<Vec<u8>> {
    let mut encoder = Encoder::new(config.clone())?;
    let serial = 0x5256_4f42 ^ (interleaved.len() as u32).rotate_left(7);
    let mut writer = OggWriter::new(Vec::new(), encoder.headers(), serial)?;
    writer.write_packets(encoder.encode(interleaved)?)?;
    writer.write_packets(encoder.finish()?)?;
    writer.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decode::decouple;

    /// Coupling is exact both ways over a range of integers, and the
    /// decoder's float decoupling agrees.
    #[test]
    fn polar_coupling_round_trips() {
        for l in -40..=40 {
            for r in -40..=40 {
                let (m, a) = couple(l, r);
                let (dl, dr) = decouple(m as f32, a as f32);
                assert_eq!((dl as i32, dr as i32), (l, r), "({l}, {r}) -> ({m}, {a})");
            }
        }
        // Decoupling as 4.3.5 writes it, case by case.
        assert_eq!(decouple(3.0, 1.0), (3.0, 2.0));
        assert_eq!(decouple(3.0, -1.0), (2.0, 3.0));
        assert_eq!(decouple(-3.0, 1.0), (-3.0, -2.0));
        assert_eq!(decouple(-3.0, -1.0), (-2.0, -3.0));
        assert_eq!(decouple(0.0, 0.0), (0.0, 0.0));
    }

    #[test]
    fn split_reassembles() {
        for q in -15592i32..=15592 {
            let passes = if q.abs() <= 16 {
                1
            } else if q.abs() <= 247 {
                2
            } else {
                3
            };
            let s = split(q, passes);
            assert_eq!(s.iter().sum::<i32>(), q);
            assert!(s[passes - 1].abs() <= 16);
        }
    }

    /// A click in a steady tone switches to short blocks around it, and
    /// only there.
    #[test]
    fn transients_switch_to_short_blocks() {
        let rate = 44100;
        let click = rate as i64;
        let tone = |i: usize| (i as f32 * 2.0 * std::f32::consts::PI * 1000.0 / rate as f32).sin() * 0.1;
        let mut with_click: Vec<f32> = (0..2 * rate as usize).map(tone).collect();
        for (k, s) in with_click[click as usize..click as usize + 40].iter_mut().enumerate() {
            *s += 0.8 * (-(k as f32) / 10.0).exp() * if k % 2 == 0 { 1.0 } else { -1.0 };
        }
        let steady: Vec<f32> = (0..2 * rate as usize).map(tone).collect();
        let config = EncoderConfig { channels: 1, ..Default::default() };
        let run = |input: &[f32]| {
            let mut e = Encoder::new(config.clone()).unwrap();
            let mut p = e.encode(input).unwrap();
            p.extend(e.finish().unwrap());
            p
        };
        // The tone's own onset at sample 0 is an attack too: look past it.
        let shorts: Vec<i64> =
            run(&with_click).iter().filter(|p| p.blocksize == 256 && p.granule > 2048).map(|p| p.granule).collect();
        assert!(!shorts.is_empty());
        assert!(shorts.iter().all(|&g| (g - click).abs() < 2048), "short blocks at {shorts:?}");
        assert!(run(&steady).iter().all(|p| p.blocksize == 2048 || p.granule <= 2048));
    }

    #[test]
    fn configurations_are_checked() {
        let bad = |f: fn(&mut EncoderConfig)| {
            let mut c = EncoderConfig::default();
            f(&mut c);
            Encoder::new(c).is_err()
        };
        assert!(bad(|c| c.channels = 0));
        assert!(bad(|c| c.channels = 9));
        assert!(bad(|c| c.sample_rate = 4000));
        assert!(bad(|c| c.quality = 11.0));
        assert!(bad(|c| c.quality = f32::NAN));
        assert!(!bad(|_| {}));
    }
}
