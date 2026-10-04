//! The audio packet decoder (specification 4.3): one Vorbis packet in,
//! PCM out, independent of the container.

use crate::bits::{BitReader, ilog};
use crate::error::{Result, invalid};
use crate::floor0::{Floor0Error, Floor0Frame};
use crate::header::{Comments, Floor, Identification, Setup, split_xiph_lacing};
use crate::mdct::Mdct;
use crate::window::Windows;

/// Inverse square polar coupling (4.3.5) of one magnitude/angle pair.
#[inline]
pub fn decouple(m: f32, a: f32) -> (f32, f32) {
    if m > 0.0 {
        if a > 0.0 { (m, m - a) } else { (m + a, m) }
    } else if a > 0.0 {
        (m, m + a)
    } else {
        (m - a, m)
    }
}

enum FloorData {
    Unused,
    Zero(Floor0Frame),
    One(Vec<i32>),
}

/// A Vorbis I decoder for one logical stream.
///
/// Feed it the stream's audio packets in order; each returns the PCM it
/// finishes, one vector per channel in the stream's channel order (4.3.9:
/// for 5.1, FL FC FR RL RR LFE). The first packet returns no samples: it
/// primes the overlap. Samples are `f32` at full scale ±1.0, unclipped.
///
/// The decoder knows nothing of granule positions: trimming the start of a
/// stream and the end of its last page is the container's job, which
/// [`OggReader`](crate::OggReader) does for Ogg.
pub struct Decoder {
    ident: Identification,
    comments: Comments,
    setup: Setup,
    strict: bool,
    windows: Windows,
    mdct: [Mdct; 2],
    /// Bark maps of each floor 0, for each block size.
    floor0_maps: Vec<Option<[Vec<i32>; 2]>>,
    mode_bits: u32,
    /// The previous block, windowed: whether it was long, and its samples.
    prev: Option<(bool, Vec<Vec<f32>>)>,
}

impl Decoder {
    /// A decoder from the three header packets.
    pub fn new(identification: &[u8], comment: &[u8], setup: &[u8]) -> Result<Self> {
        let ident = Identification::read(identification)?;
        let comments = Comments::read_lenient(comment)?;
        let setup = Setup::read(setup, ident.channels)?;
        let bs = [ident.blocksize[0] as usize, ident.blocksize[1] as usize];
        let floor0_maps = setup
            .floors
            .iter()
            .map(|f| match f {
                Floor::Zero(f0) => Some([f0.bark_map(bs[0] / 2), f0.bark_map(bs[1] / 2)]),
                Floor::One(_) => None,
            })
            .collect();
        let mode_bits = ilog(setup.modes.len() as i64 - 1);
        Ok(Decoder {
            windows: Windows::new(bs),
            mdct: [Mdct::new(bs[0]), Mdct::new(bs[1])],
            ident,
            comments,
            setup,
            strict: false,
            floor0_maps,
            mode_bits,
            prev: None,
        })
    }

    /// A decoder from the three headers packed with Xiph lacing, as
    /// Matroska's `CodecPrivate` carries them.
    pub fn from_xiph_lacing(codec_private: &[u8]) -> Result<Self> {
        let [a, b, c] = split_xiph_lacing(codec_private)?;
        Decoder::new(a, b, c)
    }

    /// Strict mode: a packet that would otherwise be ignored or decoded in
    /// part (a header-type packet among the audio, a truncated packet, an
    /// undecodable floor) is an error instead. For checking an encoder; a
    /// conforming stream may truncate packets on purpose (2.1.8), so
    /// players should leave this off.
    pub fn set_strict(&mut self, strict: bool) {
        self.strict = strict;
    }

    /// The identification header.
    pub fn identification(&self) -> &Identification {
        &self.ident
    }

    /// The comment header.
    pub fn comments(&self) -> &Comments {
        &self.comments
    }

    /// The setup header.
    pub fn setup(&self) -> &Setup {
        &self.setup
    }

    /// Channels.
    pub fn channels(&self) -> usize {
        self.ident.channels as usize
    }

    /// Sample rate in Hz.
    pub fn sample_rate(&self) -> u32 {
        self.ident.sample_rate
    }

    /// Forget the previous block (after a seek): the next packet primes
    /// the overlap again and returns no samples.
    pub fn reset(&mut self) {
        self.prev = None;
    }

    /// The block size of an audio packet (from its mode), without decoding
    /// it; `None` for a packet that is not a valid audio packet.
    pub fn packet_blocksize(&self, packet: &[u8]) -> Option<usize> {
        let mut r = BitReader::new(packet);
        if r.read(1).ok()? != 0 {
            return None;
        }
        let mode = self
            .setup
            .modes
            .get(r.read(self.mode_bits).ok()? as usize)?;
        Some(self.ident.blocksize[mode.blockflag as usize] as usize)
    }

    fn discard(&self, why: &str) -> Result<Vec<Vec<f32>>> {
        if self.strict {
            Err(invalid(why.to_string()))
        } else {
            Ok(vec![Vec::new(); self.channels()])
        }
    }

    /// Decode one audio packet; returns the finished samples, one vector
    /// per channel (empty for the first packet, and for a packet that is
    /// ignored outside strict mode).
    pub fn decode(&mut self, packet: &[u8]) -> Result<Vec<Vec<f32>>> {
        let ch = self.channels();
        let mut r = BitReader::new(packet);
        // 4.3.1: packet type, mode, window flags. End-of-packet here
        // discards the packet.
        let Ok(packet_type) = r.read(1) else {
            return self.discard("empty or truncated audio packet");
        };
        if packet_type != 0 {
            return self.discard("non-audio packet among audio packets");
        }
        let Ok(mode_number) = r.read(self.mode_bits) else {
            return self.discard("audio packet truncated in its mode");
        };
        let Some(mode) = self.setup.modes.get(mode_number as usize).cloned() else {
            return self.discard("audio packet names a missing mode");
        };
        let long = mode.blockflag;
        let n = self.ident.blocksize[long as usize] as usize;
        let (prev_flag, next_flag) = if long {
            match (r.read_flag(), r.read_flag()) {
                (Ok(p), Ok(nx)) => (p, nx),
                _ => return self.discard("audio packet truncated in its window flags"),
            }
        } else {
            (false, false)
        };
        let setup = &self.setup;
        let mapping = &setup.mappings[mode.mapping as usize];
        let books = &setup.codebooks;
        let half = n / 2;

        // 4.3.2: floors, in channel order.
        let mut floors: Vec<FloorData> = Vec::with_capacity(ch);
        let mut lost = false;
        for i in 0..ch {
            let floor_number = mapping.submap_floor[mapping.mux[i] as usize] as usize;
            let data = match &setup.floors[floor_number] {
                Floor::Zero(f) => match f.decode(&mut r, books) {
                    Ok(Some(frame)) => FloorData::Zero(frame),
                    Ok(None) => FloorData::Unused,
                    Err(Floor0Error::EndOfPacket) => {
                        lost = true;
                        FloorData::Unused
                    }
                    Err(Floor0Error::Undecodable(why)) => {
                        if self.strict {
                            return Err(invalid(why));
                        }
                        lost = true;
                        FloorData::Unused
                    }
                },
                Floor::One(f) => match f.decode(&mut r, books) {
                    Ok(Some(y)) => FloorData::One(y),
                    Ok(None) => FloorData::Unused,
                    Err(_) => {
                        lost = true;
                        FloorData::Unused
                    }
                },
            };
            floors.push(data);
        }
        if lost && self.strict {
            return Err(invalid("audio packet ends inside a floor"));
        }

        let mut spectra = vec![vec![0f32; half]; ch];
        if !lost {
            // 4.3.3: nonzero vector propagate.
            let mut no_residue: Vec<bool> = floors
                .iter()
                .map(|f| matches!(f, FloorData::Unused))
                .collect();
            for &(m, a) in &mapping.coupling {
                if !no_residue[m as usize] || !no_residue[a as usize] {
                    no_residue[m as usize] = false;
                    no_residue[a as usize] = false;
                }
            }
            // 4.3.4: residues, by submap.
            let mut truncated = false;
            for (s, &residue_number) in mapping.submap_residue.iter().enumerate() {
                let members: Vec<usize> =
                    (0..ch).filter(|&j| mapping.mux[j] as usize == s).collect();
                if members.is_empty() {
                    continue;
                }
                let flags: Vec<bool> = members.iter().map(|&j| no_residue[j]).collect();
                let mut vectors = vec![vec![0f32; half]; members.len()];
                if !truncated
                    && setup.residues[residue_number as usize]
                        .decode(books, &mut r, &flags, &mut vectors)
                        .is_err()
                {
                    truncated = true;
                }
                for (v, &j) in vectors.into_iter().zip(&members) {
                    spectra[j] = v;
                }
            }
            if truncated && self.strict {
                return Err(invalid("audio packet ends inside a residue"));
            }
            // 4.3.5: inverse coupling, last step first.
            for &(m, a) in mapping.coupling.iter().rev() {
                let (m, a) = (m as usize, a as usize);
                let (lo, hi) = (m.min(a), m.max(a));
                let (left, right) = spectra.split_at_mut(hi);
                let (vm, va) = if m < a {
                    (&mut left[lo], &mut right[0])
                } else {
                    (&mut right[0], &mut left[lo])
                };
                for (x, y) in vm.iter_mut().zip(va.iter_mut()) {
                    let (nm, na) = decouple(*x, *y);
                    *x = nm;
                    *y = na;
                }
            }
            // 4.3.6: floor curves times residues.
            for (i, floor) in floors.iter().enumerate() {
                let floor_number = mapping.submap_floor[mapping.mux[i] as usize] as usize;
                match (floor, &setup.floors[floor_number]) {
                    (FloorData::Unused, _) => spectra[i].fill(0.0),
                    (FloorData::One(y), Floor::One(f)) => {
                        let curve = f.synthesize(y, half);
                        for (s, c) in spectra[i].iter_mut().zip(&curve) {
                            *s *= *c;
                        }
                    }
                    (FloorData::Zero(frame), Floor::Zero(f)) => {
                        let map = &self.floor0_maps[floor_number]
                            .as_ref()
                            .expect("floor 0 map")[long as usize];
                        let mut curve = vec![0f32; half];
                        f.synthesize(frame, map, &mut curve);
                        for (s, c) in spectra[i].iter_mut().zip(&curve) {
                            *s *= *c;
                        }
                    }
                    _ => unreachable!("floor data matches its floor type"),
                }
            }
        }

        // 4.3.7: inverse MDCT and window.
        let mut blocks = Vec::with_capacity(ch);
        for spectrum in &spectra {
            let mut y = vec![0f32; n];
            if spectrum.iter().any(|&v| v != 0.0) {
                self.mdct[long as usize].inverse(spectrum, &mut y);
                self.windows.apply(&mut y, long, prev_flag, next_flag);
            }
            blocks.push(y);
        }

        // 4.3.8: overlap-add with the previous block's right half.
        let out = match &self.prev {
            None => vec![Vec::new(); ch],
            Some((prev_long, prev)) => {
                let np = self.ident.blocksize[*prev_long as usize] as usize;
                let len = np / 4 + n / 4;
                let shift = n as i64 / 4 - np as i64 / 4;
                prev.iter()
                    .zip(&blocks)
                    .map(|(p, c)| {
                        (0..len)
                            .map(|i| {
                                let a = p.get(np / 2 + i).copied().unwrap_or(0.0);
                                let j = i as i64 + shift;
                                let b = if j >= 0 {
                                    c.get(j as usize).copied().unwrap_or(0.0)
                                } else {
                                    0.0
                                };
                                a + b
                            })
                            .collect()
                    })
                    .collect()
            }
        };
        self.prev = Some((long, blocks));
        Ok(out)
    }

    /// As [`decode`](Self::decode), interleaved.
    pub fn decode_interleaved(&mut self, packet: &[u8]) -> Result<Vec<f32>> {
        Ok(interleave(&self.decode(packet)?))
    }
}

/// Interleave per-channel vectors of equal length.
pub fn interleave(planar: &[Vec<f32>]) -> Vec<f32> {
    let len = planar.first().map_or(0, |c| c.len());
    let mut out = Vec::with_capacity(len * planar.len());
    for i in 0..len {
        for c in planar {
            out.push(c[i]);
        }
    }
    out
}

impl std::fmt::Debug for Decoder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Decoder")
            .field("ident", &self.ident)
            .field("strict", &self.strict)
            .finish_non_exhaustive()
    }
}
