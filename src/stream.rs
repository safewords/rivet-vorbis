//! Ogg Vorbis files (specification appendix A): headers from the first
//! pages, audio from the rest, and the granule positions that trim the
//! start and the end.

use std::io::Read;

use crate::decode::Decoder;
use crate::error::{Error, Result, invalid};
use crate::header::{Comments, Identification};
use crate::ogg::{Packet, PacketReader};

/// Decoded PCM from an Ogg Vorbis stream.
#[derive(Clone, Debug, PartialEq)]
pub struct Block {
    /// One vector per channel, in the stream's channel order.
    pub samples: Vec<Vec<f32>>,
    /// Sample rate in Hz.
    pub sample_rate: u32,
    /// Granule position (sample number) of the first sample, when the
    /// stream's granule positions have placed it.
    pub position: Option<i64>,
    /// The serial number of the logical stream (it changes at a chain link).
    pub serial: u32,
}

struct Link {
    serial: u32,
    decoder: Decoder,
    /// Position after the last sample returned.
    position: Option<i64>,
    /// Blocks of the page in progress.
    group: Vec<Vec<Vec<f32>>>,
    ended: bool,
}

/// Decodes the Vorbis stream of an Ogg file: the first Vorbis logical
/// stream, then each one chained after it. Other multiplexed streams are
/// skipped. Leading samples before granule zero and trailing samples past
/// the last page's granule position are trimmed (appendix A.2).
pub struct OggReader<R: Read> {
    packets: PacketReader<R>,
    link: Option<Link>,
    /// Headers of a link being set up: serial and the packets so far.
    pending: Option<(u32, Vec<Vec<u8>>)>,
    strict: bool,
    out: std::collections::VecDeque<Block>,
    done: bool,
}

impl<R: Read> OggReader<R> {
    /// Read from `inner`.
    pub fn new(inner: R) -> Self {
        OggReader { packets: PacketReader::new(inner), link: None, pending: None, strict: false, out: Default::default(), done: false }
    }

    /// Strict mode: Ogg damage, packets the decoder would skip or decode in
    /// part, and granule positions that disagree with the decoded length
    /// are errors (see [`Decoder::set_strict`]).
    pub fn set_strict(&mut self, strict: bool) {
        self.strict = strict;
        self.packets.set_strict(strict);
    }

    /// The identification header of the current link, once read.
    pub fn identification(&self) -> Option<&Identification> {
        self.link.as_ref().map(|l| l.decoder.identification())
    }

    /// The comment header of the current link, once read.
    pub fn comments(&self) -> Option<&Comments> {
        self.link.as_ref().map(|l| l.decoder.comments())
    }

    /// Read until the headers of the first link are in, so that
    /// [`identification`](Self::identification) answers.
    pub fn read_headers(&mut self) -> Result<()> {
        while self.link.is_none() && !self.done {
            self.step()?;
        }
        if self.link.is_none() {
            return Err(invalid("no Vorbis stream in the Ogg data"));
        }
        Ok(())
    }

    /// The next block of trimmed PCM, or `None` at the end.
    pub fn next_block(&mut self) -> Result<Option<Block>> {
        loop {
            if let Some(b) = self.out.pop_front() {
                return Ok(Some(b));
            }
            if self.done {
                return Ok(None);
            }
            self.step()?;
        }
    }

    fn step(&mut self) -> Result<()> {
        let Some(packet) = self.packets.next_packet()? else {
            self.done = true;
            // A stream cut short without an EOS page: release what is held.
            if let Some(link) = self.link.as_mut() {
                let blocks = std::mem::take(&mut link.group);
                let serial = link.serial;
                let rate = link.decoder.sample_rate();
                let mut position = link.position;
                for samples in blocks {
                    let len = samples.first().map_or(0, |c| c.len()) as i64;
                    self.out.push_back(Block { samples, sample_rate: rate, position, serial });
                    position = position.map(|p| p + len);
                }
            }
            return Ok(());
        };
        self.take(packet)
    }

    fn take(&mut self, packet: Packet) -> Result<()> {
        // A new logical stream: Vorbis if its first packet is an
        // identification header, and only if no link is live.
        if packet.bos {
            let live = self.link.as_ref().is_some_and(|l| !l.ended);
            if !live && packet.data.len() >= 7 && packet.data[0] == 1 && &packet.data[1..7] == b"vorbis" {
                self.pending = Some((packet.serial, vec![packet.data]));
                self.link = None;
            }
            return Ok(());
        }
        if let Some((serial, headers)) = self.pending.as_mut() {
            if packet.serial != *serial {
                return Ok(());
            }
            headers.push(packet.data);
            if headers.len() == 3 {
                let (serial, headers) = self.pending.take().expect("pending headers");
                let mut decoder = Decoder::new(&headers[0], &headers[1], &headers[2])?;
                decoder.set_strict(self.strict);
                self.link = Some(Link { serial, decoder, position: None, group: Vec::new(), ended: false });
            }
            return Ok(());
        }
        let strict = self.strict;
        let Some(link) = self.link.as_mut() else { return Ok(()) };
        if packet.serial != link.serial || link.ended {
            return Ok(());
        }
        if !packet.data.is_empty() {
            let samples = link.decoder.decode(&packet.data)?;
            if samples.first().is_some_and(|c| !c.is_empty()) {
                link.group.push(samples);
            }
        }
        if packet.granule.is_none() && !packet.eos {
            return Ok(());
        }
        // The end of a page: place (and trim) its samples.
        let mut blocks = std::mem::take(&mut link.group);
        let total: i64 = blocks.iter().map(|b| b[0].len() as i64).sum();
        let rate = link.decoder.sample_rate();
        let serial = link.serial;
        let start = match (packet.granule, link.position) {
            (Some(g), None) if !packet.eos => {
                let start = g - total;
                if start < 0 {
                    trim_front(&mut blocks, (-start) as usize);
                    0
                } else {
                    start
                }
            }
            (Some(g), pos) if packet.eos => {
                let start = pos.unwrap_or(0);
                let keep = (g - start).max(0);
                if keep < total {
                    trim_back(&mut blocks, (total - keep) as usize);
                } else if keep > total && strict && pos.is_some() {
                    return Err(invalid("last page's granule position is past the decoded end"));
                }
                start
            }
            (Some(g), Some(pos)) => {
                if strict && g != pos + total {
                    return Err(invalid(format!("granule position {g} where the decoded length gives {}", pos + total)));
                }
                pos
            }
            (None, pos) => pos.unwrap_or(0),
            (Some(_), None) => unreachable!(),
        };
        let mut position = start;
        for samples in blocks {
            let len = samples[0].len() as i64;
            if len > 0 {
                self.out.push_back(Block { samples, sample_rate: rate, position: Some(position), serial });
            }
            position += len;
        }
        link.position = Some(position);
        if packet.eos {
            link.ended = true;
        }
        Ok(())
    }
}

fn trim_front(blocks: &mut Vec<Vec<Vec<f32>>>, mut n: usize) {
    while n > 0 && !blocks.is_empty() {
        let len = blocks[0][0].len();
        if len <= n {
            blocks.remove(0);
            n -= len;
        } else {
            for c in blocks[0].iter_mut() {
                c.drain(..n);
            }
            n = 0;
        }
    }
}

fn trim_back(blocks: &mut Vec<Vec<Vec<f32>>>, mut n: usize) {
    while n > 0 && !blocks.is_empty() {
        let last = blocks.len() - 1;
        let len = blocks[last][0].len();
        if len <= n {
            blocks.pop();
            n -= len;
        } else {
            for c in blocks[last].iter_mut() {
                c.truncate(len - n);
            }
            n = 0;
        }
    }
}

/// A whole Ogg Vorbis file, decoded.
#[derive(Clone, Debug, PartialEq)]
pub struct Decoded {
    /// Sample rate in Hz.
    pub sample_rate: u32,
    /// One vector per channel, in Vorbis channel order.
    pub samples: Vec<Vec<f32>>,
    /// The comment header of the first link.
    pub comments: Comments,
    /// The identification header of the first link.
    pub identification: Identification,
}

/// Decode a whole Ogg Vorbis file. Chained links are concatenated; a link
/// with a different channel count or sample rate is refused.
pub fn decode_ogg(bytes: &[u8]) -> Result<Decoded> {
    decode_ogg_with(bytes, false)
}

/// As [`decode_ogg`], in strict mode (see [`OggReader::set_strict`]).
pub fn decode_ogg_strict(bytes: &[u8]) -> Result<Decoded> {
    decode_ogg_with(bytes, true)
}

fn decode_ogg_with(bytes: &[u8], strict: bool) -> Result<Decoded> {
    let mut reader = OggReader::new(bytes);
    reader.set_strict(strict);
    reader.read_headers()?;
    let identification = reader.identification().cloned().expect("headers read");
    let comments = reader.comments().cloned().unwrap_or_default();
    let mut samples = vec![Vec::new(); identification.channels as usize];
    while let Some(block) = reader.next_block()? {
        if block.samples.len() != samples.len() || block.sample_rate != identification.sample_rate {
            return Err(Error::Unsupported("chained links with different channel counts or sample rates".into()));
        }
        for (s, b) in samples.iter_mut().zip(block.samples) {
            s.extend_from_slice(&b);
        }
    }
    Ok(Decoded { sample_rate: identification.sample_rate, samples, comments, identification })
}
