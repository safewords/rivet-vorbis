//! The three header packets (specification 4.2 and section 5): the
//! identification header, the comment header and the setup header, read
//! and written; and the Xiph lacing that packs them into one buffer for
//! Matroska `CodecPrivate`.

use crate::bits::{BitReader, BitWriter, ilog};
use crate::codebook::Codebook;
use crate::error::{Error, Result, invalid};
use crate::floor0::Floor0;
use crate::floor1::Floor1;
use crate::residue::Residue;

/// Check the common header (4.2.1): the packet type byte and "vorbis".
fn check_common(packet: &[u8], packet_type: u8, name: &str) -> Result<()> {
    if packet.len() < 7 || packet[0] != packet_type || &packet[1..7] != b"vorbis" {
        return Err(invalid(format!("not a Vorbis {name} header")));
    }
    Ok(())
}

fn write_common(w: &mut BitWriter, packet_type: u8) {
    w.write(packet_type as u32, 8);
    for &b in b"vorbis" {
        w.write(b as u32, 8);
    }
}

/// The identification header (4.2.2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Identification {
    /// `audio_channels`, at least 1.
    pub channels: u8,
    /// `audio_sample_rate`, in Hz.
    pub sample_rate: u32,
    /// `bitrate_maximum`, bits per second; meaningful only when positive.
    pub bitrate_maximum: i32,
    /// `bitrate_nominal`.
    pub bitrate_nominal: i32,
    /// `bitrate_minimum`.
    pub bitrate_minimum: i32,
    /// `blocksize_0` and `blocksize_1`, in samples: powers of two from 64
    /// to 8192, the first not above the second.
    pub blocksize: [u16; 2],
}

impl Identification {
    /// Decode the first header packet.
    pub fn read(packet: &[u8]) -> Result<Self> {
        check_common(packet, 1, "identification")?;
        let mut r = BitReader::new(&packet[7..]);
        let eop = |_| invalid("identification header is too short");
        let version = r.read(32).map_err(eop)?;
        let channels = r.read(8).map_err(eop)? as u8;
        let sample_rate = r.read(32).map_err(eop)?;
        let bitrate_maximum = r.read(32).map_err(eop)? as i32;
        let bitrate_nominal = r.read(32).map_err(eop)? as i32;
        let bitrate_minimum = r.read(32).map_err(eop)? as i32;
        let b0 = r.read(4).map_err(eop)?;
        let b1 = r.read(4).map_err(eop)?;
        let framing = r.read_flag().map_err(eop)?;
        if version != 0 {
            return Err(Error::Unsupported(format!("Vorbis version {version}")));
        }
        if channels == 0 || sample_rate == 0 {
            return Err(invalid("zero channels or zero sample rate"));
        }
        if !(6..=13).contains(&b0) || !(6..=13).contains(&b1) || b0 > b1 {
            return Err(invalid("block sizes outside 64..=8192 or out of order"));
        }
        if !framing {
            return Err(invalid("identification header framing bit unset"));
        }
        Ok(Identification { channels, sample_rate, bitrate_maximum, bitrate_nominal, bitrate_minimum, blocksize: [1 << b0, 1 << b1] })
    }

    /// Encode the first header packet (30 bytes).
    pub fn write(&self) -> Vec<u8> {
        let mut w = BitWriter::new();
        write_common(&mut w, 1);
        w.write(0, 32);
        w.write(self.channels as u32, 8);
        w.write(self.sample_rate, 32);
        w.write(self.bitrate_maximum as u32, 32);
        w.write(self.bitrate_nominal as u32, 32);
        w.write(self.bitrate_minimum as u32, 32);
        w.write(self.blocksize[0].trailing_zeros(), 4);
        w.write(self.blocksize[1].trailing_zeros(), 4);
        w.write_flag(true);
        w.into_bytes()
    }
}

/// The comment header (section 5): a vendor string and `NAME=value`
/// comments, eight-bit clean (decoded as UTF-8, lossily).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Comments {
    /// The vendor string.
    pub vendor: String,
    /// The user comments, each normally `NAME=value`.
    pub comments: Vec<String>,
}

impl Comments {
    /// Decode the second header packet. The packet type and signature
    /// must be right; an end-of-packet inside the body (a non-fatal
    /// condition by 4.2) is an error here, which a lenient caller may
    /// ignore with [`read_lenient`](Self::read_lenient).
    pub fn read(packet: &[u8]) -> Result<Self> {
        check_common(packet, 3, "comment")?;
        let body = &packet[7..];
        let mut pos = 0usize;
        let u32_at = |pos: &mut usize| -> Result<u32> {
            let b = body.get(*pos..*pos + 4).ok_or_else(|| invalid("comment header ends early"))?;
            *pos += 4;
            Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        };
        let vendor_len = u32_at(&mut pos)? as usize;
        let vendor = body.get(pos..pos.saturating_add(vendor_len)).ok_or_else(|| invalid("comment header ends early"))?;
        pos += vendor_len;
        let vendor = String::from_utf8_lossy(vendor).into_owned();
        let count = u32_at(&mut pos)?;
        let mut comments = Vec::new();
        for _ in 0..count {
            let len = u32_at(&mut pos)? as usize;
            let c = body.get(pos..pos.saturating_add(len)).ok_or_else(|| invalid("comment header ends early"))?;
            pos += len;
            comments.push(String::from_utf8_lossy(c).into_owned());
        }
        match body.get(pos) {
            Some(b) if b & 1 == 1 => Ok(Comments { vendor, comments }),
            _ => Err(invalid("comment header framing bit unset or missing")),
        }
    }

    /// As [`read`](Self::read), but a malformed body yields whatever was
    /// read before it went wrong (4.2: end-of-packet in the comment header
    /// is non-fatal). Only a packet that is not a comment header fails.
    pub fn read_lenient(packet: &[u8]) -> Result<Self> {
        check_common(packet, 3, "comment")?;
        Ok(Self::read(packet).unwrap_or_default())
    }

    /// Encode the second header packet.
    pub fn write(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.push(3);
        out.extend_from_slice(b"vorbis");
        out.extend_from_slice(&(self.vendor.len() as u32).to_le_bytes());
        out.extend_from_slice(self.vendor.as_bytes());
        out.extend_from_slice(&(self.comments.len() as u32).to_le_bytes());
        for c in &self.comments {
            out.extend_from_slice(&(c.len() as u32).to_le_bytes());
            out.extend_from_slice(c.as_bytes());
        }
        out.push(1);
        out
    }

    /// The values of every comment named `name` (field names compare
    /// case-insensitively, as section 5.2.2 has them).
    pub fn get<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a str> + 'a {
        self.comments.iter().filter_map(move |c| {
            let (k, v) = c.split_once('=')?;
            k.eq_ignore_ascii_case(name).then_some(v)
        })
    }
}

/// A floor configuration.
#[derive(Clone, Debug, PartialEq)]
pub enum Floor {
    /// Type 0: LSP.
    Zero(Floor0),
    /// Type 1: piecewise linear.
    One(Floor1),
}

/// A mapping configuration (type 0, 4.2.4).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mapping {
    /// Coupling steps: (magnitude channel, angle channel).
    pub coupling: Vec<(u8, u8)>,
    /// `vorbis_mapping_mux`: the submap of each channel.
    pub mux: Vec<u8>,
    /// Floor number of each submap.
    pub submap_floor: Vec<u8>,
    /// Residue number of each submap.
    pub submap_residue: Vec<u8>,
}

/// A mode configuration (4.2.4).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mode {
    /// Long block (`blocksize_1`) when set.
    pub blockflag: bool,
    /// `vorbis_mode_mapping`.
    pub mapping: u8,
}

/// The setup header (4.2.4).
#[derive(Clone, Debug, PartialEq)]
pub struct Setup {
    /// The codebooks, 1 to 256.
    pub codebooks: Vec<Codebook>,
    /// The floors, 1 to 64.
    pub floors: Vec<Floor>,
    /// The residues, 1 to 64.
    pub residues: Vec<Residue>,
    /// The mappings, 1 to 64.
    pub mappings: Vec<Mapping>,
    /// The modes, 1 to 64.
    pub modes: Vec<Mode>,
}

impl Setup {
    /// Decode the third header packet for a stream of `channels` channels.
    pub fn read(packet: &[u8], channels: u8) -> Result<Self> {
        check_common(packet, 5, "setup")?;
        let mut r = BitReader::new(&packet[7..]);
        let eop = |_| invalid("setup header ends early");
        let count = r.read(8).map_err(eop)? + 1;
        let mut codebooks = Vec::with_capacity(count as usize);
        for _ in 0..count {
            codebooks.push(Codebook::read(&mut r)?);
        }
        let times = r.read(6).map_err(eop)? + 1;
        for _ in 0..times {
            if r.read(16).map_err(eop)? != 0 {
                return Err(invalid("nonzero time-domain transform type"));
            }
        }
        let count = r.read(6).map_err(eop)? + 1;
        let mut floors = Vec::new();
        for _ in 0..count {
            floors.push(match r.read(16).map_err(eop)? {
                0 => Floor::Zero(Floor0::read(&mut r, &codebooks)?),
                1 => Floor::One(Floor1::read(&mut r, &codebooks)?),
                t => return Err(invalid(format!("floor type {t}"))),
            });
        }
        let count = r.read(6).map_err(eop)? + 1;
        let mut residues = Vec::new();
        for _ in 0..count {
            let t = r.read(16).map_err(eop)?;
            if t > 2 {
                return Err(invalid(format!("residue type {t}")));
            }
            residues.push(Residue::read(&mut r, t as u16, &codebooks)?);
        }
        let count = r.read(6).map_err(eop)? + 1;
        let mut mappings = Vec::new();
        let chbits = ilog(channels as i64 - 1);
        for _ in 0..count {
            if r.read(16).map_err(eop)? != 0 {
                return Err(invalid("mapping type other than 0"));
            }
            let submaps = if r.read_flag().map_err(eop)? { r.read(4).map_err(eop)? + 1 } else { 1 };
            let mut coupling = Vec::new();
            if r.read_flag().map_err(eop)? {
                let steps = r.read(8).map_err(eop)? + 1;
                for _ in 0..steps {
                    let m = r.read(chbits).map_err(eop)?;
                    let a = r.read(chbits).map_err(eop)?;
                    if m == a || m >= channels as u32 || a >= channels as u32 {
                        return Err(invalid("invalid coupling step"));
                    }
                    coupling.push((m as u8, a as u8));
                }
            }
            if r.read(2).map_err(eop)? != 0 {
                return Err(invalid("mapping reserved field nonzero"));
            }
            let mut mux = vec![0u8; channels as usize];
            if submaps > 1 {
                for m in mux.iter_mut() {
                    *m = r.read(4).map_err(eop)? as u8;
                    if *m as u32 >= submaps {
                        return Err(invalid("mapping mux names a missing submap"));
                    }
                }
            }
            let mut submap_floor = Vec::new();
            let mut submap_residue = Vec::new();
            for _ in 0..submaps {
                r.read(8).map_err(eop)?;
                let f = r.read(8).map_err(eop)?;
                let res = r.read(8).map_err(eop)?;
                if f as usize >= floors.len() || res as usize >= residues.len() {
                    return Err(invalid("submap names a missing floor or residue"));
                }
                submap_floor.push(f as u8);
                submap_residue.push(res as u8);
            }
            mappings.push(Mapping { coupling, mux, submap_floor, submap_residue });
        }
        let count = r.read(6).map_err(eop)? + 1;
        let mut modes = Vec::new();
        for _ in 0..count {
            let blockflag = r.read_flag().map_err(eop)?;
            let windowtype = r.read(16).map_err(eop)?;
            let transformtype = r.read(16).map_err(eop)?;
            let mapping = r.read(8).map_err(eop)?;
            if windowtype != 0 || transformtype != 0 || mapping as usize >= mappings.len() {
                return Err(invalid("invalid mode"));
            }
            modes.push(Mode { blockflag, mapping: mapping as u8 });
        }
        if !r.read_flag().map_err(eop)? {
            return Err(invalid("setup header framing bit unset"));
        }
        Ok(Setup { codebooks, floors, residues, mappings, modes })
    }

    /// Encode the third header packet for a stream of `channels` channels.
    pub fn write(&self, channels: u8) -> Vec<u8> {
        let mut w = BitWriter::new();
        write_common(&mut w, 5);
        w.write(self.codebooks.len() as u32 - 1, 8);
        for b in &self.codebooks {
            b.write(&mut w);
        }
        w.write(0, 6);
        w.write(0, 16);
        w.write(self.floors.len() as u32 - 1, 6);
        for f in &self.floors {
            match f {
                Floor::Zero(f) => {
                    w.write(0, 16);
                    f.write(&mut w);
                }
                Floor::One(f) => {
                    w.write(1, 16);
                    f.write(&mut w);
                }
            }
        }
        w.write(self.residues.len() as u32 - 1, 6);
        for r in &self.residues {
            w.write(r.residue_type as u32, 16);
            r.write(&mut w);
        }
        w.write(self.mappings.len() as u32 - 1, 6);
        let chbits = ilog(channels as i64 - 1);
        for m in &self.mappings {
            w.write(0, 16);
            let submaps = m.submap_floor.len() as u32;
            w.write_flag(submaps > 1);
            if submaps > 1 {
                w.write(submaps - 1, 4);
            }
            w.write_flag(!m.coupling.is_empty());
            if !m.coupling.is_empty() {
                w.write(m.coupling.len() as u32 - 1, 8);
                for &(mag, ang) in &m.coupling {
                    w.write(mag as u32, chbits);
                    w.write(ang as u32, chbits);
                }
            }
            w.write(0, 2);
            if submaps > 1 {
                for &x in &m.mux {
                    w.write(x as u32, 4);
                }
            }
            for (f, r) in m.submap_floor.iter().zip(&m.submap_residue) {
                w.write(0, 8);
                w.write(*f as u32, 8);
                w.write(*r as u32, 8);
            }
        }
        w.write(self.modes.len() as u32 - 1, 6);
        for m in &self.modes {
            w.write_flag(m.blockflag);
            w.write(0, 16);
            w.write(0, 16);
            w.write(m.mapping as u32, 8);
        }
        w.write_flag(true);
        w.into_bytes()
    }
}

/// Split Xiph-laced headers (Matroska `CodecPrivate`): a count less one,
/// the lacing of every packet but the last, then the packets. Vorbis has
/// exactly three.
pub fn split_xiph_lacing(bytes: &[u8]) -> Result<[&[u8]; 3]> {
    let (&count, mut rest) = bytes.split_first().ok_or_else(|| invalid("empty Xiph-laced header buffer"))?;
    if count != 2 {
        return Err(invalid(format!("Xiph lacing holds {} packets, not 3", count as usize + 1)));
    }
    let mut lens = [0usize; 2];
    for len in lens.iter_mut() {
        loop {
            let (&b, r) = rest.split_first().ok_or_else(|| invalid("Xiph lacing ends inside a length"))?;
            rest = r;
            *len += b as usize;
            if b != 255 {
                break;
            }
        }
    }
    if lens[0] + lens[1] > rest.len() {
        return Err(invalid("Xiph lacing lengths exceed the buffer"));
    }
    let (a, rest) = rest.split_at(lens[0]);
    let (b, c) = rest.split_at(lens[1]);
    Ok([a, b, c])
}

/// Pack the three headers with Xiph lacing (Matroska `CodecPrivate`).
pub fn xiph_lacing(headers: [&[u8]; 3]) -> Vec<u8> {
    let mut out = vec![2u8];
    for h in &headers[..2] {
        let mut n = h.len();
        while n >= 255 {
            out.push(255);
            n -= 255;
        }
        out.push(n as u8);
    }
    for h in headers {
        out.extend_from_slice(h);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identification_round_trips_and_is_30_bytes() {
        let id = Identification { channels: 2, sample_rate: 44100, bitrate_maximum: 0, bitrate_nominal: 128000, bitrate_minimum: 0, blocksize: [256, 2048] };
        let bytes = id.write();
        assert_eq!(bytes.len(), 30);
        assert_eq!(Identification::read(&bytes).unwrap(), id);
        // Framing bit cleared.
        let mut bad = bytes.clone();
        bad[29] = 0;
        assert!(Identification::read(&bad).is_err());
        // blocksize_0 > blocksize_1.
        let mut bad = bytes.clone();
        assert_eq!(bytes[28], 0xb8);
        bad[28] = 0x8b; // b0 = 11, b1 = 8
        assert!(Identification::read(&bad).is_err());
    }

    #[test]
    fn comments_round_trip() {
        let c = Comments { vendor: "rivet-vorbis".into(), comments: vec!["TITLE=Gr\u{fc}n".into(), "ARTIST=a".into(), "artist=b".into()] };
        let bytes = c.write();
        let back = Comments::read(&bytes).unwrap();
        assert_eq!(back, c);
        assert_eq!(back.get("Artist").collect::<Vec<_>>(), vec!["a", "b"]);
        // Truncated: an error strictly, empty leniently.
        assert!(Comments::read(&bytes[..bytes.len() - 3]).is_err());
        assert_eq!(Comments::read_lenient(&bytes[..bytes.len() - 3]).unwrap(), Comments::default());
        assert!(Comments::read_lenient(&[5, b'v']).is_err());
    }

    #[test]
    fn xiph_lacing_round_trips() {
        let a = vec![1u8; 30];
        let b = vec![2u8; 300];
        let c = vec![3u8; 5];
        let laced = xiph_lacing([&a, &b, &c]);
        assert_eq!(&laced[..4], &[2, 30, 255, 45]);
        let [x, y, z] = split_xiph_lacing(&laced).unwrap();
        assert_eq!((x, y, z), (&a[..], &b[..], &c[..]));
        assert!(split_xiph_lacing(&[1, 5, 5]).is_err());
        assert!(split_xiph_lacing(&[2, 30, 19, 0, 0]).is_err());
    }
}
