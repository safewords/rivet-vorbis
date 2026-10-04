//! The Ogg encapsulation format (RFC 3533): pages, their CRC, and packets
//! split into lacing segments across pages. Any codec's packets can pass
//! through; [`OggReader`](crate::OggReader) and
//! [`OggWriter`](crate::OggWriter) put Vorbis on top.

use std::collections::{HashMap, VecDeque};
use std::io::{Read, Write};

use crate::error::{Result, invalid};

/// Page header flag: the page starts with the continuation of a packet.
pub const FLAG_CONTINUED: u8 = 0x01;
/// Page header flag: the first page of a logical stream.
pub const FLAG_BOS: u8 = 0x02;
/// Page header flag: the last page of a logical stream.
pub const FLAG_EOS: u8 = 0x04;

/// The CRC of RFC 3533: polynomial 0x04c11db7, initial value 0, no
/// reflection, no final XOR, computed over the page with its CRC field zeroed.
pub fn crc32(data: &[u8]) -> u32 {
    static TABLE: std::sync::OnceLock<[u32; 256]> = std::sync::OnceLock::new();
    let table = TABLE.get_or_init(|| {
        let mut t = [0u32; 256];
        for (i, e) in t.iter_mut().enumerate() {
            let mut r = (i as u32) << 24;
            for _ in 0..8 {
                r = if r & 0x8000_0000 != 0 {
                    (r << 1) ^ 0x04c1_1db7
                } else {
                    r << 1
                };
            }
            *e = r;
        }
        t
    });
    let mut crc = 0u32;
    for &b in data {
        crc = (crc << 8) ^ table[((crc >> 24) as u8 ^ b) as usize];
    }
    crc
}

/// One Ogg page.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Page {
    /// `header_type` flags ([`FLAG_CONTINUED`], [`FLAG_BOS`], [`FLAG_EOS`]).
    pub flags: u8,
    /// Granule position; -1 when no packet ends on the page.
    pub granule: i64,
    /// Bitstream serial number.
    pub serial: u32,
    /// Page sequence number.
    pub sequence: u32,
    /// The segment table (lacing values).
    pub lacing: Vec<u8>,
    /// The page body.
    pub body: Vec<u8>,
}

impl Page {
    /// The page as bytes, with its CRC.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(27 + self.lacing.len() + self.body.len());
        out.extend_from_slice(b"OggS");
        out.push(0);
        out.push(self.flags);
        out.extend_from_slice(&self.granule.to_le_bytes());
        out.extend_from_slice(&self.serial.to_le_bytes());
        out.extend_from_slice(&self.sequence.to_le_bytes());
        out.extend_from_slice(&[0, 0, 0, 0]);
        out.push(self.lacing.len() as u8);
        out.extend_from_slice(&self.lacing);
        out.extend_from_slice(&self.body);
        let crc = crc32(&out);
        out[22..26].copy_from_slice(&crc.to_le_bytes());
        out
    }
}

/// Reads pages from a byte stream, resynchronising on the capture pattern.
pub struct PageReader<R: Read> {
    inner: R,
    buf: Vec<u8>,
    pos: usize,
    eof: bool,
    strict: bool,
    /// Bytes skipped looking for pages, and pages dropped for a bad CRC.
    pub skipped_bytes: u64,
    /// Pages dropped because their CRC did not match.
    pub bad_pages: u64,
}

impl<R: Read> PageReader<R> {
    /// Read pages from `inner`.
    pub fn new(inner: R) -> Self {
        PageReader {
            inner,
            buf: Vec::new(),
            pos: 0,
            eof: false,
            strict: false,
            skipped_bytes: 0,
            bad_pages: 0,
        }
    }

    /// Strict mode: garbage between pages or a CRC mismatch is an error
    /// rather than skipped.
    pub fn set_strict(&mut self, strict: bool) {
        self.strict = strict;
    }

    /// Ensure `n` bytes are buffered from `pos`; false at end of input.
    fn fill(&mut self, n: usize) -> Result<bool> {
        while self.buf.len() - self.pos < n {
            if self.eof {
                return Ok(false);
            }
            if self.pos > 0 && self.pos >= self.buf.len() / 2 {
                self.buf.drain(..self.pos);
                self.pos = 0;
            }
            let mut chunk = [0u8; 65536];
            let got = self.inner.read(&mut chunk)?;
            if got == 0 {
                self.eof = true;
            } else {
                self.buf.extend_from_slice(&chunk[..got]);
            }
        }
        Ok(true)
    }

    /// The next page, or `None` at the end of the input.
    pub fn next_page(&mut self) -> Result<Option<Page>> {
        loop {
            if !self.fill(27)? {
                let left = (self.buf.len() - self.pos) as u64;
                if left > 0 && self.strict {
                    return Err(invalid("trailing bytes after the last Ogg page"));
                }
                self.skipped_bytes += left;
                self.pos = self.buf.len();
                return Ok(None);
            }
            let h = &self.buf[self.pos..];
            if &h[..4] != b"OggS" || h[4] != 0 {
                if self.strict {
                    return Err(invalid("expected an Ogg page"));
                }
                self.pos += 1;
                self.skipped_bytes += 1;
                continue;
            }
            let segments = h[26] as usize;
            if !self.fill(27 + segments)? {
                if self.strict {
                    return Err(invalid("Ogg page truncated in its header"));
                }
                self.skipped_bytes += (self.buf.len() - self.pos) as u64;
                self.pos = self.buf.len();
                return Ok(None);
            }
            let body_len: usize = self.buf[self.pos + 27..self.pos + 27 + segments]
                .iter()
                .map(|&b| b as usize)
                .sum();
            let total = 27 + segments + body_len;
            if !self.fill(total)? {
                if self.strict {
                    return Err(invalid("Ogg page truncated in its body"));
                }
                self.skipped_bytes += (self.buf.len() - self.pos) as u64;
                self.pos = self.buf.len();
                return Ok(None);
            }
            let raw = &self.buf[self.pos..self.pos + total];
            let stored = u32::from_le_bytes([raw[22], raw[23], raw[24], raw[25]]);
            let mut check = raw.to_vec();
            check[22..26].fill(0);
            if crc32(&check) != stored {
                if self.strict {
                    return Err(invalid("Ogg page CRC mismatch"));
                }
                self.bad_pages += 1;
                self.pos += 1;
                self.skipped_bytes += 1;
                continue;
            }
            let page = Page {
                flags: raw[5],
                granule: i64::from_le_bytes(raw[6..14].try_into().expect("8 bytes")),
                serial: u32::from_le_bytes(raw[14..18].try_into().expect("4 bytes")),
                sequence: u32::from_le_bytes(raw[18..22].try_into().expect("4 bytes")),
                lacing: raw[27..27 + segments].to_vec(),
                body: raw[27 + segments..].to_vec(),
            };
            self.pos += total;
            return Ok(Some(page));
        }
    }
}

/// One packet reassembled from pages.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Packet {
    /// The packet's bytes.
    pub data: Vec<u8>,
    /// The logical stream it belongs to.
    pub serial: u32,
    /// The page's granule position, on the last packet completed on a page
    /// (and `None` on the others).
    pub granule: Option<i64>,
    /// The first packet of its logical stream (completed on a BOS page).
    pub bos: bool,
    /// The last packet completed on its stream's EOS page.
    pub eos: bool,
}

/// Reassembles packets from pages, for every logical stream interleaved.
pub struct PacketReader<R: Read> {
    pages: PageReader<R>,
    partial: HashMap<u32, Vec<u8>>,
    /// Next expected page sequence number, per stream.
    sequence: HashMap<u32, u32>,
    queue: VecDeque<Packet>,
    strict: bool,
    /// Packets lost to gaps in the page sequence or broken continuations.
    pub lost_packets: u64,
}

impl<R: Read> PacketReader<R> {
    /// Read packets from `inner`.
    pub fn new(inner: R) -> Self {
        PacketReader {
            pages: PageReader::new(inner),
            partial: HashMap::new(),
            sequence: HashMap::new(),
            queue: VecDeque::new(),
            strict: false,
            lost_packets: 0,
        }
    }

    /// Strict mode: any damage (garbage, CRC mismatch, a gap in the page
    /// sequence, a broken continuation) is an error.
    pub fn set_strict(&mut self, strict: bool) {
        self.strict = strict;
        self.pages.set_strict(strict);
    }

    /// The page reader underneath, for its damage counters.
    pub fn pages(&self) -> &PageReader<R> {
        &self.pages
    }

    /// The next packet of any stream, or `None` at the end of the input.
    /// A packet left incomplete at the end of the input is dropped.
    pub fn next_packet(&mut self) -> Result<Option<Packet>> {
        while self.queue.is_empty() {
            let Some(page) = self.pages.next_page()? else {
                return Ok(None);
            };
            self.take_page(page)?;
        }
        Ok(self.queue.pop_front())
    }

    fn take_page(&mut self, page: Page) -> Result<()> {
        let serial = page.serial;
        if page.flags & FLAG_BOS != 0 {
            self.partial.remove(&serial);
        } else if let Some(&want) = self.sequence.get(&serial)
            && want != page.sequence
        {
            if self.strict {
                return Err(invalid("gap in the Ogg page sequence"));
            }
            if self.partial.remove(&serial).is_some() {
                self.lost_packets += 1;
            }
        }
        self.sequence.insert(serial, page.sequence.wrapping_add(1));
        let mut partial = self.partial.remove(&serial);
        let continued = page.flags & FLAG_CONTINUED != 0;
        // A continuation with nothing to continue: drop the fragment.
        let mut skipping = continued && partial.is_none();
        if !continued && partial.take().is_some() {
            if self.strict {
                return Err(invalid("Ogg packet left unfinished by the next page"));
            }
            self.lost_packets += 1;
        }
        if skipping && self.strict {
            return Err(invalid("Ogg page continues a packet that never began"));
        }
        let mut offset = 0usize;
        let mut completed: Vec<Packet> = Vec::new();
        let mut current = partial.unwrap_or_default();
        for &lace in &page.lacing {
            let seg = &page.body[offset..offset + lace as usize];
            offset += lace as usize;
            if !skipping {
                current.extend_from_slice(seg);
            }
            if lace < 255 {
                if skipping {
                    skipping = false;
                    self.lost_packets += 1;
                } else {
                    completed.push(Packet {
                        data: std::mem::take(&mut current),
                        serial,
                        granule: None,
                        bos: false,
                        eos: false,
                    });
                }
            }
        }
        let bos = page.flags & FLAG_BOS != 0;
        let eos = page.flags & FLAG_EOS != 0;
        if !current.is_empty() || (page.lacing.last() == Some(&255) && !skipping) {
            self.partial.insert(serial, current);
        }
        if let Some(first) = completed.first_mut() {
            first.bos = bos;
        }
        if let Some(last) = completed.last_mut() {
            if page.granule != -1 {
                last.granule = Some(page.granule);
            }
            last.eos = eos;
        } else if eos {
            // An EOS page that completes nothing: report the end with an
            // empty packet carrying the granule, so readers can trim.
            completed.push(Packet {
                data: Vec::new(),
                serial,
                granule: (page.granule != -1).then_some(page.granule),
                bos,
                eos,
            });
        }
        self.queue.extend(completed);
        Ok(())
    }
}

/// Writes packets of one logical stream as pages.
pub struct PacketWriter<W: Write> {
    inner: W,
    serial: u32,
    sequence: u32,
    lacing: Vec<u8>,
    body: Vec<u8>,
    granule: i64,
    /// The page under construction starts mid-packet.
    continued: bool,
    /// No page written yet (the next is BOS).
    first: bool,
    /// Flush a page once its body reaches this many bytes.
    pub page_target: usize,
}

impl<W: Write> PacketWriter<W> {
    /// Write a logical stream with serial number `serial` to `inner`.
    pub fn new(inner: W, serial: u32) -> Self {
        PacketWriter {
            inner,
            serial,
            sequence: 0,
            lacing: Vec::new(),
            body: Vec::new(),
            granule: -1,
            continued: false,
            first: true,
            page_target: 4096,
        }
    }

    fn emit(&mut self, eos: bool) -> Result<()> {
        let mut flags = 0;
        if self.continued {
            flags |= FLAG_CONTINUED;
        }
        if self.first {
            flags |= FLAG_BOS;
        }
        if eos {
            flags |= FLAG_EOS;
        }
        let page = Page {
            flags,
            granule: self.granule,
            serial: self.serial,
            sequence: self.sequence,
            lacing: std::mem::take(&mut self.lacing),
            body: std::mem::take(&mut self.body),
        };
        self.inner.write_all(&page.to_bytes())?;
        self.sequence += 1;
        self.first = false;
        self.granule = -1;
        Ok(())
    }

    /// Append a packet ending at `granule`. `flush` ends the page after it
    /// (as the Vorbis headers require); `eos` ends the stream after it.
    pub fn write_packet(
        &mut self,
        data: &[u8],
        granule: i64,
        flush: bool,
        eos: bool,
    ) -> Result<()> {
        let mut rest = data;
        let mut pushed = false;
        loop {
            if self.lacing.len() == 255 {
                // Full segment table: end the page here; the next one
                // continues this packet if part of it is already out.
                self.emit(false)?;
                self.continued = pushed;
            }
            pushed = true;
            let take = rest.len().min(255);
            self.lacing.push(take as u8);
            self.body.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if take < 255 {
                break;
            }
        }
        self.granule = granule;
        if flush || eos || self.body.len() >= self.page_target {
            self.emit(eos)?;
            self.continued = false;
        }
        Ok(())
    }

    /// End the stream with an EOS page: the pending page if there is one,
    /// otherwise an empty page at `granule`.
    pub fn finish(&mut self, granule: i64) -> Result<()> {
        if self.lacing.is_empty() {
            self.granule = granule;
        }
        self.emit(true)?;
        self.continued = false;
        self.inner.flush()?;
        Ok(())
    }

    /// The underlying writer.
    pub fn into_inner(self) -> W {
        self.inner
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The CRC of RFC 3533 is the plain (unreflected) CRC-32 with
    /// polynomial 0x04c11db7: its check value for "123456789" is 0x89a1897f.
    #[test]
    fn crc_check_value() {
        assert_eq!(crc32(b"123456789"), 0x89a1_897f);
        assert_eq!(crc32(b""), 0);
    }

    #[test]
    fn packets_round_trip_across_pages() {
        let sizes = [
            0usize,
            1,
            254,
            255,
            256,
            510,
            70_000,
            3,
            255 * 255,
            255 * 255 + 1,
            9,
        ];
        let packets: Vec<Vec<u8>> = sizes
            .iter()
            .enumerate()
            .map(|(i, &n)| (0..n).map(|j| (i * 31 + j) as u8).collect())
            .collect();
        let mut w = PacketWriter::new(Vec::new(), 0x1234_5678);
        w.page_target = 1000;
        for (i, p) in packets.iter().enumerate() {
            w.write_packet(p, i as i64 * 100, i == 0, i == packets.len() - 1)
                .unwrap();
        }
        let bytes = w.into_inner();
        let mut r = PacketReader::new(&bytes[..]);
        r.set_strict(true);
        let mut got = Vec::new();
        while let Some(p) = r.next_packet().unwrap() {
            assert_eq!(p.serial, 0x1234_5678);
            got.push(p);
        }
        assert_eq!(got.len(), packets.len());
        for (i, (g, p)) in got.iter().zip(&packets).enumerate() {
            assert_eq!(&g.data, p, "packet {i}");
        }
        assert!(got[0].bos);
        assert_eq!(got[0].granule, Some(0));
        let last = got.last().unwrap();
        assert!(last.eos);
        assert_eq!(last.granule, Some((packets.len() as i64 - 1) * 100));
        // Pages wholly inside one long packet carry granule -1.
        let mut pr = PageReader::new(&bytes[..]);
        let mut minus_one = 0;
        while let Some(page) = pr.next_page().unwrap() {
            if page.granule == -1 {
                minus_one += 1;
                assert!(page.lacing.iter().all(|&l| l == 255));
            }
        }
        assert!(minus_one > 0);
    }

    #[test]
    fn damage_is_skipped_or_refused() {
        let mut w = PacketWriter::new(Vec::new(), 7);
        for i in 0..4 {
            w.write_packet(&[i as u8; 100], i, true, i == 3).unwrap();
        }
        let mut bytes = w.into_inner();
        // Corrupt the second page's body.
        let page_len = 27 + 1 + 100;
        bytes[page_len + 40] ^= 0xff;
        // Leniently: three packets survive, the bad page counted.
        let mut r = PacketReader::new(&bytes[..]);
        let mut n = 0;
        while r.next_packet().unwrap().is_some() {
            n += 1;
        }
        assert_eq!(n, 3);
        assert_eq!(r.pages().bad_pages, 1);
        // Strictly: an error.
        let mut r = PacketReader::new(&bytes[..]);
        r.set_strict(true);
        let mut result = Ok(None);
        for _ in 0..4 {
            result = r.next_packet();
            if result.is_err() {
                break;
            }
        }
        assert!(result.is_err());
    }
}
