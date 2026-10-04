//! Residues (specification section 8): the fine spectral detail,
//! partitioned, classified, and coded in up to eight VQ passes. Types 0
//! and 1 differ in how a partition's scalars are interleaved; type 2 codes
//! all channels as one interleaved vector.

use crate::bits::{BitReader, BitWriter, EndOfPacket};
use crate::codebook::Codebook;
use crate::error::{Result, invalid};

/// A residue configuration (8.6.1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Residue {
    /// 0, 1 or 2.
    pub residue_type: u16,
    /// `residue_begin`.
    pub begin: u32,
    /// `residue_end`.
    pub end: u32,
    /// `residue_partition_size`.
    pub partition_size: u32,
    /// `residue_classifications`, 1 to 64.
    pub classifications: u8,
    /// `residue_classbook`.
    pub classbook: u8,
    /// `residue_cascade`: per classification, which passes code values.
    pub cascade: Vec<u8>,
    /// `residue_books`: per classification, per pass, a book or -1.
    pub books: Vec<[i16; 8]>,
}

impl Residue {
    /// Check a configuration against the stream's codebooks, as header
    /// decode does.
    pub fn validate(&self, codebooks: &[Codebook]) -> Result<()> {
        if self.residue_type > 2 {
            return Err(invalid("residue type above 2"));
        }
        if self.partition_size == 0 || !(1..=64).contains(&self.classifications) {
            return Err(invalid(
                "residue partition size or classifications out of range",
            ));
        }
        if self.cascade.len() != self.classifications as usize
            || self.books.len() != self.classifications as usize
        {
            return Err(invalid(
                "residue cascade does not match its classifications",
            ));
        }
        let Some(cb) = codebooks.get(self.classbook as usize) else {
            return Err(invalid("residue classbook does not exist"));
        };
        if cb.dimensions == 0 {
            return Err(invalid("residue classbook with zero dimensions"));
        }
        let mut words: u64 = 1;
        for _ in 0..cb.dimensions {
            words = words.saturating_mul(self.classifications as u64);
        }
        if words > cb.entries() as u64 {
            return Err(invalid(
                "residue classifications^dimensions exceed the classbook's entries",
            ));
        }
        for (c, books) in self.books.iter().enumerate() {
            for (pass, &b) in books.iter().enumerate() {
                let coded = self.cascade[c] >> pass & 1 == 1;
                if coded != (b >= 0) {
                    return Err(invalid("residue book list does not match the cascade"));
                }
                if b >= 0 {
                    let Some(book) = codebooks.get(b as usize) else {
                        return Err(invalid("residue book does not exist"));
                    };
                    if !book.has_lookup() {
                        return Err(invalid("residue book without a value mapping"));
                    }
                }
            }
        }
        Ok(())
    }

    /// Header decode (8.6.1), after the 16-bit type.
    pub(crate) fn read(
        r: &mut BitReader,
        residue_type: u16,
        codebooks: &[Codebook],
    ) -> Result<Self> {
        let eop = |_| invalid("setup header ends inside a residue configuration");
        let begin = r.read(24).map_err(eop)?;
        let end = r.read(24).map_err(eop)?;
        let partition_size = r.read(24).map_err(eop)? + 1;
        let classifications = r.read(6).map_err(eop)? as u8 + 1;
        let classbook = r.read(8).map_err(eop)? as u8;
        let mut cascade = Vec::new();
        for _ in 0..classifications {
            let low = r.read(3).map_err(eop)?;
            let high = if r.read_flag().map_err(eop)? {
                r.read(5).map_err(eop)?
            } else {
                0
            };
            cascade.push((high * 8 + low) as u8);
        }
        let mut books = Vec::new();
        for &c in &cascade {
            let mut row = [-1i16; 8];
            for (j, slot) in row.iter_mut().enumerate() {
                if c >> j & 1 == 1 {
                    *slot = r.read(8).map_err(eop)? as i16;
                }
            }
            books.push(row);
        }
        let residue = Residue {
            residue_type,
            begin,
            end,
            partition_size,
            classifications,
            classbook,
            cascade,
            books,
        };
        residue.validate(codebooks)?;
        Ok(residue)
    }

    pub(crate) fn write(&self, w: &mut BitWriter) {
        w.write(self.begin, 24);
        w.write(self.end, 24);
        w.write(self.partition_size - 1, 24);
        w.write(self.classifications as u32 - 1, 6);
        w.write(self.classbook as u32, 8);
        for &c in &self.cascade {
            w.write(c as u32 & 7, 3);
            let high = c as u32 >> 3;
            w.write_flag(high > 0);
            if high > 0 {
                w.write(high, 5);
            }
        }
        for row in &self.books {
            for &b in row {
                if b >= 0 {
                    w.write(b as u32, 8);
                }
            }
        }
    }

    /// Packet decode (8.6.2): `out` holds one zeroed vector of `n` values
    /// per channel of the submap; vectors flagged in `do_not_decode` are
    /// left zero (for type 2, only if every one is flagged). An
    /// end-of-packet returns what was decoded up to that point.
    pub(crate) fn decode(
        &self,
        books: &[Codebook],
        r: &mut BitReader,
        do_not_decode: &[bool],
        out: &mut [Vec<f32>],
    ) -> std::result::Result<(), EndOfPacket> {
        if out.is_empty() {
            return Ok(());
        }
        let n = out[0].len();
        if self.residue_type == 2 {
            if do_not_decode.iter().all(|&d| d) {
                return Ok(());
            }
            let ch = out.len();
            let mut v = vec![vec![0f32; n * ch]];
            let result = self.decode_format(books, r, &[false], &mut v, true);
            for (i, x) in v[0].iter().enumerate() {
                out[i % ch][i / ch] = *x;
            }
            result
        } else {
            self.decode_format(books, r, do_not_decode, out, self.residue_type == 1)
        }
    }

    /// The common decode of 8.6.2, for format 0 (`in_order` false) or
    /// format 1 (`in_order` true).
    fn decode_format(
        &self,
        books: &[Codebook],
        r: &mut BitReader,
        do_not_decode: &[bool],
        out: &mut [Vec<f32>],
        in_order: bool,
    ) -> std::result::Result<(), EndOfPacket> {
        let actual = out[0].len();
        let begin = (self.begin as usize).min(actual);
        let end = (self.end as usize).min(actual);
        let psize = self.partition_size as usize;
        let to_read = end.saturating_sub(begin) / psize;
        if to_read == 0 {
            return Ok(());
        }
        let classbook = &books[self.classbook as usize];
        let per_word = classbook.dimensions as usize;
        let classes = self.classifications as u32;
        let ch = out.len();
        let mut classifications = vec![vec![0u8; to_read + per_word]; ch];
        for pass in 0..8 {
            let mut partition = 0;
            while partition < to_read {
                if pass == 0 {
                    for j in 0..ch {
                        if do_not_decode[j] {
                            continue;
                        }
                        let mut temp = classbook.decode_scalar(r)?;
                        for i in (0..per_word).rev() {
                            classifications[j][i + partition] = (temp % classes) as u8;
                            temp /= classes;
                        }
                    }
                }
                for _ in 0..per_word {
                    if partition >= to_read {
                        break;
                    }
                    for j in 0..ch {
                        if do_not_decode[j] {
                            continue;
                        }
                        let class = classifications[j][partition] as usize;
                        let book = self.books[class][pass];
                        if book < 0 {
                            continue;
                        }
                        let book = &books[book as usize];
                        let offset = begin + partition * psize;
                        if in_order {
                            decode_partition_1(book, r, &mut out[j], offset, psize)?;
                        } else {
                            decode_partition_0(book, r, &mut out[j], offset, psize)?;
                        }
                    }
                    partition += 1;
                }
            }
        }
        Ok(())
    }
}

/// Format 0 (8.6.3): a vector's scalars are `step` apart.
fn decode_partition_0(
    book: &Codebook,
    r: &mut BitReader,
    v: &mut [f32],
    offset: usize,
    psize: usize,
) -> std::result::Result<(), EndOfPacket> {
    let dim = book.dimensions as usize;
    let step = psize / dim;
    let mut tmp = vec![0f32; dim];
    for i in 0..step {
        tmp.fill(0.0);
        book.decode_vector_add(r, &mut tmp)?;
        for (j, t) in tmp.iter().enumerate() {
            if let Some(x) = v.get_mut(offset + i + j * step) {
                *x += *t;
            }
        }
    }
    Ok(())
}

/// Format 1 (8.6.4): a vector's scalars are consecutive.
fn decode_partition_1(
    book: &Codebook,
    r: &mut BitReader,
    v: &mut [f32],
    offset: usize,
    psize: usize,
) -> std::result::Result<(), EndOfPacket> {
    let dim = book.dimensions as usize;
    let mut i = 0;
    while i < psize {
        let at = offset + i;
        if at + dim <= v.len() {
            book.decode_vector_add(r, &mut v[at..at + dim])?;
        } else {
            let mut tmp = vec![0f32; dim];
            book.decode_vector_add(r, &mut tmp)?;
            for (k, t) in tmp.iter().enumerate() {
                if let Some(x) = v.get_mut(at + k) {
                    *x += *t;
                }
            }
        }
        i += dim;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codebook::{VqLookup, float32_pack};

    /// A scalar book of values -8..=7 (16 entries, 4 bits each), and a
    /// 2-dimensional lattice book of {-1, 0, 1} (9 entries).
    fn books() -> Vec<Codebook> {
        let scalar = Codebook::new(
            1,
            vec![4; 16],
            Some(VqLookup {
                lookup_type: 1,
                minimum: float32_pack(-8.0),
                delta: float32_pack(1.0),
                value_bits: 4,
                sequence_p: false,
                multiplicands: (0..16).collect(),
            }),
        )
        .unwrap();
        let pair = Codebook::new(
            2,
            vec![3, 3, 3, 3, 3, 3, 3, 4, 4],
            Some(VqLookup {
                lookup_type: 1,
                minimum: float32_pack(-1.0),
                delta: float32_pack(1.0),
                value_bits: 2,
                sequence_p: false,
                multiplicands: vec![0, 1, 2],
            }),
        )
        .unwrap();
        // Classbook: 2 classifications, 2 per word: 4 entries.
        let classbook = Codebook::new(2, vec![2, 2, 2, 2], None).unwrap();
        // A 4-dimensional type-2 book whose entries are arbitrary vectors.
        let quad = Codebook::new(
            4,
            vec![1, 1],
            Some(VqLookup {
                lookup_type: 2,
                minimum: float32_pack(0.0),
                delta: float32_pack(1.0),
                value_bits: 4,
                sequence_p: false,
                multiplicands: vec![1, 2, 3, 4, 5, 6, 7, 8],
            }),
        )
        .unwrap();
        vec![scalar, pair, classbook, quad]
    }

    fn config(
        residue_type: u16,
        begin: u32,
        end: u32,
        psize: u32,
        pass_books: [i16; 8],
    ) -> Residue {
        let mut cascade1 = 0u8;
        for (j, &b) in pass_books.iter().enumerate() {
            if b >= 0 {
                cascade1 |= 1 << j;
            }
        }
        Residue {
            residue_type,
            begin,
            end,
            partition_size: psize,
            classifications: 2,
            classbook: 2,
            cascade: vec![0, cascade1],
            books: vec![[-1; 8], pass_books],
        }
    }

    /// Write classwords and partition vectors exactly as 8.6.2 lays them
    /// out, by hand, for one or more channels: `classes[ch][partition]`,
    /// `entries[pass][ch][partition]` (codebook entries per partition).
    fn pack(
        books: &[Codebook],
        res: &Residue,
        classes: &[Vec<u8>],
        entries: &[Vec<Vec<Vec<u32>>>],
    ) -> Vec<u8> {
        let mut w = BitWriter::new();
        let parts = classes[0].len();
        for (pass, per_ch) in entries.iter().enumerate() {
            let mut p = 0;
            while p < parts {
                if pass == 0 {
                    for c in classes {
                        let word = c[p] as u32 * 2 + c.get(p + 1).copied().unwrap_or(0) as u32;
                        books[2].write_entry(&mut w, word);
                    }
                }
                for _ in 0..2 {
                    if p >= parts {
                        break;
                    }
                    for (ch, c) in classes.iter().enumerate() {
                        if res.books[c[p] as usize][pass] >= 0 {
                            for &e in &per_ch[ch][p] {
                                books[res.books[c[p] as usize][pass] as usize]
                                    .write_entry(&mut w, e);
                            }
                        }
                    }
                    p += 1;
                }
            }
        }
        w.into_bytes()
    }

    /// Format 1: scalars in order; begin/end bound the coded range.
    #[test]
    fn residue_1_in_order() {
        let books = books();
        let res = config(1, 2, 10, 4, [1, -1, -1, -1, -1, -1, -1, -1]);
        res.validate(&books).unwrap();
        // Two partitions [2,6) and [6,10); the second classified 0 (silent).
        // Partition 0: pair entries 0 (-1,-1) and 5 (1, 0).
        let bytes = pack(
            &books,
            &res,
            &[vec![1, 0]],
            &[vec![vec![vec![0, 5], vec![]]]],
        );
        let mut out = vec![vec![0f32; 12]];
        res.decode(&books, &mut BitReader::new(&bytes), &[false], &mut out)
            .unwrap();
        assert_eq!(
            out[0],
            vec![0.0, 0.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        );
    }

    /// Format 0 interleaves by the codebook dimension: with a 4-dimensional
    /// book on an 8-value partition, the first vector fills 0, 2, 4, 6.
    #[test]
    fn residue_0_interleaves() {
        let books = books();
        let res = config(0, 0, 8, 8, [3, -1, -1, -1, -1, -1, -1, -1]);
        let bytes = pack(&books, &res, &[vec![1]], &[vec![vec![vec![0, 1]]]]);
        let mut out = vec![vec![0f32; 8]];
        res.decode(&books, &mut BitReader::new(&bytes), &[false], &mut out)
            .unwrap();
        // Entry 0 = (1,2,3,4) at 0,2,4,6; entry 1 = (5,6,7,8) at 1,3,5,7.
        assert_eq!(out[0], vec![1.0, 5.0, 2.0, 6.0, 3.0, 7.0, 4.0, 8.0]);
    }

    /// Several passes add; a "do not decode" vector stays zero and costs
    /// no bits.
    #[test]
    fn passes_accumulate_and_flags_skip() {
        let books = books();
        let res = config(1, 0, 4, 4, [1, 0, -1, -1, -1, -1, -1, -1]);
        // Channel 1 is "do not decode": only channel 0 is packed.
        let bytes = pack(
            &books,
            &res,
            &[vec![1]],
            &[vec![vec![vec![8, 0]]], vec![vec![vec![15, 0, 9, 8]]]],
        );
        let mut out = vec![vec![0f32; 4], vec![0f32; 4]];
        res.decode(
            &books,
            &mut BitReader::new(&bytes),
            &[false, true],
            &mut out,
        )
        .unwrap();
        // Pass 0: (1,1), (-1,-1); pass 1 scalars 7, -8, 1, 0.
        assert_eq!(out[0], vec![8.0, -7.0, 0.0, -1.0]);
        assert_eq!(out[1], vec![0.0; 4]);
    }

    /// Format 2 decodes one interleaved vector and deinterleaves it.
    #[test]
    fn residue_2_interleaves_channels() {
        let books = books();
        let res = config(2, 0, 8, 4, [1, -1, -1, -1, -1, -1, -1, -1]);
        // One vector of 8 values = 2 channels x 4.
        let bytes = pack(
            &books,
            &res,
            &[vec![1, 1]],
            &[vec![vec![vec![0, 2], vec![6, 4]]]],
        );
        let mut out = vec![vec![0f32; 4], vec![0f32; 4]];
        // One flagged channel does not stop the decode.
        res.decode(
            &books,
            &mut BitReader::new(&bytes),
            &[true, false],
            &mut out,
        )
        .unwrap();
        // Interleaved: (-1,-1),(1,-1),(-1,1),(0,0) -> ch0 = -1,1,-1,0; ch1 = -1,-1,1,0.
        assert_eq!(out[0], vec![-1.0, 1.0, -1.0, 0.0]);
        assert_eq!(out[1], vec![-1.0, -1.0, 1.0, 0.0]);
        // All flagged: nothing is read.
        let mut out = vec![vec![0f32; 4], vec![0f32; 4]];
        let mut r = BitReader::new(&bytes);
        res.decode(&books, &mut r, &[true, true], &mut out).unwrap();
        assert_eq!(r.remaining(), bytes.len() * 8);
    }

    /// An end-of-packet keeps what was decoded before it.
    #[test]
    fn truncation_keeps_the_decoded_prefix() {
        let books = books();
        let res = config(1, 0, 8, 4, [0, -1, -1, -1, -1, -1, -1, -1]);
        let bytes = pack(
            &books,
            &res,
            &[vec![1, 1]],
            &[vec![vec![vec![9, 10, 11, 12], vec![13, 14, 15, 0]]]],
        );
        let cut = &bytes[..3];
        let mut out = vec![vec![0f32; 8]];
        assert!(
            res.decode(&books, &mut BitReader::new(cut), &[false], &mut out)
                .is_err()
        );
        assert_eq!(&out[0][..4], &[1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn configurations_round_trip_and_are_checked() {
        let books = books();
        let res = config(2, 0, 8, 4, [1, 0, -1, -1, -1, -1, -1, 3]);
        res.validate(&books).unwrap();
        let mut w = BitWriter::new();
        res.write(&mut w);
        let bytes = w.into_bytes();
        let back = Residue::read(&mut BitReader::new(&bytes), 2, &books).unwrap();
        assert_eq!(back, res);
        // A book without a lookup table is refused.
        assert!(
            config(1, 0, 8, 4, [2, -1, -1, -1, -1, -1, -1, -1])
                .validate(&books)
                .is_err()
        );
        // Three classifications cannot fit two per word in four entries.
        let mut bad = res.clone();
        bad.classifications = 3;
        bad.cascade.push(0);
        bad.books.push([-1; 8]);
        assert!(bad.validate(&books).is_err());
    }
}
