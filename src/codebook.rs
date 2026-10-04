//! Codebooks (specification section 3): the packed format, the Huffman
//! decision tree the codeword lengths define, and the VQ lookup tables.

use crate::bits::{BitReader, BitWriter, EndOfPacket, ilog};
use crate::error::{Result, invalid};

/// Codebooks larger than this many scalars compute their VQ vectors on
/// demand rather than unpacking them at setup (a guard against headers
/// that declare millions of entries).
const PRECOMPUTE_LIMIT: u64 = 1 << 20;

/// Bits the first-level decode table covers.
const TABLE_BITS: u32 = 10;

/// A child slot holding a leaf (an entry number) rather than a node index.
const LEAF: u32 = 0x8000_0000;

/// The VQ lookup table of a codebook (lookup types 1 and 2, 3.2.1).
#[derive(Clone, Debug, PartialEq)]
pub struct VqLookup {
    /// 1 (lattice: values permuted from a short list) or 2 (one value per
    /// scalar of every entry).
    pub lookup_type: u8,
    /// `codebook_minimum_value`, packed (see [`float32_unpack`]).
    pub minimum: u32,
    /// `codebook_delta_value`, packed.
    pub delta: u32,
    /// Bits per multiplicand, 1 to 16.
    pub value_bits: u8,
    /// `codebook_sequence_p`: each scalar adds the one before it.
    pub sequence_p: bool,
    /// `codebook_multiplicands`.
    pub multiplicands: Vec<u32>,
}

/// One codebook: codeword lengths (0 marks an unused entry), a dimension,
/// and optionally a VQ lookup table.
#[derive(Clone, Debug)]
pub struct Codebook {
    /// `codebook_dimensions`: scalars per VQ vector.
    pub dimensions: u16,
    /// Codeword length of each entry, 1 to 32, or 0 for an unused entry of
    /// a sparse codebook. Its length is `codebook_entries`.
    pub lengths: Vec<u8>,
    /// The VQ lookup table; `None` is lookup type 0.
    pub lookup: Option<VqLookup>,
    huffman: Huffman,
    /// Unpacked VQ vectors, `entries * dimensions` scalars, when small enough.
    values: Option<Vec<f32>>,
    /// `lookup1_values` or the multiplicand count.
    lookup_values: u32,
}

impl PartialEq for Codebook {
    fn eq(&self, other: &Self) -> bool {
        self.dimensions == other.dimensions
            && self.lengths == other.lengths
            && self.lookup == other.lookup
    }
}

/// `float32_unpack` (9.2.2): a 21-bit mantissa, a 10-bit exponent biased by
/// 788 and a sign bit.
pub fn float32_unpack(x: u32) -> f32 {
    let mantissa = (x & 0x1f_ffff) as f64;
    let exponent = ((x & 0x7fe0_0000) >> 21) as i32;
    let m = if x & 0x8000_0000 != 0 {
        -mantissa
    } else {
        mantissa
    };
    (m * 2f64.powi(exponent - 788)) as f32
}

/// The inverse of [`float32_unpack`] for any finite value: exact for every
/// value with a 21-bit mantissa in range (integers below 2^21 among them).
pub fn float32_pack(v: f64) -> u32 {
    if v == 0.0 || !v.is_finite() {
        return 0;
    }
    let sign = if v < 0.0 { 0x8000_0000 } else { 0 };
    let m = v.abs();
    let mut e = m.log2().floor() as i32;
    let mut mant = (m / 2f64.powi(e - 20)).round() as u64;
    if mant >= 1 << 21 {
        mant >>= 1;
        e += 1;
    }
    let exp = (e - 20 + 788).clamp(0, 1023) as u32;
    sign | (exp << 21) | (mant as u32 & 0x1f_ffff)
}

/// `lookup1_values` (9.2.3): the greatest `r` with `r^dimensions <= entries`.
pub fn lookup1_values(entries: u32, dimensions: u16) -> u32 {
    if dimensions == 0 {
        return 0;
    }
    let pow = |r: u64| -> u64 {
        let mut acc: u64 = 1;
        for _ in 0..dimensions {
            acc = acc.saturating_mul(r);
            if acc > u32::MAX as u64 {
                return u64::MAX;
            }
        }
        acc
    };
    let mut r = (entries as f64).powf(1.0 / dimensions as f64).floor() as u64;
    while r > 0 && pow(r) > entries as u64 {
        r -= 1;
    }
    while pow(r + 1) <= entries as u64 {
        r += 1;
    }
    r as u32
}

impl Codebook {
    /// Build a codebook from its fields, checking them as header decode
    /// does (a complete Huffman tree; a lookup table of the right size).
    pub fn new(dimensions: u16, lengths: Vec<u8>, lookup: Option<VqLookup>) -> Result<Self> {
        if lengths.len() >= 1 << 24 {
            return Err(invalid("codebook has more than 2^24 - 1 entries"));
        }
        if lengths.iter().any(|&l| l > 32) {
            return Err(invalid("codeword longer than 32 bits"));
        }
        let entries = lengths.len() as u32;
        let mut lookup_values = 0;
        if let Some(l) = &lookup {
            if dimensions == 0 {
                return Err(invalid("VQ codebook with zero dimensions"));
            }
            if !(1..=16).contains(&l.value_bits) {
                return Err(invalid("VQ value bits outside 1..=16"));
            }
            let want = match l.lookup_type {
                1 => lookup1_values(entries, dimensions) as u64,
                2 => entries as u64 * dimensions as u64,
                t => return Err(invalid(format!("codebook lookup type {t}"))),
            };
            if l.multiplicands.len() as u64 != want {
                return Err(invalid(
                    "VQ multiplicand count does not match the lookup type",
                ));
            }
            if l.multiplicands.iter().any(|&m| m >> l.value_bits != 0) {
                return Err(invalid("VQ multiplicand wider than its value bits"));
            }
            lookup_values = want as u32;
        }
        let huffman = Huffman::build(&lengths)?;
        let mut book = Codebook {
            dimensions,
            lengths,
            lookup,
            huffman,
            values: None,
            lookup_values,
        };
        if book.lookup.is_some() && entries as u64 * dimensions as u64 <= PRECOMPUTE_LIMIT {
            let d = dimensions as usize;
            let mut values = vec![0f32; entries as usize * d];
            for e in 0..entries {
                book.compute_vector(e, &mut values[e as usize * d..(e as usize + 1) * d]);
            }
            book.values = Some(values);
        }
        Ok(book)
    }

    /// `codebook_entries`.
    pub fn entries(&self) -> u32 {
        self.lengths.len() as u32
    }

    /// Decode one codebook from the setup header (3.2.1).
    pub(crate) fn read(r: &mut BitReader) -> Result<Self> {
        let eop = |_| invalid("setup header ends inside a codebook");
        let sync = r.read(24).map_err(eop)?;
        if sync != 0x56_4342 {
            return Err(invalid("codebook sync pattern missing"));
        }
        let dimensions = r.read(16).map_err(eop)? as u16;
        let entries = r.read(24).map_err(eop)?;
        let ordered = r.read_flag().map_err(eop)?;
        let mut lengths;
        if !ordered {
            let sparse = r.read_flag().map_err(eop)?;
            // Each entry costs at least one bit: refuse to allocate for a
            // count the packet cannot hold.
            if entries as usize > r.remaining() {
                return Err(invalid("setup header ends inside a codebook"));
            }
            lengths = vec![0u8; entries as usize];
            for l in lengths.iter_mut() {
                let used = if sparse {
                    r.read_flag().map_err(eop)?
                } else {
                    true
                };
                if used {
                    *l = r.read(5).map_err(eop)? as u8 + 1;
                }
            }
        } else {
            lengths = vec![0u8; entries as usize];
            let mut current = 0u32;
            let mut length = r.read(5).map_err(eop)? + 1;
            while current < entries {
                let number = r.read(ilog((entries - current) as i64)).map_err(eop)?;
                if current + number > entries {
                    return Err(invalid("ordered codebook lists more entries than it has"));
                }
                if number > 0 && length > 32 {
                    return Err(invalid("codeword longer than 32 bits"));
                }
                for l in &mut lengths[current as usize..(current + number) as usize] {
                    *l = length as u8;
                }
                current += number;
                length += 1;
            }
        }
        let lookup_type = r.read(4).map_err(eop)? as u8;
        let lookup = match lookup_type {
            0 => None,
            1 | 2 => {
                let minimum = r.read(32).map_err(eop)?;
                let delta = r.read(32).map_err(eop)?;
                let value_bits = r.read(4).map_err(eop)? as u8 + 1;
                let sequence_p = r.read_flag().map_err(eop)?;
                if dimensions == 0 {
                    return Err(invalid("VQ codebook with zero dimensions"));
                }
                let count = if lookup_type == 1 {
                    lookup1_values(entries, dimensions) as u64
                } else {
                    entries as u64 * dimensions as u64
                };
                if count * value_bits as u64 > r.remaining() as u64 {
                    return Err(invalid(
                        "setup header ends inside a codebook's lookup table",
                    ));
                }
                let mut multiplicands = Vec::with_capacity(count as usize);
                for _ in 0..count {
                    multiplicands.push(r.read(value_bits as u32).map_err(eop)?);
                }
                Some(VqLookup {
                    lookup_type,
                    minimum,
                    delta,
                    value_bits,
                    sequence_p,
                    multiplicands,
                })
            }
            t => return Err(invalid(format!("codebook lookup type {t} is reserved"))),
        };
        Codebook::new(dimensions, lengths, lookup)
    }

    /// Pack this codebook as the setup header carries it, choosing the
    /// ordered or sparse length encoding where it applies.
    pub(crate) fn write(&self, w: &mut BitWriter) {
        let entries = self.entries();
        w.write(0x56_4342, 24);
        w.write(self.dimensions as u32, 16);
        w.write(entries, 24);
        let all_used = self.lengths.iter().all(|&l| l > 0);
        let ordered = all_used && entries > 0 && self.lengths.windows(2).all(|p| p[0] <= p[1]);
        w.write_flag(ordered);
        if ordered {
            let mut current = 0u32;
            let mut length = self.lengths[0] as u32;
            w.write(length - 1, 5);
            while current < entries {
                let number = self.lengths[current as usize..]
                    .iter()
                    .take_while(|&&l| l as u32 == length)
                    .count() as u32;
                w.write(number, ilog((entries - current) as i64));
                current += number;
                length += 1;
            }
        } else {
            w.write_flag(!all_used);
            for &l in &self.lengths {
                if !all_used {
                    w.write_flag(l > 0);
                }
                if l > 0 {
                    w.write(l as u32 - 1, 5);
                }
            }
        }
        match &self.lookup {
            None => w.write(0, 4),
            Some(l) => {
                w.write(l.lookup_type as u32, 4);
                w.write(l.minimum, 32);
                w.write(l.delta, 32);
                w.write(l.value_bits as u32 - 1, 4);
                w.write_flag(l.sequence_p);
                for &m in &l.multiplicands {
                    w.write(m, l.value_bits as u32);
                }
            }
        }
    }

    /// Read one codeword and return its entry number (scalar context).
    #[inline]
    pub(crate) fn decode_scalar(&self, r: &mut BitReader) -> std::result::Result<u32, EndOfPacket> {
        self.huffman.decode(r)
    }

    /// Whether the book can be used in VQ context.
    pub fn has_lookup(&self) -> bool {
        self.lookup.is_some()
    }

    /// Read one codeword in VQ context and add its vector to `out` (whose
    /// length is the dimension). Calls on a book without a lookup table are
    /// the caller's error; it checks [`has_lookup`](Self::has_lookup).
    #[inline]
    pub(crate) fn decode_vector_add(
        &self,
        r: &mut BitReader,
        out: &mut [f32],
    ) -> std::result::Result<(), EndOfPacket> {
        let e = self.huffman.decode(r)? as usize;
        let d = self.dimensions as usize;
        if let Some(values) = &self.values {
            for (o, v) in out.iter_mut().zip(&values[e * d..(e + 1) * d]) {
                *o += *v;
            }
        } else {
            let mut tmp = vec![0f32; d];
            self.compute_vector(e as u32, &mut tmp);
            for (o, v) in out.iter_mut().zip(&tmp) {
                *o += *v;
            }
        }
        Ok(())
    }

    /// The VQ vector of `entry` (3.2.1, "VQ lookup table vector
    /// representation"); empty for a book without a lookup table.
    pub fn vector(&self, entry: u32) -> Vec<f32> {
        if self.lookup.is_none() {
            return Vec::new();
        }
        let mut out = vec![0f32; self.dimensions as usize];
        self.compute_vector(entry, &mut out);
        out
    }

    fn compute_vector(&self, entry: u32, out: &mut [f32]) {
        let Some(l) = &self.lookup else { return };
        let minimum = float32_unpack(l.minimum);
        let delta = float32_unpack(l.delta);
        let mut last = 0f32;
        match l.lookup_type {
            1 => {
                let mut divisor: u64 = 1;
                for o in out.iter_mut() {
                    let offset = ((entry as u64 / divisor) % self.lookup_values as u64) as usize;
                    let v = l.multiplicands[offset] as f32 * delta + minimum + last;
                    *o = v;
                    if l.sequence_p {
                        last = v;
                    }
                    divisor = divisor.saturating_mul(self.lookup_values as u64);
                }
            }
            _ => {
                let base = entry as usize * self.dimensions as usize;
                for (i, o) in out.iter_mut().enumerate() {
                    let v = l.multiplicands[base + i] as f32 * delta + minimum + last;
                    *o = v;
                    if l.sequence_p {
                        last = v;
                    }
                }
            }
        }
    }

    /// The codeword of `entry`: its bits with the first in the most
    /// significant used position, and its length (0 for an unused entry).
    pub fn codeword(&self, entry: u32) -> (u32, u8) {
        (
            self.huffman.codes[entry as usize],
            self.lengths[entry as usize],
        )
    }

    /// Write the codeword of `entry`.
    #[inline]
    pub(crate) fn write_entry(&self, w: &mut BitWriter, entry: u32) {
        let (code, len) = self.codeword(entry);
        debug_assert!(len > 0, "writing an unused codebook entry");
        w.write_codeword(code, len as u32);
    }
}

/// The decision tree of 3.2.1, with a lookup table for the short codewords.
#[derive(Clone, Debug, Default)]
struct Huffman {
    /// Codeword per entry (first bit most significant), 0 for unused.
    codes: Vec<u32>,
    /// Indexed by the next `table_bits` stream bits (first bit in bit 0):
    /// `entry << 6 | length`, or 0 when the codeword is longer.
    table: Vec<u32>,
    table_bits: u32,
    /// Internal nodes: child slots are 0 (absent), `LEAF | entry`, or a node index.
    nodes: Vec<[u32; 2]>,
    /// The entry of a single-entry codebook (the 2015 errata).
    single: Option<u32>,
    /// No used entries: nothing can be read.
    empty: bool,
}

impl Huffman {
    /// Assign every used entry, in order, the lowest-valued free codeword
    /// of its length; refuse an over- or under-specified tree.
    fn build(lengths: &[u8]) -> Result<Self> {
        let used: Vec<usize> = (0..lengths.len()).filter(|&i| lengths[i] > 0).collect();
        let mut codes = vec![0u32; lengths.len()];
        match used.len() {
            // A codebook with no used entry at all is, strictly, an
            // underspecified tree; streams carry them (the Xiph.Org
            // one-entry-codebook vector does) in books nothing reads.
            // Accept it; reading from it is an end-of-packet.
            0 => {
                return Ok(Huffman {
                    codes,
                    empty: true,
                    ..Default::default()
                });
            }
            1 => {
                // Errata 20150226: a single used entry must declare length 1;
                // it reads one bit, whatever its value.
                if lengths[used[0]] != 1 {
                    return Err(invalid(
                        "single-entry codebook whose codeword length is not 1",
                    ));
                }
                return Ok(Huffman {
                    codes,
                    single: Some(used[0] as u32),
                    ..Default::default()
                });
            }
            _ => {}
        }
        const INF: u8 = u8::MAX;
        // Per internal node: its children and the shallowest depth at which
        // its subtree still has a free slot.
        let mut nodes: Vec<[u32; 2]> = vec![[0, 0]];
        let mut shallow: Vec<u8> = vec![1];
        let mut path: Vec<(usize, u8)> = Vec::with_capacity(33);
        for &e in &used {
            let len = lengths[e];
            if shallow[0] > len {
                return Err(invalid("overspecified Huffman tree"));
            }
            path.clear();
            let mut node = 0usize;
            let mut depth = 0u8;
            let mut code = 0u32;
            loop {
                // Leftmost child with a free slot at depth `len`.
                let fits = |slot: u32| -> bool {
                    if slot == 0 {
                        true
                    } else if slot & LEAF != 0 {
                        false
                    } else {
                        shallow[slot as usize] <= len
                    }
                };
                let b = if fits(nodes[node][0]) { 0 } else { 1 };
                debug_assert!(fits(nodes[node][b]));
                path.push((node, depth));
                code = (code << 1) | b as u32;
                if depth + 1 == len {
                    debug_assert_eq!(nodes[node][b], 0);
                    nodes[node][b] = LEAF | e as u32;
                    break;
                }
                if nodes[node][b] == 0 {
                    nodes.push([0, 0]);
                    shallow.push(depth + 2);
                    nodes[node][b] = (nodes.len() - 1) as u32;
                }
                node = nodes[node][b] as usize;
                depth += 1;
            }
            codes[e] = code;
            for &(n, d) in path.iter().rev() {
                let mut s = INF;
                for &slot in &nodes[n] {
                    let v = if slot == 0 {
                        d + 1
                    } else if slot & LEAF != 0 {
                        INF
                    } else {
                        shallow[slot as usize]
                    };
                    s = s.min(v);
                }
                shallow[n] = s;
            }
        }
        if shallow[0] != INF {
            return Err(invalid("underspecified Huffman tree"));
        }
        let max_len = used.iter().map(|&e| lengths[e] as u32).max().unwrap_or(1);
        let table_bits = max_len.min(TABLE_BITS);
        let mut table = vec![0u32; 1 << table_bits];
        for &e in &used {
            let len = lengths[e] as u32;
            if len > table_bits {
                continue;
            }
            let reversed = codes[e].reverse_bits() >> (32 - len);
            for fill in 0..(1u32 << (table_bits - len)) {
                table[(reversed | (fill << len)) as usize] = ((e as u32) << 6) | len;
            }
        }
        Ok(Huffman {
            codes,
            table,
            table_bits,
            nodes,
            single: None,
            empty: false,
        })
    }

    #[inline]
    fn decode(&self, r: &mut BitReader) -> std::result::Result<u32, EndOfPacket> {
        if let Some(e) = self.single {
            r.read(1)?;
            return Ok(e);
        }
        if self.empty {
            r.set_eop();
            return Err(EndOfPacket);
        }
        let rem = r.remaining();
        let t = self.table[r.peek(self.table_bits) as usize];
        let len = t & 63;
        if len != 0 && len as usize <= rem {
            r.skip(len);
            return Ok(t >> 6);
        }
        let mut node = 0usize;
        loop {
            let b = r.read(1)? as usize;
            let slot = self.nodes[node][b];
            if slot & LEAF != 0 {
                return Ok(slot & !LEAF);
            }
            if slot == 0 {
                // Unreachable for a complete tree; treated as end of packet.
                r.set_eop();
                return Err(EndOfPacket);
            }
            node = slot as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn code_string(book: &Codebook, e: u32) -> String {
        let (code, len) = book.codeword(e);
        (0..len)
            .rev()
            .map(|i| if code >> i & 1 == 1 { '1' } else { '0' })
            .collect()
    }

    /// The worked example of 3.2.1.
    #[test]
    fn spec_huffman_example() {
        let book = Codebook::new(1, vec![2, 4, 4, 4, 4, 2, 3, 3], None).unwrap();
        let want = ["00", "0100", "0101", "0110", "0111", "10", "110", "111"];
        for (e, w) in want.iter().enumerate() {
            assert_eq!(code_string(&book, e as u32), *w, "entry {e}");
        }
    }

    #[test]
    fn under_and_overspecified_trees_are_refused() {
        // The spec's example without entry seven: underspecified.
        assert!(Codebook::new(1, vec![2, 4, 4, 4, 4, 2, 3], None).is_err());
        // With a ninth codeword: overspecified.
        assert!(Codebook::new(1, vec![2, 4, 4, 4, 4, 2, 3, 3, 3], None).is_err());
        assert!(Codebook::new(1, vec![1, 1, 1], None).is_err());
        // No used entry: accepted (streams carry such books), unreadable.
        let empty = Codebook::new(1, vec![0, 0], None).unwrap();
        assert_eq!(
            empty.decode_scalar(&mut BitReader::new(&[0xff])),
            Err(EndOfPacket)
        );
    }

    #[test]
    fn lowest_free_codeword_may_lie_left_of_a_shorter_one() {
        // Lengths 2, 1, 2: "00", then "1", then "01".
        let book = Codebook::new(1, vec![2, 1, 2], None).unwrap();
        assert_eq!(code_string(&book, 0), "00");
        assert_eq!(code_string(&book, 1), "1");
        assert_eq!(code_string(&book, 2), "01");
    }

    #[test]
    fn sparse_entries_get_no_codeword() {
        let book = Codebook::new(1, vec![0, 1, 0, 2, 2], None).unwrap();
        assert_eq!(code_string(&book, 1), "0");
        assert_eq!(code_string(&book, 3), "10");
        assert_eq!(code_string(&book, 4), "11");
        assert_eq!(book.codeword(0).1, 0);
    }

    #[test]
    fn single_entry_codebooks_follow_the_errata() {
        let book = Codebook::new(1, vec![0, 0, 1, 0], None).unwrap();
        // Either bit value decodes to the one entry, sinking one bit.
        let mut r = BitReader::new(&[0b10]);
        assert_eq!(book.decode_scalar(&mut r), Ok(2));
        assert_eq!(book.decode_scalar(&mut r), Ok(2));
        assert_eq!(r.remaining(), 6);
        assert!(Codebook::new(1, vec![2], None).is_err());
    }

    /// Every codeword decodes to its own entry, through the table and
    /// through the tree walk (lengths past the table width).
    #[test]
    fn every_codeword_decodes_to_its_entry() {
        // A complete tree with lengths 1..=20 and two of length 20.
        let mut lengths: Vec<u8> = (1..=20).collect();
        lengths.push(20);
        let book = Codebook::new(1, lengths.clone(), None).unwrap();
        let mut w = BitWriter::new();
        for e in 0..lengths.len() as u32 {
            book.write_entry(&mut w, e);
        }
        let bytes = w.into_bytes();
        let mut r = BitReader::new(&bytes);
        for e in 0..lengths.len() as u32 {
            assert_eq!(book.decode_scalar(&mut r), Ok(e));
        }
    }

    #[test]
    fn float32_unpack_and_pack() {
        // 1.0: mantissa 1, exponent 788.
        assert_eq!(float32_unpack(788 << 21 | 1), 1.0);
        assert_eq!(float32_unpack(0x8000_0000 | 788 << 21 | 3), -3.0);
        assert_eq!(float32_unpack((788 - 2) << 21 | 1), 0.25);
        for v in [1.0, -1.0, 33.0, -231.0, 0.5, 1e-3, -7.25, 123456.0, 0.0] {
            let p = float32_pack(v);
            let back = float32_unpack(p) as f64;
            assert!((back - v).abs() <= v.abs() * 1e-6, "{v} -> {back}");
        }
        assert_eq!(float32_unpack(float32_pack(-16.0)), -16.0);
    }

    #[test]
    fn lookup1_values_is_the_integer_root() {
        assert_eq!(lookup1_values(81, 4), 3);
        assert_eq!(lookup1_values(80, 4), 2);
        assert_eq!(lookup1_values(625, 4), 5);
        assert_eq!(lookup1_values(1089, 2), 33);
        assert_eq!(lookup1_values(1, 1), 1);
        assert_eq!(lookup1_values(16_777_215, 1), 16_777_215);
        assert_eq!(lookup1_values(0, 3), 0);
    }

    fn lookup(t: u8, min: f64, delta: f64, seq: bool, mult: Vec<u32>) -> Option<VqLookup> {
        Some(VqLookup {
            lookup_type: t,
            minimum: float32_pack(min),
            delta: float32_pack(delta),
            value_bits: 8,
            sequence_p: seq,
            multiplicands: mult,
        })
    }

    /// Lookup type 1: the entry number's base-`lookup_values` digits, least
    /// significant first, select each scalar.
    #[test]
    fn lookup_type_1_vectors() {
        // 9 entries, 2 dimensions: 3 values {-1, 0, 1}.
        let book = Codebook::new(
            2,
            vec![4, 4, 4, 4, 4, 4, 4, 3, 3]
                .into_iter()
                .map(|l| l as u8)
                .collect(),
            lookup(1, -1.0, 1.0, false, vec![0, 1, 2]),
        );
        // Lengths 4x7 + 3x2 is not complete; use a complete set instead.
        assert!(book.is_err());
        let lengths = vec![3, 3, 3, 3, 3, 3, 3, 4, 4];
        let book = Codebook::new(2, lengths, lookup(1, -1.0, 1.0, false, vec![0, 1, 2])).unwrap();
        for e in 0..9u32 {
            let v = book.vector(e);
            assert_eq!(
                v,
                vec![(e % 3) as f32 - 1.0, (e / 3) as f32 - 1.0],
                "entry {e}"
            );
        }
        // sequence_p: each scalar adds the previous one.
        let lengths = vec![3, 3, 3, 3, 3, 3, 3, 4, 4];
        let book = Codebook::new(2, lengths, lookup(1, 0.5, 2.0, true, vec![0, 1, 3])).unwrap();
        // Entry 7: digits (1, 2) -> values 2.5, 6.5 + 2.5.
        assert_eq!(book.vector(7), vec![2.5, 9.0]);
    }

    /// Lookup type 2: each entry's scalars are read in order.
    #[test]
    fn lookup_type_2_vectors() {
        let mult = vec![1, 2, 3, 4, 5, 6];
        let book = Codebook::new(3, vec![1, 1], lookup(2, -2.0, 0.5, false, mult.clone())).unwrap();
        assert_eq!(book.vector(0), vec![-1.5, -1.0, -0.5]);
        assert_eq!(book.vector(1), vec![0.0, 0.5, 1.0]);
        let book = Codebook::new(3, vec![1, 1], lookup(2, -2.0, 0.5, true, mult)).unwrap();
        assert_eq!(book.vector(0), vec![-1.5, -2.5, -3.0]);
        // Wrong multiplicand count is refused.
        assert!(Codebook::new(3, vec![1, 1], lookup(2, 0.0, 1.0, false, vec![1, 2])).is_err());
    }

    /// Packing and unpacking agree, through each length encoding.
    #[test]
    fn codebooks_round_trip_through_the_packed_format() {
        let books = vec![
            // ordered
            Codebook::new(1, vec![1, 2, 3, 3], None).unwrap(),
            // ordered with a skipped length
            Codebook::new(1, vec![1, 3, 3, 3, 3], None).unwrap(),
            // plain (not ordered)
            Codebook::new(
                2,
                vec![3, 3, 3, 3, 3, 3, 3, 4, 4].into_iter().rev().collect(),
                lookup(1, -1.0, 1.0, false, vec![0, 1, 2]),
            )
            .unwrap(),
            // sparse
            Codebook::new(
                3,
                vec![0, 1, 0, 2, 2],
                lookup(2, -1.0, 0.25, true, (0..15).collect()),
            )
            .unwrap(),
            // single entry
            Codebook::new(1, vec![0, 1], None).unwrap(),
        ];
        for book in books {
            let mut w = BitWriter::new();
            book.write(&mut w);
            let bytes = w.into_bytes();
            let mut r = BitReader::new(&bytes);
            let back = Codebook::read(&mut r).unwrap();
            assert_eq!(back, book);
        }
    }
}
