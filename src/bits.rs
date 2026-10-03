//! The Vorbis bit-packing convention (specification section 2): fields of
//! 0 to 32 bits, packed least-significant bit first into each byte, bytes
//! in order.

/// The end-of-packet condition (specification 2.1.8). Not always an error:
/// audio packets may be truncated on purpose, and the decoder decides what
/// running out of bits means at each point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct EndOfPacket;

/// Reads fields from one packet.
pub(crate) struct BitReader<'a> {
    data: &'a [u8],
    /// Bit position of the next read.
    pos: usize,
    /// Set once a read has run past the end; every later read fails too
    /// (2.1.8, 2.1.9).
    eop: bool,
}

impl<'a> BitReader<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        BitReader { data, pos: 0, eop: false }
    }

    /// Bits not yet read.
    pub(crate) fn remaining(&self) -> usize {
        if self.eop { 0 } else { self.data.len() * 8 - self.pos }
    }

    /// Read an unsigned field of `n` bits, `n` at most 32.
    pub(crate) fn read(&mut self, n: u32) -> Result<u32, EndOfPacket> {
        debug_assert!(n <= 32);
        if self.eop {
            return Err(EndOfPacket);
        }
        if n == 0 {
            return Ok(0);
        }
        let end = self.pos + n as usize;
        if end > self.data.len() * 8 {
            self.set_eop();
            return Err(EndOfPacket);
        }
        let v = self.peek_raw(n);
        self.pos = end;
        Ok(v)
    }

    /// One bit, as a flag.
    pub(crate) fn read_flag(&mut self) -> Result<bool, EndOfPacket> {
        Ok(self.read(1)? == 1)
    }

    /// The next `n` bits (at most 32) without consuming them, the first bit
    /// of the stream in bit 0; bits past the end read as zero. Callers check
    /// [`remaining`](Self::remaining) before trusting the high bits.
    pub(crate) fn peek(&self, n: u32) -> u32 {
        if self.eop || n == 0 { 0 } else { self.peek_raw(n) }
    }

    fn peek_raw(&self, n: u32) -> u32 {
        let byte = self.pos / 8;
        let shift = self.pos % 8;
        let mut v: u64 = 0;
        for i in 0..5 {
            if let Some(&b) = self.data.get(byte + i) {
                v |= (b as u64) << (8 * i);
            }
        }
        v >>= shift;
        let mask = if n == 32 { u32::MAX as u64 } else { (1u64 << n) - 1 };
        (v & mask) as u32
    }

    /// Consume `n` bits already known to be present (after a [`peek`](Self::peek)).
    pub(crate) fn skip(&mut self, n: u32) {
        debug_assert!(self.pos + n as usize <= self.data.len() * 8);
        self.pos += n as usize;
    }

    /// Mark the packet exhausted (a read that ran out of bits).
    pub(crate) fn set_eop(&mut self) {
        self.pos = self.data.len() * 8;
        self.eop = true;
    }
}

/// Writes fields into a packet; the unused bits of the last byte are zero
/// (2.1.8).
#[derive(Default)]
pub(crate) struct BitWriter {
    data: Vec<u8>,
    /// Bits used in the last byte (0 means the next bit starts a new byte).
    used: u32,
}

impl BitWriter {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Write the low `n` bits of `value`, `n` at most 32.
    pub(crate) fn write(&mut self, value: u32, n: u32) {
        debug_assert!(n <= 32);
        let mut v = if n == 32 { value as u64 } else { (value as u64) & ((1u64 << n) - 1) };
        let mut left = n;
        while left > 0 {
            if self.used == 0 {
                self.data.push(0);
            }
            let room = 8 - self.used;
            let take = room.min(left);
            let last = self.data.len() - 1;
            self.data[last] |= ((v & ((1u64 << take) - 1)) as u8) << self.used;
            v >>= take;
            left -= take;
            self.used = (self.used + take) % 8;
        }
    }

    pub(crate) fn write_flag(&mut self, flag: bool) {
        self.write(flag as u32, 1);
    }

    /// Write a Huffman codeword: its bits go out first-bit-first (3.2.1),
    /// `code` holding the first bit in its most significant used position.
    pub(crate) fn write_codeword(&mut self, code: u32, len: u32) {
        let reversed = if len == 0 { 0 } else { code.reverse_bits() >> (32 - len) };
        self.write(reversed, len);
    }

    /// Bits written so far.
    pub(crate) fn bit_len(&self) -> usize {
        if self.used == 0 {
            self.data.len() * 8
        } else {
            (self.data.len() - 1) * 8 + self.used as usize
        }
    }

    #[cfg(test)]
    pub(crate) fn bytes(&self) -> &[u8] {
        &self.data
    }

    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.data
    }
}

/// `ilog` (9.2.1): the position of the highest set bit, 0 for zero and
/// negative values.
pub(crate) fn ilog(x: i64) -> u32 {
    if x <= 0 { 0 } else { 64 - (x as u64).leading_zeros() }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The coding example of 2.1.6 and the decoding example of 2.1.7.
    #[test]
    fn spec_bitpacking_example() {
        let mut w = BitWriter::new();
        w.write(12, 4);
        assert_eq!(w.bytes(), &[0b0000_1100]);
        w.write((-1i32) as u32, 3);
        assert_eq!(w.bytes(), &[0b0111_1100]);
        w.write(17, 7);
        assert_eq!(w.bytes(), &[0b1111_1100, 0b0000_1000]);
        assert_eq!(w.bit_len(), 14);
        w.write(6969, 13);
        let bytes = w.into_bytes();
        assert_eq!(bytes, vec![0b1111_1100, 0b0100_1000, 0b1100_1110, 0b0000_0110]);

        let mut r = BitReader::new(&bytes);
        assert_eq!(r.read(2), Ok(0b00));
        assert_eq!(r.read(2), Ok(0b11));
        let mut r = BitReader::new(&bytes);
        assert_eq!(r.read(4), Ok(12));
        assert_eq!(r.read(3), Ok(7));
        assert_eq!(r.read(7), Ok(17));
        assert_eq!(r.read(13), Ok(6969));
        // Five zero bits of padding remain, then end-of-packet.
        assert_eq!(r.remaining(), 5);
        assert_eq!(r.read(5), Ok(0));
        // 2.1.9: reading zero bits at the end succeeds...
        assert_eq!(r.read(0), Ok(0));
        // ...reading past it does not, and the condition is sticky.
        assert_eq!(r.read(1), Err(EndOfPacket));
        assert_eq!(r.read(0), Err(EndOfPacket));
    }

    #[test]
    fn partial_read_past_end_is_end_of_packet() {
        let mut r = BitReader::new(&[0xff]);
        assert_eq!(r.read(6), Ok(0x3f));
        assert_eq!(r.read(3), Err(EndOfPacket));
        assert_eq!(r.read(1), Err(EndOfPacket));
    }

    #[test]
    fn thirty_two_bit_fields_round_trip() {
        let mut w = BitWriter::new();
        w.write(1, 1);
        w.write(0xdead_beef, 32);
        w.write(0x8000_0001, 32);
        let bytes = w.into_bytes();
        let mut r = BitReader::new(&bytes);
        assert_eq!(r.read(1), Ok(1));
        assert_eq!(r.read(32), Ok(0xdead_beef));
        assert_eq!(r.read(32), Ok(0x8000_0001));
    }

    #[test]
    fn codewords_go_out_first_bit_first() {
        // Codeword "001" (first bit 0, then 0, then 1).
        let mut w = BitWriter::new();
        w.write_codeword(0b001, 3);
        let bytes = w.into_bytes();
        let mut r = BitReader::new(&bytes);
        assert_eq!(r.read(1), Ok(0));
        assert_eq!(r.read(1), Ok(0));
        assert_eq!(r.read(1), Ok(1));
    }

    #[test]
    fn ilog_matches_the_spec_examples() {
        assert_eq!(ilog(0), 0);
        assert_eq!(ilog(1), 1);
        assert_eq!(ilog(2), 2);
        assert_eq!(ilog(3), 2);
        assert_eq!(ilog(4), 3);
        assert_eq!(ilog(7), 3);
        assert_eq!(ilog(-1), 0);
    }
}
