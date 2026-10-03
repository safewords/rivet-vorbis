//! Floor type 0 (specification section 6): the spectral envelope as the
//! frequency response of an LSP filter, on a Bark-warped frequency axis.

use crate::bits::{BitReader, BitWriter, EndOfPacket, ilog};
use crate::codebook::Codebook;
use crate::error::{Result, invalid};

/// A floor 0 configuration (6.2.1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Floor0 {
    /// `floor0_order`: LSP coefficients per frame.
    pub order: u8,
    /// `floor0_rate`.
    pub rate: u16,
    /// `floor0_bark_map_size`.
    pub bark_map_size: u16,
    /// `floor0_amplitude_bits`.
    pub amplitude_bits: u8,
    /// `floor0_amplitude_offset`.
    pub amplitude_offset: u8,
    /// `floor0_book_list`, 1 to 16 books.
    pub books: Vec<u8>,
}

/// One frame's decoded floor 0 data.
#[derive(Clone, Debug, PartialEq)]
pub struct Floor0Frame {
    /// The nonzero amplitude.
    pub amplitude: u32,
    /// The LSP coefficients (at least `order` of them; extras are unused).
    pub coefficients: Vec<f32>,
}

/// `bark(x)` of 6.2.3 (as corrected by the 2015 errata).
pub fn bark(x: f64) -> f64 {
    13.1 * (0.00074 * x).atan() + 2.24 * (0.0000000185 * x * x).atan() + 0.0001 * x
}

impl Floor0 {
    /// Header decode (6.2.1).
    pub(crate) fn read(r: &mut BitReader, codebooks: &[Codebook]) -> Result<Self> {
        let eop = |_| invalid("setup header ends inside a floor 0 configuration");
        let order = r.read(8).map_err(eop)? as u8;
        let rate = r.read(16).map_err(eop)? as u16;
        let bark_map_size = r.read(16).map_err(eop)? as u16;
        let amplitude_bits = r.read(6).map_err(eop)? as u8;
        let amplitude_offset = r.read(8).map_err(eop)? as u8;
        let count = r.read(4).map_err(eop)? + 1;
        let mut books = Vec::new();
        for _ in 0..count {
            books.push(r.read(8).map_err(eop)? as u8);
        }
        if books.iter().any(|&b| b as usize >= codebooks.len()) {
            return Err(invalid("floor 0 book does not exist"));
        }
        if bark_map_size == 0 {
            return Err(invalid("floor 0 bark map size of zero"));
        }
        Ok(Floor0 { order, rate, bark_map_size, amplitude_bits, amplitude_offset, books })
    }

    pub(crate) fn write(&self, w: &mut BitWriter) {
        w.write(self.order as u32, 8);
        w.write(self.rate as u32, 16);
        w.write(self.bark_map_size as u32, 16);
        w.write(self.amplitude_bits as u32, 6);
        w.write(self.amplitude_offset as u32, 8);
        w.write(self.books.len() as u32 - 1, 4);
        for &b in &self.books {
            w.write(b as u32, 8);
        }
    }

    /// Packet decode (6.2.2): `None` for an unused floor. A book number past
    /// the list is an undecodable packet, reported as `Err(Some(..))`.
    pub(crate) fn decode(&self, r: &mut BitReader, books: &[Codebook]) -> std::result::Result<Option<Floor0Frame>, Floor0Error> {
        let amplitude = r.read(self.amplitude_bits as u32)?;
        if amplitude == 0 {
            return Ok(None);
        }
        let booknumber = r.read(ilog(self.books.len() as i64))? as usize;
        if booknumber >= self.books.len() {
            return Err(Floor0Error::Undecodable("floor 0 book number past the book list"));
        }
        let book = &books[self.books[booknumber] as usize];
        if !book.has_lookup() || book.dimensions == 0 {
            return Err(Floor0Error::Undecodable("floor 0 book has no VQ lookup"));
        }
        let dim = book.dimensions as usize;
        let mut coefficients: Vec<f32> = Vec::with_capacity(self.order as usize + dim);
        // The running `last` carries from one vector to the next, so the
        // LSP angles accumulate across vectors (see docs/PROVENANCE.md on
        // the reading of step 11's "continue at step 6").
        let mut last = 0f32;
        while coefficients.len() < self.order as usize {
            let start = coefficients.len();
            coefficients.resize(start + dim, 0.0);
            book.decode_vector_add(r, &mut coefficients[start..])?;
            for c in &mut coefficients[start..] {
                *c += last;
            }
            last = coefficients[start + dim - 1];
        }
        Ok(Some(Floor0Frame { amplitude, coefficients }))
    }

    /// The Bark map of 6.2.3 for a vector of `n` values.
    pub fn bark_map(&self, n: usize) -> Vec<i32> {
        let rate = self.rate as f64;
        let size = self.bark_map_size as f64;
        let scale = size / bark(0.5 * rate);
        (0..n)
            .map(|i| {
                let foobar = (bark(rate * i as f64 / (2.0 * n as f64)) * scale).floor();
                (foobar as i64).min(self.bark_map_size as i64 - 1) as i32
            })
            .collect()
    }

    /// `p + q` of 6.2.3 at angle `omega`: the squared magnitude response of
    /// the LSP polynomial.
    pub fn p_plus_q(&self, coefficients: &[f32], omega: f64) -> f64 {
        let order = self.order as usize;
        let cw = omega.cos();
        let term = |j: usize| {
            let d = (coefficients[j] as f64).cos() - cw;
            4.0 * d * d
        };
        let (mut p, mut q);
        if order % 2 == 1 {
            p = 1.0 - cw * cw;
            for j in 0..=(order.saturating_sub(3)) / 2 {
                if 2 * j + 1 < order {
                    p *= term(2 * j + 1);
                }
            }
            q = 0.25;
            for j in 0..=(order - 1) / 2 {
                q *= term(2 * j);
            }
        } else {
            p = (1.0 - cw) / 2.0;
            q = (1.0 + cw) / 2.0;
            for j in 0..order / 2 {
                p *= term(2 * j + 1);
                q *= term(2 * j);
            }
        }
        p + q
    }

    /// Curve computation (6.2.3) into `out` (length `n`, the map's length).
    pub fn synthesize(&self, frame: &Floor0Frame, map: &[i32], out: &mut [f32]) {
        let n = map.len();
        let max = ((1u64 << self.amplitude_bits) - 1) as f64;
        let offset = self.amplitude_offset as f64;
        let mut i = 0;
        while i < n {
            let omega = std::f64::consts::PI * map[i] as f64 / self.bark_map_size as f64;
            let pq = self.p_plus_q(&frame.coefficients, omega);
            let value = (0.11512925 * (frame.amplitude as f64 * offset / (max * pq.sqrt()) - offset)).exp() as f32;
            // A coefficient set whose response has a zero (p + q = 0) has no
            // finite floor there; silence it rather than emit infinity.
            let value = if value.is_finite() { value } else { 0.0 };
            let condition = map[i];
            while i < n && map[i] == condition {
                out[i] = value;
                i += 1;
            }
        }
    }
}

/// Why a floor 0 decode stopped.
pub(crate) enum Floor0Error {
    EndOfPacket,
    Undecodable(&'static str),
}

impl From<EndOfPacket> for Floor0Error {
    fn from(_: EndOfPacket) -> Self {
        Floor0Error::EndOfPacket
    }
}

// Index loops read plainest against the formulas they check.
#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Multiply polynomials in z^-1.
    fn mul(a: &[f64], b: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        out
    }

    /// |A(e^{i w})|^2 for A = (P + Q)/2 built from the LSP angles: P has
    /// the odd-indexed angles and a (1 - z^-1) factor (even order) or
    /// (1 - z^-2) (odd order), Q the even-indexed ones (and a (1 + z^-1)
    /// factor for even order). An independent route to 6.2.3's p + q.
    fn lpc_power(lsp: &[f64], w: f64) -> f64 {
        let order = lsp.len();
        let mut p = if order.is_multiple_of(2) { vec![1.0, -1.0] } else { vec![1.0, 0.0, -1.0] };
        let mut q = if order.is_multiple_of(2) { vec![1.0, 1.0] } else { vec![1.0] };
        for (j, &a) in lsp.iter().enumerate() {
            let f = [1.0, -2.0 * a.cos(), 1.0];
            if j % 2 == 1 { p = mul(&p, &f) } else { q = mul(&q, &f) }
        }
        let len = p.len().max(q.len());
        p.resize(len, 0.0);
        q.resize(len, 0.0);
        let (mut re, mut im) = (0.0, 0.0);
        for k in 0..len {
            let c = (p[k] + q[k]) / 2.0;
            re += c * (w * k as f64).cos();
            im -= c * (w * k as f64).sin();
        }
        re * re + im * im
    }

    #[test]
    fn p_plus_q_is_the_lpc_power_response() {
        for order in [1usize, 2, 3, 4, 7, 10, 16, 21] {
            let lsp: Vec<f64> = (0..order).map(|j| PI * (j as f64 + 0.7) / (order as f64 + 1.0)).collect();
            let coefficients: Vec<f32> = lsp.iter().map(|&v| v as f32).collect();
            let fl = Floor0 { order: order as u8, rate: 44100, bark_map_size: 256, amplitude_bits: 6, amplitude_offset: 100, books: vec![0] };
            let lsp32: Vec<f64> = coefficients.iter().map(|&v| v as f64).collect();
            for k in 0..50 {
                let w = PI * (k as f64 + 0.31) / 50.0;
                let want = lpc_power(&lsp32, w);
                let got = fl.p_plus_q(&coefficients, w);
                assert!((got - want).abs() <= 1e-9 * want.max(1e-12) + 1e-12, "order {order} w {w}: {got} vs {want}");
            }
        }
    }

    #[test]
    fn bark_map_is_monotonic_and_bounded() {
        let fl = Floor0 { order: 16, rate: 44100, bark_map_size: 256, amplitude_bits: 6, amplitude_offset: 100, books: vec![0] };
        for n in [128, 1024] {
            let map = fl.bark_map(n);
            assert_eq!(map[0], 0);
            assert!(map.windows(2).all(|w| w[0] <= w[1]));
            assert!(*map.last().unwrap() <= 255);
            assert!(*map.last().unwrap() >= 250);
        }
    }

    /// The curve is exp(.11512925 (A * off / ((2^bits - 1) sqrt(p + q)) - off)),
    /// constant over runs of equal map values.
    #[test]
    fn curve_follows_the_formula() {
        let fl = Floor0 { order: 4, rate: 22050, bark_map_size: 128, amplitude_bits: 6, amplitude_offset: 80, books: vec![0] };
        let frame = Floor0Frame { amplitude: 40, coefficients: vec![0.3, 0.9, 1.7, 2.6] };
        let map = fl.bark_map(256);
        let mut out = vec![0f32; 256];
        fl.synthesize(&frame, &map, &mut out);
        for i in 0..256 {
            let w = PI * map[i] as f64 / 128.0;
            let pq = lpc_power(&[0.3f32 as f64, 0.9f32 as f64, 1.7f32 as f64, 2.6f32 as f64], w);
            let want = (0.11512925 * (40.0 * 80.0 / (63.0 * pq.sqrt()) - 80.0)).exp();
            assert!((out[i] as f64 / want - 1.0).abs() < 1e-5, "bin {i}");
        }
    }
}
