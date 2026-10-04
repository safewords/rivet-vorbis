//! The MDCT of Vorbis I (4.3.7), forward for the encoder and inverse for
//! the decoder, through a DCT-IV computed with a complex FFT of a quarter
//! of the block size.
//!
//! With `N` the block size and `M = N/2` coefficients,
//!
//! ```text
//! forward: X[k] = sum_{n<N}  x[n] cos(2 pi / N (n + 1/2 + N/4)(k + 1/2))
//! inverse: y[n] = sum_{k<M}  X[k] cos(2 pi / N (n + 1/2 + N/4)(k + 1/2))
//! ```
//!
//! both unnormalised. The decoder applies the inverse as defined; the
//! encoder scales its forward transform by `4/N`, which makes the windowed
//! lapped pair an identity (see `forward`).

use std::f64::consts::PI;

/// Transforms for one block size.
pub(crate) struct Mdct {
    n: usize,
    /// Pre-twiddle `exp(-i pi j / M)`, `j < M/2`.
    pre: Vec<(f64, f64)>,
    /// Post-twiddle `exp(-i pi (p + 1/4) / M)`, `p < M/2`.
    post: Vec<(f64, f64)>,
    fft: Fft,
}

impl Mdct {
    /// Block size `n`: a power of two, at least 16.
    pub(crate) fn new(n: usize) -> Self {
        assert!(n.is_power_of_two() && n >= 16);
        let m = n / 2;
        let q = m / 2;
        let pre = (0..q)
            .map(|j| {
                let a = -PI * j as f64 / m as f64;
                (a.cos(), a.sin())
            })
            .collect();
        let post = (0..q)
            .map(|p| {
                let a = -PI * (p as f64 + 0.25) / m as f64;
                (a.cos(), a.sin())
            })
            .collect();
        Mdct {
            n,
            pre,
            post,
            fft: Fft::new(q),
        }
    }

    /// DCT-IV of `u` (length M) into `out`, unnormalised:
    /// `C[k] = sum u[n] cos(pi/M (n + 1/2)(k + 1/2))`.
    fn dct4(&self, u: &[f64], out: &mut [f64]) {
        let m = self.n / 2;
        let q = m / 2;
        let (mut zr, mut zi) = (vec![0f64; q], vec![0f64; q]);
        for j in 0..q {
            let (re, im) = (u[2 * j], u[m - 1 - 2 * j]);
            let (c, s) = self.pre[j];
            zr[j] = re * c - im * s;
            zi[j] = re * s + im * c;
        }
        self.fft.forward(&mut zr, &mut zi);
        for p in 0..q {
            let (re, im) = (zr[p], zi[p]);
            let (c, s) = self.post[p];
            let (sr, si) = (re * c - im * s, re * s + im * c);
            out[2 * p] = sr;
            out[m - 1 - 2 * p] = -si;
        }
    }

    /// Inverse MDCT: `spectrum` (N/2 values) to `out` (N samples), as the
    /// definition above, unwindowed.
    pub(crate) fn inverse(&self, spectrum: &[f32], out: &mut [f32]) {
        let n = self.n;
        let m = n / 2;
        let h = m / 2;
        let u: Vec<f64> = spectrum.iter().map(|&v| v as f64).collect();
        let mut v = vec![0f64; m];
        self.dct4(&u, &mut v);
        // y = [v2, -reverse(v), -v1] for v = [v1 v2] (halves of M/2).
        for i in 0..h {
            out[i] = v[h + i] as f32;
        }
        for i in 0..m {
            out[h + i] = -v[m - 1 - i] as f32;
        }
        for i in 0..h {
            out[h + m + i] = -v[i] as f32;
        }
    }

    /// Forward MDCT of `input` (N samples, already windowed) to `out` (N/2
    /// values), scaled by `4/N` so that the inverse, windowed and lapped,
    /// restores the input.
    pub(crate) fn forward(&self, input: &[f32], out: &mut [f32]) {
        let n = self.n;
        let m = n / 2;
        let h = m / 2;
        // Fold x = [a b c d] (quarters) into (-c_r - d, a - b_r).
        let a = |i: usize| input[i] as f64;
        let b = |i: usize| input[h + i] as f64;
        let c = |i: usize| input[m + i] as f64;
        let d = |i: usize| input[m + h + i] as f64;
        let mut u = vec![0f64; m];
        for i in 0..h {
            u[i] = -c(h - 1 - i) - d(i);
            u[h + i] = a(i) - b(h - 1 - i);
        }
        let mut v = vec![0f64; m];
        self.dct4(&u, &mut v);
        let scale = 4.0 / n as f64;
        for (o, x) in out.iter_mut().zip(&v) {
            *o = (*x * scale) as f32;
        }
    }
}

/// Iterative radix-2 complex FFT, `X[p] = sum z[j] exp(-2 pi i j p / L)`,
/// on separate real and imaginary arrays. Each stage's twiddles are laid
/// out contiguously (the same values as `exp(-2 pi i k / L)` at the
/// stage's stride), so the butterflies of a group run side by side; every
/// butterfly is the same arithmetic in the same order, so the result does
/// not depend on how wide the vectors are.
struct Fft {
    len: usize,
    rev: Vec<u32>,
    /// Per stage (butterfly span 2, 4, …, len): its `span / 2` twiddles,
    /// real parts then imaginary parts.
    stages: Vec<(Vec<f64>, Vec<f64>)>,
}

impl Fft {
    fn new(len: usize) -> Self {
        let bits = len.trailing_zeros();
        let rev = (0..len as u32)
            .map(|i| {
                if bits == 0 {
                    0
                } else {
                    i.reverse_bits() >> (32 - bits)
                }
            })
            .collect();
        let twiddle: Vec<(f64, f64)> = (0..len / 2)
            .map(|k| {
                let a = -2.0 * PI * k as f64 / len as f64;
                (a.cos(), a.sin())
            })
            .collect();
        let mut stages = Vec::new();
        let mut size = 2;
        while size <= len {
            let step = len / size;
            let (wr, wi) = (0..size / 2).map(|k| twiddle[k * step]).unzip();
            stages.push((wr, wi));
            size *= 2;
        }
        Fft { len, rev, stages }
    }

    fn forward(&self, re: &mut [f64], im: &mut [f64]) {
        let n = self.len;
        for i in 0..n {
            let j = self.rev[i] as usize;
            if j > i {
                re.swap(i, j);
                im.swap(i, j);
            }
        }
        let mut stages = self.stages.iter();
        if let [s1, s2, ..] = self.stages.as_slice() {
            // The two shortest spans, fused (the same butterflies in the
            // same order, four points at a time).
            first_stages(re, im, s1, s2);
            stages.nth(1);
        }
        for (wr, wi) in stages {
            butterflies(re, im, wr, wi);
        }
    }
}

crate::simd::multiversion! {
/// The stages of span 2 and 4 on each group of four points.
fn first_stages(re: &mut [f64], im: &mut [f64], s1: &(Vec<f64>, Vec<f64>), s2: &(Vec<f64>, Vec<f64>)) {
    let (w1r, w1i) = (s1.0[0], s1.1[0]);
    let (w2r, w2i): (&[f64; 2], &[f64; 2]) =
        (s2.0[..2].try_into().expect("two twiddles"), s2.1[..2].try_into().expect("two twiddles"));
    for (r, i) in re.as_chunks_mut::<4>().0.iter_mut().zip(im.as_chunks_mut::<4>().0) {
        // Span 2: (0, 1) and (2, 3).
        for p in [0, 2] {
            let (a, b) = (r[p + 1], i[p + 1]);
            let (tr, ti) = (a * w1r - b * w1i, a * w1i + b * w1r);
            let (cr, ci) = (r[p], i[p]);
            r[p] = cr + tr;
            i[p] = ci + ti;
            r[p + 1] = cr - tr;
            i[p + 1] = ci - ti;
        }
        // Span 4: (0, 2) and (1, 3).
        for k in 0..2 {
            let (a, b) = (r[k + 2], i[k + 2]);
            let (tr, ti) = (a * w2r[k] - b * w2i[k], a * w2i[k] + b * w2r[k]);
            let (cr, ci) = (r[k], i[k]);
            r[k] = cr + tr;
            i[k] = ci + ti;
            r[k + 2] = cr - tr;
            i[k + 2] = ci - ti;
        }
    }
}
}

crate::simd::multiversion! {
/// One radix-2 stage: for every group of `2 · half` points (`half` the
/// twiddle count), `t = w_k · z[k + half]`, `z[k] ± t`.
fn butterflies(re: &mut [f64], im: &mut [f64], wr: &[f64], wi: &[f64]) {
    let half = wr.len();
    for (gr, gi) in re.chunks_exact_mut(2 * half).zip(im.chunks_exact_mut(2 * half)) {
        let (ur, xr) = gr.split_at_mut(half);
        let (ui, xi) = gi.split_at_mut(half);
        for k in 0..half {
            let (a, b) = (xr[k], xi[k]);
            let (tr, ti) = (a * wr[k] - b * wi[k], a * wi[k] + b * wr[k]);
            let (cr, ci) = (ur[k], ui[k]);
            ur[k] = cr + tr;
            ui[k] = ci + ti;
            xr[k] = cr - tr;
            xi[k] = ci - ti;
        }
    }
}
}

// Index loops read plainest against the formulas they check.
#[cfg(test)]
#[allow(clippy::needless_range_loop)]
mod tests {
    use super::*;

    fn phase(n: usize, i: usize, k: usize) -> f64 {
        2.0 * PI / n as f64 * (i as f64 + 0.5 + n as f64 / 4.0) * (k as f64 + 0.5)
    }

    fn noise(len: usize, seed: u64) -> Vec<f32> {
        let mut s = seed;
        (0..len)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((s >> 33) as f64 / (1u64 << 31) as f64 * 2.0 - 1.0) as f32
            })
            .collect()
    }

    /// The FFT as first written: complex pairs, twiddles indexed at the
    /// stage's stride, one butterfly at a time.
    fn fft_pairs(z: &mut [(f64, f64)]) {
        let n = z.len();
        let bits = n.trailing_zeros();
        for i in 0..n {
            let j = if bits == 0 {
                0
            } else {
                (i as u32).reverse_bits() as usize >> (32 - bits)
            };
            if j > i {
                z.swap(i, j);
            }
        }
        let twiddle: Vec<(f64, f64)> = (0..n / 2)
            .map(|k| {
                let a = -2.0 * PI * k as f64 / n as f64;
                (a.cos(), a.sin())
            })
            .collect();
        let mut size = 2;
        while size <= n {
            let half = size / 2;
            let step = n / size;
            for start in (0..n).step_by(size) {
                for k in 0..half {
                    let (wr, wi) = twiddle[k * step];
                    let (xr, xi) = z[start + k + half];
                    let t = (xr * wr - xi * wi, xr * wi + xi * wr);
                    let u = z[start + k];
                    z[start + k] = (u.0 + t.0, u.1 + t.1);
                    z[start + k + half] = (u.0 - t.0, u.1 - t.1);
                }
            }
            size *= 2;
        }
    }

    #[test]
    fn the_split_fft_is_the_paired_one_to_the_bit() {
        for len in [1usize, 2, 4, 8, 16, 64, 512, 2048] {
            let x = noise(2 * len, len as u64);
            let mut pairs: Vec<(f64, f64)> = (0..len)
                .map(|i| (x[2 * i] as f64, x[2 * i + 1] as f64))
                .collect();
            let (mut re, mut im): (Vec<f64>, Vec<f64>) = pairs.iter().copied().unzip();
            fft_pairs(&mut pairs);
            Fft::new(len).forward(&mut re, &mut im);
            for (i, &(a, b)) in pairs.iter().enumerate() {
                assert_eq!(
                    (re[i].to_bits(), im[i].to_bits()),
                    (a.to_bits(), b.to_bits()),
                    "len {len} bin {i}"
                );
            }
        }
    }

    /// The fast inverse against the definition, summed in f64.
    #[test]
    fn inverse_matches_the_definition() {
        for n in [64usize, 128, 256, 512, 2048, 8192] {
            let mdct = Mdct::new(n);
            let x = noise(n / 2, n as u64);
            let mut fast = vec![0f32; n];
            mdct.inverse(&x, &mut fast);
            let mut worst = 0f64;
            for i in 0..n {
                let reference: f64 = (0..n / 2).map(|k| x[k] as f64 * phase(n, i, k).cos()).sum();
                worst = worst.max((reference - fast[i] as f64).abs());
            }
            // Relative to the sqrt(N/2) scale of a sum of N/2 unit terms.
            assert!(worst < 1e-5 * (n as f64).sqrt(), "n {n}: worst {worst}");
        }
    }

    #[test]
    fn forward_matches_the_definition() {
        for n in [64usize, 256, 2048] {
            let mdct = Mdct::new(n);
            let x = noise(n, 7 + n as u64);
            let mut fast = vec![0f32; n / 2];
            mdct.forward(&x, &mut fast);
            for k in 0..n / 2 {
                let reference: f64 = (0..n)
                    .map(|i| x[i] as f64 * phase(n, i, k).cos())
                    .sum::<f64>()
                    * 4.0
                    / n as f64;
                assert!((reference - fast[k] as f64).abs() < 1e-5, "n {n} k {k}");
            }
        }
    }

    /// Windowed with a power-complementary window and lapped by half a
    /// block, forward then inverse restores the signal (TDAC).
    #[test]
    fn lapped_transform_reconstructs() {
        let n = 256;
        let mdct = Mdct::new(n);
        let w: Vec<f32> = (0..n)
            .map(|i| {
                let s = ((i as f64 + 0.5) / n as f64 * PI).sin();
                (0.5 * PI * s * s).sin() as f32
            })
            .collect();
        let signal = noise(n * 4, 3);
        let mut out = vec![0f32; n * 4];
        let mut start = 0;
        while start + n <= signal.len() {
            let block: Vec<f32> = (0..n).map(|i| signal[start + i] * w[i]).collect();
            let mut spec = vec![0f32; n / 2];
            mdct.forward(&block, &mut spec);
            let mut y = vec![0f32; n];
            mdct.inverse(&spec, &mut y);
            for i in 0..n {
                out[start + i] += y[i] * w[i];
            }
            start += n / 2;
        }
        for i in n / 2..signal.len() - n / 2 {
            assert!(
                (out[i] - signal[i]).abs() < 1e-5,
                "sample {i}: {} vs {}",
                out[i],
                signal[i]
            );
        }
    }
}
