//! The Vorbis window (4.3.1): the slope `sin(pi/2 sin^2((x + 1/2)/n pi/2))`
//! on each side, with a long block's sides shortened to the short slope
//! where it laps a short block.

use std::f64::consts::FRAC_PI_2;

/// The rising slope of length `len` (half of the lapping block size).
pub(crate) fn slope(len: usize) -> Vec<f32> {
    (0..len)
        .map(|i| {
            let s = ((i as f64 + 0.5) / len as f64 * FRAC_PI_2).sin();
            (FRAC_PI_2 * s * s).sin() as f32
        })
        .collect()
}

/// The slopes of one stream's two block sizes.
pub(crate) struct Windows {
    blocksize: [usize; 2],
    slopes: [Vec<f32>; 2],
}

impl Windows {
    pub(crate) fn new(blocksize: [usize; 2]) -> Self {
        Windows { blocksize, slopes: [slope(blocksize[0] / 2), slope(blocksize[1] / 2)] }
    }

    /// The window as 4.3.1 builds it, `[left_start, left_end)` rising and
    /// `[right_start, right_end)` falling: returns those bounds and the slope
    /// used on each side.
    fn layout(&self, long: bool, prev_long: bool, next_long: bool) -> (usize, usize, usize, usize, &[f32], &[f32]) {
        let n = self.blocksize[long as usize];
        let bs0 = self.blocksize[0];
        let (ls, le, lslope) = if long && !prev_long {
            (n / 4 - bs0 / 4, n / 4 + bs0 / 4, &self.slopes[0])
        } else {
            (0, n / 2, &self.slopes[long as usize])
        };
        let (rs, re, rslope) = if long && !next_long {
            (n * 3 / 4 - bs0 / 4, n * 3 / 4 + bs0 / 4, &self.slopes[0])
        } else {
            (n / 2, n, &self.slopes[long as usize])
        };
        (ls, le, rs, re, lslope, rslope)
    }

    /// Multiply `buf` (one block) by its window. The flags are those of a
    /// long block; for a short block they are ignored.
    pub(crate) fn apply(&self, buf: &mut [f32], long: bool, prev_long: bool, next_long: bool) {
        let (ls, le, rs, re, lslope, rslope) = self.layout(long, prev_long, next_long);
        let n = buf.len();
        buf[..ls].fill(0.0);
        for (b, w) in buf[ls..le].iter_mut().zip(lslope) {
            *b *= *w;
        }
        // [le, rs) is one.
        for (b, w) in buf[rs..re].iter_mut().zip(rslope.iter().rev()) {
            *b *= *w;
        }
        buf[re..n].fill(0.0);
    }

    /// The window itself, for inspection and tests.
    pub(crate) fn shape(&self, long: bool, prev_long: bool, next_long: bool) -> Vec<f32> {
        let mut w = vec![1f32; self.blocksize[long as usize]];
        self.apply(&mut w, long, prev_long, next_long);
        w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Each side follows the formula of 4.3.1 literally.
    #[test]
    fn windows_follow_the_spec_formula() {
        let windows = Windows::new([256, 2048]);
        for (long, prev, next) in [(false, false, false), (true, true, true), (true, false, true), (true, true, false), (true, false, false)] {
            let n = if long { 2048 } else { 256 };
            let w = windows.shape(long, prev, next);
            let (ls, le, ln) = if long && !prev { (n / 4 - 64, n / 4 + 64, 128) } else { (0, n / 2, n / 2) };
            let (rs, re, rn) = if long && !next { (n * 3 / 4 - 64, n * 3 / 4 + 64, 128) } else { (n / 2, n, n / 2) };
            for i in 0..n {
                let want = if i < ls {
                    0.0
                } else if i < le {
                    let s = ((i - ls) as f64 + 0.5) / ln as f64 * PI / 2.0;
                    (PI / 2.0 * s.sin().powi(2)).sin()
                } else if i < rs {
                    1.0
                } else if i < re {
                    let s = ((i - rs) as f64 + 0.5) / rn as f64 * PI / 2.0 + PI / 2.0;
                    (PI / 2.0 * s.sin().powi(2)).sin()
                } else {
                    0.0
                };
                assert!((w[i] as f64 - want).abs() < 1e-6, "({long},{prev},{next}) i {i}: {} vs {want}", w[i]);
            }
        }
    }

    /// Princen-Bradley: where two windows lap, w^2 + w'^2 = 1, for every
    /// pair of block sizes the stream can lap.
    #[test]
    fn lapping_windows_are_power_complementary() {
        let windows = Windows::new([256, 2048]);
        // (previous block, its next flag) lapped with (current block, its prev flag).
        let cases = [(false, false), (true, true), (true, false), (false, true)];
        for (prev_long, cur_long) in cases {
            let np = if prev_long { 2048 } else { 256 };
            let nc = if cur_long { 2048 } else { 256 };
            let wp = windows.shape(prev_long, true, cur_long);
            let wc = windows.shape(cur_long, prev_long, true);
            // The current block's start sits at the previous block's 3/4
            // point less nc/4; the finished range runs from the previous
            // block's centre to the current one's.
            let offset = (np * 3 / 4) as i64 - (nc / 4) as i64;
            for i in np / 2..np * 3 / 4 + nc / 4 {
                let a = if i < np { wp[i] as f64 } else { 0.0 };
                let j = i as i64 - offset;
                let b = if j >= 0 && (j as usize) < nc { wc[j as usize] as f64 } else { 0.0 };
                assert!((a * a + b * b - 1.0).abs() < 1e-6, "{prev_long}->{cur_long} at {i}: {a} {b}");
            }
        }
    }
}
