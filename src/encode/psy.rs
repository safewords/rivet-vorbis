//! The psychoacoustic model: how much noise each MDCT bin can hide.
//!
//! A simultaneous-masking model on the MDCT spectrum itself, from the
//! textbook ingredients: band energies on a Bark scale, a tonality
//! estimate per band (peak against neighbourhood), the masker-dependent
//! offset (about 14.5 + z dB below a tonal masker, 5.5 dB below a noise
//! masker), a spreading function falling 25 dB per Bark below a masker and
//! 10 dB per Bark above it, and the absolute threshold of hearing (Terhardt's
//! approximation) with full scale at 96 dB SPL.

/// The Bark scale (Zwicker and Terhardt's approximation).
pub(crate) fn bark(f: f64) -> f64 {
    13.0 * (0.00076 * f).atan() + 3.5 * (f / 7500.0).powi(2).atan()
}

/// The threshold in quiet, dB SPL (Terhardt).
pub(crate) fn ath_db(f: f64) -> f64 {
    let k = (f.max(20.0)) / 1000.0;
    3.64 * k.powf(-0.8) - 6.5 * (-0.6 * (k - 3.3).powi(2)).exp() + 1e-3 * k.powi(4)
}

/// The model for one block size at one sample rate.
pub(crate) struct Psy {
    half: usize,
    /// First bin of each band, and one past the last.
    band_start: Vec<usize>,
    /// Bark at each band's centre.
    band_bark: Vec<f64>,
    /// Threshold in quiet per bin, as power (full scale 1.0 = 96 dB).
    ath: Vec<f64>,
}

impl Psy {
    pub(crate) fn new(half: usize, rate: u32) -> Self {
        let bin_hz = rate as f64 / 2.0 / half as f64;
        let mut band_start = vec![0usize];
        let mut start_bark = bark(0.0);
        for k in 1..half {
            if bark(k as f64 * bin_hz) - start_bark >= 0.5 {
                band_start.push(k);
                start_bark = bark(k as f64 * bin_hz);
            }
        }
        band_start.push(half);
        let band_bark = band_start.windows(2).map(|w| bark((w[0] + w[1]) as f64 / 2.0 * bin_hz)).collect();
        let ath = (0..half)
            .map(|k| {
                let f = (k as f64 + 0.5) * bin_hz;
                10f64.powf((ath_db(f).min(110.0) - 96.0) / 10.0)
            })
            .collect();
        Psy { half, band_start, band_bark, ath }
    }

    /// The noise power each bin can carry unheard, scaled by `adjust_db`
    /// (positive allows more noise).
    pub(crate) fn mask(&self, x: &[f32], adjust_db: f64) -> Vec<f64> {
        let half = self.half;
        let power: Vec<f64> = x.iter().map(|&v| (v as f64) * (v as f64)).collect();
        let bands = self.band_bark.len();
        let mut energy = vec![0f64; bands];
        let mut offset = vec![0f64; bands];
        for b in 0..bands {
            let (s, e) = (self.band_start[b], self.band_start[b + 1]);
            energy[b] = power[s..e].iter().sum();
            // Tonality: the band's strongest bin against the mean power of
            // a neighbourhood around it.
            let (peak_k, peak) = (s..e).map(|k| (k, power[k])).fold((s, 0.0), |a, b| if b.1 > a.1 { b } else { a });
            let lo = peak_k.saturating_sub(8);
            let hi = (peak_k + 9).min(half);
            let mean = power[lo..hi].iter().sum::<f64>() / (hi - lo) as f64;
            let ratio_db = 10.0 * ((peak + 1e-30) / (mean + 1e-30)).log10();
            let tonal = ((ratio_db - 5.0) / 7.0).clamp(0.0, 1.0);
            offset[b] = tonal * (14.5 + self.band_bark[b]) + (1.0 - tonal) * 5.5;
        }
        let mut threshold = vec![0f64; bands];
        for (b, t) in threshold.iter_mut().enumerate() {
            let mut sum = 0.0;
            for j in 0..bands {
                let dz = self.band_bark[b] - self.band_bark[j];
                let atten = if dz >= 0.0 { 10.0 * dz } else { -25.0 * dz };
                let db = atten + offset[j];
                if db < 120.0 {
                    sum += energy[j] * 10f64.powf(-db / 10.0);
                }
            }
            *t = sum;
        }
        let scale = 10f64.powf(adjust_db / 10.0);
        let mut out = vec![0f64; half];
        for (b, t) in threshold.iter().enumerate() {
            let (s, e) = (self.band_start[b], self.band_start[b + 1]);
            let per_bin = t / (e - s) as f64;
            for (o, a) in out[s..e].iter_mut().zip(&self.ath[s..e]) {
                *o = per_bin.max(*a) * scale;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn threshold_in_quiet_has_its_known_shape() {
        // Most sensitive near 3.3 kHz, rising steeply at both ends.
        assert!(ath_db(3300.0) < 0.0);
        assert!(ath_db(100.0) > 20.0);
        assert!(ath_db(16000.0) > 40.0);
    }

    /// A tone masks its neighbours more than far bins, and less noise is
    /// allowed under a tone than under noise of the same energy.
    #[test]
    fn masking_spreads_and_tones_mask_less() {
        let psy = Psy::new(1024, 44100);
        let mut tone = vec![0f32; 1024];
        tone[200] = 0.5;
        let m = psy.mask(&tone, 0.0);
        assert!(m[205] > m[400] * 10.0);
        assert!(m[195] > m[100]);
        let mut noise = vec![0f32; 1024];
        let mut s = 1u64;
        for v in noise[190..211].iter_mut() {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            *v = ((s >> 40) as f32 / (1u64 << 24) as f32 - 0.5) * 0.25;
        }
        let e_noise: f64 = noise.iter().map(|&v| (v as f64).powi(2)).sum();
        let scale = (0.25 / e_noise).sqrt() as f32;
        for v in noise.iter_mut() {
            *v *= scale;
        }
        let mn = psy.mask(&noise, 0.0);
        assert!(mn[200] > m[200] * 3.0, "noise {} tone {}", mn[200], m[200]);
    }
}
