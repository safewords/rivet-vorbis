//! Sample-wise comparison against Xiph.Org's reference decoder, `oggdec`
//! (vorbis-tools), run by `tools/make-references.sh` in Docker strictly as
//! a black-box binary. Its 16-bit output for the 29 test vectors and for
//! several files from this crate's encoder lives (uncommitted) in
//! `tests/references/`; `tests/references/SHA256SUMS` (committed) pins
//! every input and output.
//!
//! Opt-in: set `VORBIS_REFERENCE=1` after running the script. Without it
//! the comparison skips.
//!
//! The reference is 16-bit, so this decoder's float output is compared
//! after clipping to the 16-bit range: the error is reported in 16-bit
//! LSBs (maximum and RMS), with the share of samples whose own rounding to
//! 16 bits lands on the reference value.

#![allow(clippy::needless_range_loop)]

use std::collections::BTreeMap;
use std::path::PathBuf;

use vorbis::{EncoderConfig, OggReader, encode_ogg};

mod common;
use common::sha256;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("references")
}

/// A deterministic test signal (as the round-trip suite's): harmonic
/// tones, low-passed noise, clicks, a stretch of silence; an LFE channel
/// gets low tones.
fn signal(rate: u32, channels: usize, seconds: f64) -> Vec<f32> {
    let n = (rate as f64 * seconds) as usize;
    let lfe = match channels {
        6 => Some(5),
        7 => Some(6),
        8 => Some(7),
        _ => None,
    };
    let mut out = vec![0f32; n * channels];
    let mut seed = 0x1234_5678_u64;
    let mut lp = vec![0f64; channels];
    let tau = 2.0 * std::f64::consts::PI;
    for i in 0..n {
        let t = i as f64 / rate as f64;
        for c in 0..channels {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let white = ((seed >> 33) as f64 / (1u64 << 31) as f64) * 2.0 - 1.0;
            lp[c] = 0.9 * lp[c] + 0.1 * white;
            let f0 = 220.0 * (1.0 + c as f64 * 0.25) * (1.0 + 0.003 * (tau * 5.0 * t).sin());
            let tone: f64 = (1..8)
                .map(|h| (tau * f0 * h as f64 * t).sin() / h as f64)
                .sum();
            let mut v = 0.18 * tone * (0.5 + 0.5 * (tau * 0.7 * t).sin().abs()) + 0.06 * lp[c];
            if lfe == Some(c) {
                v = 0.3 * (tau * 50.0 * t).sin() + 0.2 * (tau * 80.0 * t).sin();
            }
            if t % 0.6 < 0.01 {
                v += 0.6 * (-(t % 0.6) * 400.0).exp() * white;
            }
            let f = i as f64 / n as f64;
            if f > 0.4 && f < 0.5 {
                v = 0.0;
            }
            out[i * channels + c] = v as f32;
        }
    }
    out
}

/// The encoder files the reference decodes: name, rate, channels, quality.
const FIXTURES: [(&str, u32, u8, f32); 6] = [
    ("enc-44100-2-q5", 44100, 2, 5.0),
    ("enc-44100-2-qm1", 44100, 2, -1.0),
    ("enc-48000-2-q10", 48000, 2, 10.0),
    ("enc-22050-1-q2", 22050, 1, 2.0),
    ("enc-48000-6-q4", 48000, 6, 4.0),
    ("enc-32000-8-q6", 32000, 8, 6.0),
];

fn fixture(rate: u32, channels: u8, quality: f32) -> Vec<u8> {
    // 2.37 s: not a whole number of blocks, so the end is trimmed.
    let input = signal(rate, channels as usize, 2.37);
    encode_ogg(
        &EncoderConfig {
            sample_rate: rate,
            channels,
            quality,
            comments: vec![("TITLE".into(), "reference".into())],
        },
        &input,
    )
    .unwrap()
}

/// With `VORBIS_WRITE_ENCODER_FIXTURES`, write the encoder files for the
/// script to decode; otherwise nothing.
#[test]
fn write_encoder_fixtures() {
    if std::env::var("VORBIS_WRITE_ENCODER_FIXTURES").is_err() {
        return;
    }
    std::fs::create_dir_all(dir()).unwrap();
    for (name, rate, ch, q) in FIXTURES {
        std::fs::write(dir().join(format!("{name}.ogg")), fixture(rate, ch, q)).unwrap();
    }
}

fn manifest() -> BTreeMap<String, String> {
    let Ok(text) = std::fs::read_to_string(dir().join("SHA256SUMS")) else {
        return BTreeMap::new();
    };
    text.lines()
        .filter_map(|l| l.split_once("  "))
        .map(|(sum, name)| (name.to_string(), sum.to_string()))
        .collect()
}

/// This decoder's output for a file, every link's channels interleaved and
/// concatenated (as oggdec writes a chained file), and the file's links'
/// channel counts.
fn decode_flat(bytes: &[u8]) -> Vec<f32> {
    let mut reader = OggReader::new(bytes);
    let mut out = Vec::new();
    while let Some(block) = reader.next_block().unwrap() {
        out.extend(vorbis::interleave(&block.samples));
    }
    out
}

struct Figures {
    ours: usize,
    reference: usize,
    max_lsb: f64,
    rms_lsb: f64,
    exact: f64,
}

fn compare(ours: &[f32], reference: &[i16]) -> Figures {
    let n = ours.len().min(reference.len());
    let mut max = 0f64;
    let mut sum = 0f64;
    let mut exact = 0usize;
    for i in 0..n {
        let x = (ours[i] as f64 * 32768.0).clamp(-32768.0, 32767.0);
        let r = reference[i] as f64;
        let e = (x - r).abs();
        max = max.max(e);
        sum += e * e;
        if x.round() == r {
            exact += 1;
        }
    }
    Figures {
        ours: ours.len(),
        reference: reference.len(),
        max_lsb: max,
        rms_lsb: (sum / n.max(1) as f64).sqrt(),
        exact: exact as f64 / n.max(1) as f64,
    }
}

#[test]
fn matches_the_reference_decoder() {
    if std::env::var("VORBIS_REFERENCE").is_err() {
        eprintln!(
            "skipping: set VORBIS_REFERENCE=1 after tools/make-references.sh to compare against oggdec"
        );
        return;
    }
    let sums = manifest();
    assert!(
        !sums.is_empty(),
        "tests/references/SHA256SUMS is missing: run tools/make-references.sh"
    );
    let fixtures: BTreeMap<&str, (u32, u8, f32)> = FIXTURES
        .iter()
        .map(|&(n, r, c, q)| (n, (r, c, q)))
        .collect();
    println!(
        "{:<30} {:>10} {:>10} {:>9} {:>9} {:>8}",
        "file", "samples", "oggdec", "max LSB", "RMS LSB", "exact"
    );
    let mut failures = Vec::new();
    for (name, sum) in &sums {
        let Some(stem) = name.strip_suffix(".raw") else {
            continue;
        };
        let raw = std::fs::read(dir().join(name))
            .unwrap_or_else(|_| panic!("{name} missing: run tools/make-references.sh"));
        assert_eq!(
            &sha256(&raw),
            sum,
            "{name}: SHA-256 differs from the manifest"
        );
        let ogg_name = format!("{stem}.ogg");
        let ogg = match fixtures.get(stem) {
            // Encoder files are regenerated: a changed encoder makes the
            // reference stale, which is an error to fix by rerunning the script.
            Some(&(rate, ch, q)) => fixture(rate, ch, q),
            None => std::fs::read(dir().join(&ogg_name)).unwrap(),
        };
        assert_eq!(
            Some(&sha256(&ogg)),
            sums.get(&ogg_name),
            "{ogg_name}: input differs from the one decoded (stale references: rerun tools/make-references.sh)"
        );
        let reference: Vec<i16> = raw
            .as_chunks::<2>()
            .0
            .iter()
            .map(|b| i16::from_le_bytes(*b))
            .collect();
        let ours = decode_flat(&ogg);
        let f = compare(&ours, &reference);
        println!(
            "{:<30} {:>10} {:>10} {:>9.3} {:>9.4} {:>7.3}%",
            ogg_name,
            f.ours,
            f.reference,
            f.max_lsb,
            f.rms_lsb,
            f.exact * 100.0
        );
        if f.ours != f.reference || f.max_lsb > 1.0 {
            failures.push(ogg_name);
        }
    }
    assert!(
        failures.is_empty(),
        "disagreements with oggdec: {failures:?}"
    );
}
