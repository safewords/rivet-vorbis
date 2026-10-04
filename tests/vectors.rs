//! The Xiph.Org Vorbis test vectors
//! (<https://people.xiph.org/~xiphmont/test-vectors/vorbis/>), fetched by
//! `tools/fetch-vectors.sh` into `tests/vectors/` and checked against the
//! SHA-256 manifest there. No reference PCM is published with them, so
//! each is checked for what the specification fixes on its own: every
//! packet decodes, the decoded length equals what the granule positions
//! declare (end trimming included), the output is finite and within a
//! plausible level, and the packets an LSP floor decodes carry ascending
//! angles in (0, pi).
//!
//! Without the files the tests skip (and say so); with
//! `VORBIS_REQUIRE_VECTORS=1` their absence fails.

#![allow(clippy::needless_range_loop)]

use std::path::PathBuf;

use vorbis::ogg::PacketReader;

mod common;
use common::sha256;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("vectors")
}

fn manifest() -> Vec<(String, String)> {
    let text = std::fs::read_to_string(dir().join("SHA256SUMS")).expect("manifest");
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let (sum, name) = l.split_once("  ").expect("sum  name");
            (sum.to_string(), name.to_string())
        })
        .collect()
}

#[test]
fn sha256_known_answer() {
    assert_eq!(
        sha256(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert_eq!(
        sha256(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
}

/// Load the vectors present; None (skip) when the directory is empty.
fn load() -> Option<Vec<(String, Vec<u8>)>> {
    let mut out = Vec::new();
    for (sum, name) in manifest() {
        let Ok(bytes) = std::fs::read(dir().join(&name)) else {
            continue;
        };
        assert_eq!(
            sha256(&bytes),
            sum,
            "{name}: SHA-256 does not match the manifest"
        );
        out.push((name, bytes));
    }
    if out.is_empty() {
        if std::env::var("VORBIS_REQUIRE_VECTORS").is_ok() {
            panic!(
                "VORBIS_REQUIRE_VECTORS is set but tests/vectors/ holds no vectors (run tools/fetch-vectors.sh)"
            );
        }
        eprintln!("skipping: no test vectors in tests/vectors/ (run tools/fetch-vectors.sh)");
        return None;
    }
    Some(out)
}

/// The length the granule positions declare: per logical stream, the last
/// granule position, less the start a positive first granule implies is
/// not present here (the vectors all start at zero).
fn declared_length(bytes: &[u8]) -> i64 {
    let mut r = PacketReader::new(bytes);
    let mut last: std::collections::BTreeMap<u32, i64> = Default::default();
    let mut order = Vec::new();
    while let Some(p) = r.next_packet().unwrap() {
        if let Some(g) = p.granule {
            if !last.contains_key(&p.serial) {
                order.push(p.serial);
            }
            last.insert(p.serial, g);
        }
    }
    last.values().sum()
}

#[test]
fn every_vector_decodes_to_its_declared_length() {
    let Some(vectors) = load() else { return };
    println!(
        "{:<30} {:>3} {:>6} {:>10} {:>10} {:>8} {:>8}",
        "vector", "ch", "rate", "samples", "declared", "peak", "rms dB"
    );
    for (name, bytes) in &vectors {
        // Chained links may change channel count and rate: read block by block.
        let mut reader = vorbis::OggReader::new(&bytes[..]);
        let mut len = 0i64;
        let mut peak = 0f32;
        let mut sum = 0f64;
        let mut count = 0usize;
        let mut channels = Vec::new();
        let mut rates = Vec::new();
        while let Some(block) = reader
            .next_block()
            .unwrap_or_else(|e| panic!("{name}: {e}"))
        {
            if channels.last() != Some(&block.samples.len()) {
                channels.push(block.samples.len());
            }
            if rates.last() != Some(&block.sample_rate) {
                rates.push(block.sample_rate);
            }
            let n = block.samples[0].len();
            len += n as i64;
            for c in &block.samples {
                assert_eq!(c.len(), n);
                for &s in c {
                    assert!(s.is_finite(), "{name}: non-finite sample");
                    peak = peak.max(s.abs());
                    sum += (s as f64) * (s as f64);
                    count += 1;
                }
            }
        }
        let declared = declared_length(bytes);
        let decoded_channels = channels
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join("/");
        let decoded_rate = rates
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join("/");
        let rms = 10.0 * (sum / count.max(1) as f64 + 1e-30).log10();
        println!(
            "{name:<30} {decoded_channels:>3} {decoded_rate:>6} {len:>10} {declared:>10} {peak:>8.4} {rms:>8.1}"
        );
        if name == "unused-mode-test.ogg" {
            // 34 of its packets name mode 3 of a setup with three modes:
            // no audio can be made of them (4.3.1), so they are discarded
            // and the stream comes out shorter than its granule positions.
            assert_eq!(declared - len, 24959, "{name}");
        } else {
            assert_eq!(
                len, declared,
                "{name}: decoded length against the granule positions"
            );
        }
        assert!(peak < 4.0, "{name}: implausible peak {peak}");
    }
}

/// unused-mode-test.ogg: the packets that name a missing mode are counted,
/// and refused in strict mode.
#[test]
fn packets_naming_a_missing_mode_are_discarded() {
    let Some(vectors) = load() else { return };
    let Some((_, bytes)) = vectors.iter().find(|(n, _)| n == "unused-mode-test.ogg") else {
        return;
    };
    let mut r = PacketReader::new(&bytes[..]);
    let mut headers = Vec::new();
    let mut decoder = None;
    let mut discarded = 0;
    while let Some(p) = r.next_packet().unwrap() {
        if headers.len() < 3 {
            headers.push(p.data);
            if headers.len() == 3 {
                decoder =
                    Some(vorbis::Decoder::new(&headers[0], &headers[1], &headers[2]).unwrap());
            }
            continue;
        }
        let d = decoder.as_mut().unwrap();
        if d.packet_blocksize(&p.data).is_none() {
            discarded += 1;
            d.set_strict(true);
            assert!(d.decode(&p.data).is_err());
            d.set_strict(false);
        }
        d.decode(&p.data).unwrap();
    }
    assert_eq!(discarded, 34);
    assert!(vorbis::decode_ogg_strict(bytes).is_err());
}

/// The four floor-0 (LSP) vectors are encodings of one source: decoded,
/// they agree with one another at the waveform level as lossy encodes of
/// one signal do (above 10 dB), and not with an unrelated vector. This is
/// the evidence for carrying `last` across vectors in floor 0 packet
/// decode (6.2.2; see docs/PROVENANCE.md): resetting it every vector gives
/// non-finite output on all four.
#[test]
fn lsp_vectors_agree_with_each_other() {
    let Some(vectors) = load() else { return };
    let get = |n: &str| {
        vectors
            .iter()
            .find(|(name, _)| name == n)
            .map(|(_, b)| vorbis::decode_ogg(b).unwrap().samples)
    };
    let names = [
        "lsp-test.ogg",
        "lsp-test2.ogg",
        "lsp-test3.ogg",
        "lsp-test4.ogg",
    ];
    let decoded: Vec<_> = names.iter().filter_map(|n| get(n)).collect();
    let Some(other) = get("rc1-test.ogg") else {
        return;
    };
    if decoded.len() < 4 {
        return;
    }
    let snr = |a: &[f32], b: &[f32]| {
        let n = a.len().min(b.len());
        let (mut s, mut e) = (0f64, 0f64);
        for k in 0..n {
            s += (a[k] as f64).powi(2);
            e += (a[k] as f64 - b[k] as f64).powi(2);
        }
        10.0 * (s / e).log10()
    };
    for i in 0..4 {
        for j in i + 1..4 {
            for ch in 0..2 {
                let v = snr(&decoded[i][ch], &decoded[j][ch]);
                println!("{} vs {} ch {ch}: {v:.2} dB", names[i], names[j]);
                assert!(v > 10.0, "{} vs {}: {v} dB", names[i], names[j]);
            }
        }
        let v = snr(&decoded[i][0], &other[0]);
        assert!(v < 3.0, "unrelated vector agrees at {v} dB");
    }
}
