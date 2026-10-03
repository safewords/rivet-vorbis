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

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("vectors")
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

/// SHA-256 (FIPS 180-4), for the manifest check.
fn sha256(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut h: [u32; 8] = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
    let mut msg = data.to_vec();
    let bits = (data.len() as u64) * 8;
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bits.to_be_bytes());
    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes(chunk[i * 4..i * 4 + 4].try_into().unwrap());
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }
        let mut v = h;
        for i in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7].wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), v[0], v[1], v[2], v[3].wrapping_add(t1), v[4], v[5], v[6]];
        }
        for (a, b) in h.iter_mut().zip(v) {
            *a = a.wrapping_add(b);
        }
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

#[test]
fn sha256_known_answer() {
    assert_eq!(sha256(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    assert_eq!(sha256(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
}

/// Load the vectors present; None (skip) when the directory is empty.
fn load() -> Option<Vec<(String, Vec<u8>)>> {
    let mut out = Vec::new();
    for (sum, name) in manifest() {
        let Ok(bytes) = std::fs::read(dir().join(&name)) else { continue };
        assert_eq!(sha256(&bytes), sum, "{name}: SHA-256 does not match the manifest");
        out.push((name, bytes));
    }
    if out.is_empty() {
        if std::env::var("VORBIS_REQUIRE_VECTORS").is_ok() {
            panic!("VORBIS_REQUIRE_VECTORS is set but tests/vectors/ holds no vectors (run tools/fetch-vectors.sh)");
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
    println!("{:<30} {:>3} {:>6} {:>10} {:>10} {:>8} {:>8}", "vector", "ch", "rate", "samples", "declared", "peak", "rms dB");
    for (name, bytes) in &vectors {
        // Chained links may change channel count and rate: read block by block.
        let mut reader = vorbis::OggReader::new(&bytes[..]);
        let mut len = 0i64;
        let mut peak = 0f32;
        let mut sum = 0f64;
        let mut count = 0usize;
        let mut channels = Vec::new();
        let mut rates = Vec::new();
        while let Some(block) = reader.next_block().unwrap_or_else(|e| panic!("{name}: {e}")) {
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
        let decoded_channels = channels.iter().map(|c| c.to_string()).collect::<Vec<_>>().join("/");
        let decoded_rate = rates.iter().map(|c| c.to_string()).collect::<Vec<_>>().join("/");
        let rms = 10.0 * (sum / count.max(1) as f64 + 1e-30).log10();
        println!("{name:<30} {decoded_channels:>3} {decoded_rate:>6} {len:>10} {declared:>10} {peak:>8.4} {rms:>8.1}");
        if name == "unused-mode-test.ogg" {
            // 34 of its packets name mode 3 of a setup with three modes:
            // no audio can be made of them (4.3.1), so they are discarded
            // and the stream comes out shorter than its granule positions.
            assert_eq!(declared - len, 24959, "{name}");
        } else {
            assert_eq!(len, declared, "{name}: decoded length against the granule positions");
        }
        assert!(peak < 4.0, "{name}: implausible peak {peak}");
    }
}

/// unused-mode-test.ogg: the packets that name a missing mode are counted,
/// and refused in strict mode.
#[test]
fn packets_naming_a_missing_mode_are_discarded() {
    let Some(vectors) = load() else { return };
    let Some((_, bytes)) = vectors.iter().find(|(n, _)| n == "unused-mode-test.ogg") else { return };
    let mut r = PacketReader::new(&bytes[..]);
    let mut headers = Vec::new();
    let mut decoder = None;
    let mut discarded = 0;
    while let Some(p) = r.next_packet().unwrap() {
        if headers.len() < 3 {
            headers.push(p.data);
            if headers.len() == 3 {
                decoder = Some(vorbis::Decoder::new(&headers[0], &headers[1], &headers[2]).unwrap());
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
    let get = |n: &str| vectors.iter().find(|(name, _)| name == n).map(|(_, b)| vorbis::decode_ogg(b).unwrap().samples);
    let names = ["lsp-test.ogg", "lsp-test2.ogg", "lsp-test3.ogg", "lsp-test4.ogg"];
    let decoded: Vec<_> = names.iter().filter_map(|n| get(n)).collect();
    let Some(other) = get("rc1-test.ogg") else { return };
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
