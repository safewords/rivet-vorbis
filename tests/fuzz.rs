//! Malformed input: valid streams with bits flipped, bytes cut and garbage
//! spliced, random audio packets, and damaged setup headers, through every
//! entry point. Errors are fine; a panic is not. Deterministic (a fixed
//! xorshift generator), no dependencies.

use vorbis::ogg::PacketReader;
use vorbis::{Decoder, EncoderConfig, decode_ogg, decode_ogg_strict, encode_ogg};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n.max(1) as u64) as usize
    }
}

fn sample_stream(channels: u8) -> Vec<u8> {
    let rate = 22050;
    let input: Vec<f32> = (0..rate as usize * channels as usize / 2)
        .map(|i| {
            let t = (i / channels as usize) as f32 / rate as f32;
            (t * 2.0 * std::f32::consts::PI * 440.0).sin() * 0.3 + if i % 5000 < 20 { 0.5 } else { 0.0 }
        })
        .collect();
    encode_ogg(&EncoderConfig { sample_rate: rate, channels, quality: 3.0, comments: vec![] }, &input).unwrap()
}

fn mutate(rng: &mut Rng, bytes: &[u8]) -> Vec<u8> {
    let mut b = bytes.to_vec();
    match rng.below(4) {
        0 => {
            for _ in 0..1 + rng.below(20) {
                let i = rng.below(b.len());
                b[i] ^= 1 << rng.below(8);
            }
        }
        1 => b.truncate(rng.below(b.len())),
        2 => {
            let at = rng.below(b.len());
            let junk: Vec<u8> = (0..rng.below(300)).map(|_| rng.next() as u8).collect();
            b.splice(at..at, junk);
        }
        _ => {
            let at = rng.below(b.len());
            let len = rng.below(b.len() - at);
            b.drain(at..at + len);
        }
    }
    b
}

#[test]
fn damaged_files_never_panic() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut streams = vec![sample_stream(1), sample_stream(2), sample_stream(6)];
    // A floor 0 stream too, when the test vectors are on disk.
    let vectors = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/vectors");
    for name in ["test-short.ogg", "lsp-test.ogg"] {
        if let Ok(b) = std::fs::read(vectors.join(name)) {
            streams.push(b[..b.len().min(60_000)].to_vec());
        }
    }
    for stream in &streams {
        assert!(decode_ogg(stream).is_ok());
        for _ in 0..150 {
            let damaged = mutate(&mut rng, stream);
            let _ = decode_ogg(&damaged);
            let _ = decode_ogg_strict(&damaged);
        }
    }
}

fn headers(stream: &[u8]) -> Vec<Vec<u8>> {
    let mut r = PacketReader::new(stream);
    (0..3).map(|_| r.next_packet().unwrap().unwrap().data).collect()
}

fn audio(stream: &[u8]) -> Vec<Vec<u8>> {
    let mut r = PacketReader::new(stream);
    let mut out = Vec::new();
    while let Some(p) = r.next_packet().unwrap() {
        out.push(p.data);
    }
    out.split_off(3)
}

#[test]
fn random_and_damaged_packets_never_panic() {
    let mut rng = Rng(42);
    for channels in [1u8, 2, 6] {
        let stream = sample_stream(channels);
        let h = headers(&stream);
        let packets = audio(&stream);
        let mut d = Decoder::new(&h[0], &h[1], &h[2]).unwrap();
        for _ in 0..1500 {
            let p: Vec<u8> = match rng.below(3) {
                0 => (0..rng.below(400)).map(|_| rng.next() as u8).collect(),
                1 => {
                    let i = rng.below(packets.len());
                    mutate(&mut rng, &packets[i])
                }
                _ => packets[rng.below(packets.len())].clone(),
            };
            d.set_strict(rng.below(2) == 0);
            if let Ok(out) = d.decode(&p) {
                assert_eq!(out.len(), channels as usize);
                assert!(out.iter().all(|c| c.iter().all(|s| s.is_finite())));
            }
        }
    }
}

#[test]
fn damaged_headers_never_panic() {
    let mut rng = Rng(7);
    let stream = sample_stream(2);
    let h = headers(&stream);
    let packets = audio(&stream);
    for _ in 0..600 {
        let which = rng.below(3);
        let mut hs = h.clone();
        hs[which] = mutate(&mut rng, &h[which]);
        if let Ok(mut d) = Decoder::new(&hs[0], &hs[1], &hs[2]) {
            for p in packets.iter().take(6) {
                let _ = d.decode(p);
            }
        }
        let _ = Decoder::from_xiph_lacing(&mutate(&mut rng, &vorbis::xiph_lacing([&h[0], &h[1], &h[2]])));
    }
}
