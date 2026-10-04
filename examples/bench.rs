//! Encoder and decoder throughput:
//! `cargo run --release --example bench -- <pcm.raw> [seconds] [runs] [file.ogg …]`.
//!
//! `pcm.raw` is 16-bit little-endian stereo PCM at 44.1 kHz. Its first
//! `seconds` (default 60) are encoded at qualities 2, 5 and 8, and each
//! stream decoded back, each the best of `runs` (default 3). Each extra
//! Ogg Vorbis file named is decoded too. Every encoded stream's and decoded
//! output's FNV-1a hash is printed, so a change to either shows. Encoding
//! runs on one thread, then on all of them (which must give the same
//! stream).

use std::time::Instant;

use vorbis::{Encoder, EncoderConfig, OggWriter, decode_ogg};

fn best<T>(runs: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut out = None;
    let mut t = f64::MAX;
    for _ in 0..runs {
        let s = Instant::now();
        let v = f();
        t = t.min(s.elapsed().as_secs_f64());
        out = Some(v);
    }
    (t, out.expect("at least one run"))
}

fn fnv(bytes: impl IntoIterator<Item = u8>) -> u64 {
    bytes.into_iter().fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3))
}

/// `vorbis::encode_ogg` on `threads` threads.
fn encode_ogg(config: &EncoderConfig, pcm: &[f32], threads: usize) -> Vec<u8> {
    let mut encoder = Encoder::new(config.clone()).expect("encoder");
    encoder.set_threads(threads);
    let serial = 0x5256_4f42 ^ (pcm.len() as u32).rotate_left(7);
    let mut writer = OggWriter::new(Vec::new(), encoder.headers(), serial).expect("writer");
    writer.write_packets(encoder.encode(pcm).expect("encode")).expect("write");
    writer.write_packets(encoder.finish().expect("finish")).expect("write");
    writer.finish().expect("finish")
}

fn pcm_hash(samples: &[Vec<f32>]) -> u64 {
    fnv(samples.iter().flatten().flat_map(|s| s.to_bits().to_le_bytes()))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: bench <pcm.raw> [seconds] [runs] [file.ogg ...]");
    let seconds: f64 = args.get(2).map_or(60.0, |s| s.parse().expect("seconds"));
    let runs: usize = args.get(3).map_or(3, |s| s.parse().expect("runs"));
    let raw = std::fs::read(path).expect("read");
    let take = ((seconds * 44_100.0) as usize * 2).min(raw.len() / 2);
    let pcm: Vec<f32> =
        raw.as_chunks::<2>().0.iter().take(take).map(|b| f32::from(i16::from_le_bytes(*b)) / 32768.0).collect();
    let dur = pcm.len() as f64 / 2.0 / 44_100.0;
    println!("{dur:.1} s of stereo 44.1 kHz");
    for quality in [2.0f32, 5.0, 8.0] {
        let config = EncoderConfig { sample_rate: 44_100, channels: 2, quality, ..EncoderConfig::default() };
        let (te, stream) = best(runs, || encode_ogg(&config, &pcm, 1));
        let (tm, threaded) = best(runs, || encode_ogg(&config, &pcm, 0));
        assert!(threaded == stream, "the threaded encoder's stream differs");
        assert!(stream == vorbis::encode_ogg(&config, &pcm).expect("encode"), "encode_ogg's stream differs");
        let (td, out) = best(runs, || decode_ogg(&stream).expect("decode"));
        println!(
            "q{quality:<4} encode {:7.1} x realtime, {:7.1} x threaded ({} kB, hash {:016x}); decode {:7.1} x realtime (output hash {:016x})",
            dur / te,
            dur / tm,
            stream.len() / 1000,
            fnv(stream.iter().copied()),
            dur / td,
            pcm_hash(&out.samples)
        );
    }
    for file in args.iter().skip(4) {
        let bytes = std::fs::read(file).expect("read ogg");
        if let Err(e) = decode_ogg(&bytes) {
            println!("decode {file}: {e}");
            continue;
        }
        let (td, out) = best(runs, || decode_ogg(&bytes).expect("decode"));
        let d = out.samples.first().map_or(0, Vec::len) as f64 / f64::from(out.sample_rate);
        println!("decode {file}: {:7.1} x realtime ({d:.1} s, output hash {:016x})", d / td, pcm_hash(&out.samples));
    }
}
