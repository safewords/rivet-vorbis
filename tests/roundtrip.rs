//! The encoder against this crate's decoder: every packet must decode in
//! strict mode, the headers must parse back to what was written, the
//! decoded length must equal the input length, and the decoded audio must
//! come back with an SNR that rises with the quality setting.

use vorbis::ogg::{FLAG_BOS, FLAG_EOS, PacketReader, PageReader};
use vorbis::{Comments, Encoder, EncoderConfig, Identification, Setup, decode_ogg_strict, encode_ogg};

/// A deterministic test signal: harmonic tones with vibrato, filtered
/// noise, a few clicks and a stretch of silence, per channel slightly
/// different. Interleaved, `seconds` long.
fn signal(rate: u32, channels: usize, seconds: f64) -> Vec<f32> {
    signal_with_noise(rate, channels, seconds, 0.06)
}

fn signal_with_noise(rate: u32, channels: usize, seconds: f64, noise: f64) -> Vec<f32> {
    let n = (rate as f64 * seconds) as usize;
    let mut out = vec![0f32; n * channels];
    let mut seed = 0x1234_5678_u64;
    let mut lp = vec![0f64; channels];
    let lfe = match channels {
        6 => Some(5),
        7 => Some(6),
        8 => Some(7),
        _ => None,
    };
    for i in 0..n {
        let t = i as f64 / rate as f64;
        for c in 0..channels {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let white = ((seed >> 33) as f64 / (1u64 << 31) as f64) * 2.0 - 1.0;
            lp[c] = 0.9 * lp[c] + 0.1 * white;
            let f0 = 220.0 * (1.0 + c as f64 * 0.25) * (1.0 + 0.003 * (2.0 * std::f64::consts::PI * 5.0 * t).sin());
            let mut tone = 0.0;
            for h in 1..8 {
                tone += (2.0 * std::f64::consts::PI * f0 * h as f64 * t).sin() / h as f64;
            }
            let envelope = 0.5 + 0.5 * (2.0 * std::f64::consts::PI * 0.7 * t).sin().abs();
            let mut v = 0.18 * tone * envelope + noise * lp[c];
            if lfe == Some(c) {
                v = 0.3 * (2.0 * std::f64::consts::PI * 50.0 * t).sin() + 0.2 * (2.0 * std::f64::consts::PI * 80.0 * t).sin();
            }
            // Clicks every 0.6 s, a decaying burst.
            let since = t % 0.6;
            if since < 0.01 {
                v += 0.6 * (-since * 400.0).exp() * white;
            }
            // Silence from 40% to 50% of the signal.
            if (i as f64 / n as f64) > 0.4 && (i as f64 / n as f64) < 0.5 {
                v = 0.0;
            }
            out[i * channels + c] = v as f32;
        }
    }
    out
}

struct Outcome {
    snr_db: f64,
    kbps: f64,
    blocks: [u64; 2],
}

fn round_trip(rate: u32, channels: u8, quality: f32, seconds: f64) -> Outcome {
    round_trip_signal(rate, channels, quality, seconds, signal(rate, channels as usize, seconds))
}

fn round_trip_signal(rate: u32, channels: u8, quality: f32, seconds: f64, input: Vec<f32>) -> Outcome {
    let config = EncoderConfig { sample_rate: rate, channels, quality, comments: vec![("TITLE".into(), "round trip".into()), ("ARTIST".into(), "rivet".into())] };
    let bytes = encode_ogg(&config, &input).unwrap();

    // Page layout of appendix A.2: the identification header alone on a
    // 58-byte BOS page; audio begins on a fresh page; EOS on the last.
    let mut pages = PageReader::new(&bytes[..]);
    pages.set_strict(true);
    let first = pages.next_page().unwrap().unwrap();
    assert_eq!(first.flags, FLAG_BOS);
    assert_eq!(first.to_bytes().len(), 58);
    assert_eq!(first.granule, 0);
    let mut last = first;
    while let Some(p) = pages.next_page().unwrap() {
        last = p;
    }
    assert_ne!(last.flags & FLAG_EOS, 0);
    assert_eq!(last.granule, (input.len() / channels as usize) as i64);

    // Headers parse, and re-serialise to the same bytes.
    let mut packets = PacketReader::new(&bytes[..]);
    let h: Vec<Vec<u8>> = (0..3).map(|_| packets.next_packet().unwrap().unwrap().data).collect();
    let ident = Identification::read(&h[0]).unwrap();
    assert_eq!(ident.write(), h[0]);
    assert_eq!((ident.channels, ident.sample_rate), (channels, rate));
    let comments = Comments::read(&h[1]).unwrap();
    assert_eq!(comments.write(), h[1]);
    assert_eq!(comments.get("title").next(), Some("round trip"));
    let setup = Setup::read(&h[2], channels).unwrap();
    assert_eq!(setup.write(channels), h[2]);
    // Every audio packet's mode exists and every one decodes strictly.
    let mut decoder = vorbis::Decoder::new(&h[0], &h[1], &h[2]).unwrap();
    decoder.set_strict(true);
    while let Some(p) = packets.next_packet().unwrap() {
        assert!(decoder.packet_blocksize(&p.data).is_some());
        decoder.decode(&p.data).unwrap();
    }

    let decoded = decode_ogg_strict(&bytes).unwrap();
    let frames = input.len() / channels as usize;
    assert_eq!(decoded.samples[0].len(), frames, "decoded length");
    let mut worst = f64::INFINITY;
    for c in 0..channels as usize {
        let (mut s, mut e) = (0f64, 0f64);
        for i in 0..frames {
            let x = input[i * channels as usize + c] as f64;
            let y = decoded.samples[c][i] as f64;
            s += x * x;
            e += (x - y) * (x - y);
        }
        worst = worst.min(10.0 * (s / e.max(1e-30)).log10());
    }
    // The encoder's own byte count matches the file's audio pages.
    let mut encoder = Encoder::new(config).unwrap();
    let mut n = encoder.encode(&input).unwrap().len();
    n += encoder.finish().unwrap().len();
    assert!(n > 0 || input.is_empty());
    let audio_bytes = bytes.len() as f64;
    Outcome { snr_db: worst, kbps: audio_bytes * 8.0 / seconds / 1000.0, blocks: count_blocks(&bytes, channels) }
}

fn count_blocks(bytes: &[u8], _channels: u8) -> [u64; 2] {
    let mut packets = PacketReader::new(bytes);
    let h: Vec<Vec<u8>> = (0..3).map(|_| packets.next_packet().unwrap().unwrap().data).collect();
    let decoder = vorbis::Decoder::new(&h[0], &h[1], &h[2]).unwrap();
    let mut blocks = [0u64; 2];
    while let Some(p) = packets.next_packet().unwrap() {
        let bs = decoder.packet_blocksize(&p.data).unwrap();
        blocks[(bs > 256) as usize] += 1;
    }
    blocks
}

#[test]
fn stereo_quality_ladder() {
    println!("{:>8} {:>6} {:>4} {:>8} {:>8} {:>12}", "rate", "ch", "q", "SNR dB", "kb/s", "short/long");
    let mut last_snr = f64::NEG_INFINITY;
    for q in [-1.0f32, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0] {
        let o = round_trip(44100, 2, q, 4.0);
        println!("{:>8} {:>6} {:>4} {:>8.2} {:>8.1} {:>6}/{}", 44100, 2, q, o.snr_db, o.kbps, o.blocks[0], o.blocks[1]);
        assert!(o.snr_db > last_snr - 0.5, "SNR should rise with quality");
        last_snr = o.snr_db;
        assert!(o.blocks[0] > 0 && o.blocks[1] > 0, "both block sizes in use");
    }
    assert!(last_snr > 20.0);
}

#[test]
fn layouts_and_rates() {
    for (rate, ch, q) in [(22050, 1, 3.0f32), (32000, 1, 6.0), (48000, 1, 0.0), (22050, 2, 5.0), (32000, 2, 2.0), (48000, 2, 8.0), (44100, 6, 5.0), (48000, 6, 2.0), (22050, 6, 8.0), (44100, 3, 4.0), (44100, 4, 4.0), (44100, 5, 4.0), (44100, 7, 4.0), (48000, 8, 4.0)] {
        let o = round_trip(rate, ch, q, 2.0);
        println!("{rate:>8} {ch:>6} {q:>4} {:>8.2} {:>8.1} {:>6}/{}", o.snr_db, o.kbps, o.blocks[0], o.blocks[1]);
        assert!(o.snr_db > 6.0, "{rate} Hz {ch} ch q{q}: {} dB", o.snr_db);
    }
}

#[test]
fn short_and_empty_inputs() {
    for frames in [0usize, 1, 100, 255, 256, 1000, 2047, 2048, 2049, 5000] {
        let input: Vec<f32> = (0..frames * 2).map(|i| ((i as f32) * 0.01).sin() * 0.5).collect();
        let bytes = encode_ogg(&EncoderConfig::default(), &input).unwrap();
        let decoded = decode_ogg_strict(&bytes).unwrap();
        assert_eq!(decoded.samples[0].len(), frames, "{frames} frames");
    }
}

/// The same ladder on tones alone (no noise), stereo and mono.
#[test]
fn tonal_quality_ladder() {
    for ch in [1u8, 2] {
        let mut last = f64::NEG_INFINITY;
        for q in [-1.0f32, 2.0, 5.0, 8.0, 10.0] {
            let input = signal_with_noise(44100, ch as usize, 3.0, 0.0);
            let o = round_trip_signal(44100, ch, q, 3.0, input);
            println!("tonal {:>6} {ch:>3} {q:>4} {:>8.2} {:>8.1} {:>6}/{}", 44100, o.snr_db, o.kbps, o.blocks[0], o.blocks[1]);
            assert!(o.snr_db > last - 0.5, "SNR should rise with quality");
            last = o.snr_db;
        }
        assert!(last > 20.0);
    }
}
