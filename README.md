# rivet-vorbis

[![CI](https://github.com/safewords/rivet-vorbis/actions/workflows/ci.yml/badge.svg)](https://github.com/safewords/rivet-vorbis/actions/workflows/ci.yml)

A **Vorbis I** decoder and encoder in Rust, with an **Ogg** (RFC 3533)
reader and writer: no C, no system libraries, no build script, nothing to
install on a build host. Written from the Vorbis I specification and
RFC 3533, not translated from any other implementation.

Written for the **[rivet](https://github.com/safewords/rivet)**
transcoder, where it replaces `lewton` on the decode side and adds Vorbis
encoding. Usable on its own by anything holding Ogg Vorbis files, or
Matroska / WebM Vorbis (Xiph-laced `CodecPrivate` plus raw packets), and
wanting PCM back — or holding PCM and wanting Vorbis.

Published as `rivet-vorbis`; **imported as `vorbis`** (`use vorbis::…`).
One dependency (`thiserror`), no features, no build script.

```toml
[dependencies]
vorbis = { package = "rivet-vorbis", git = "https://github.com/safewords/rivet-vorbis", branch = "develop" }
```

## What it decodes

All of Vorbis I:

| | |
|---|---|
| **Headers** | identification, comment (read leniently, as 4.2 allows), setup; Xiph lacing for Matroska `CodecPrivate` |
| **Codebooks** | ordered, sparse and plain length lists; lookup types 0, 1 and 2, `sequence_p`; single-entry books per the 2015 errata; over- and under-specified Huffman trees refused |
| **Floors** | type 0 (LSP) and type 1 (piecewise linear) |
| **Residues** | types 0, 1 and 2, every cascade depth |
| **Mappings** | any number of submaps and coupling steps (square polar) |
| **Blocks** | block sizes 64 to 8192, long/short window shapes from the window flags |
| **Ogg** | CRC-checked pages, resynchronisation, packets across pages, chained links, multiplexed streams (other codecs skipped), granule positions: leading samples before zero and trailing samples past the last page's position trimmed |

Output is `f32` at full scale ±1.0, unclipped (Vorbis can exceed full
scale), one vector per channel or interleaved, in Vorbis channel order
(specification 4.3.9) — for 5.1, FL C FR RL RR LFE. Truncated packets
decode as far as they go, as 2.1.8 intends; packets that cannot be decoded
(a mode that does not exist) are skipped. `Decoder::set_strict(true)` and
`decode_ogg_strict` turn every such leniency into an error instead.

## What it encodes

Vorbis I at 8 to 192 kHz, mono to 7.1 (Vorbis channel order), variable
bit rate from one quality value, −1 to 10. Ogg files from `encode_ogg` or
`OggWriter`, or raw packets plus `codec_private()` for Matroska.

- **Blocks**: 256 and 2048 samples, switched per block by a transient
  detector (a high-passed energy jump against the four 64-sample segments
  before it).
- **Psychoacoustics**: band energies on a Bark scale, tonality from each
  band's peak against its neighbourhood, masker-dependent offsets
  (14.5 + z dB tonal, 5.5 dB noise), a spreading function (−25 dB/Bark
  below, −10 dB/Bark above), the threshold in quiet; quality moves the
  mask, the dead zone and the low-pass.
- **Floor 1** fitted to the mask: 35 posts on long blocks, 20 on short,
  coded coarse-to-fine by geometric bisection, three posts per partition
  with a master book and three subclass books.
- **Residue 2** over all channels: eight classes by partition peak, up to
  three cascade passes, lattice VQ books of dimension 1 to 4.
- **Coupling**: square polar on the front pair, and on the rear / side
  pairs of 4.0 to 7.1. The LFE is low-passed at 150 Hz.
- **Codebooks** designed at start-up from modelled symbol statistics
  (Huffman over Laplacian probabilities whose scales were chosen by
  measuring this encoder's output) and sent in the setup header.


## Speed

On a Ryzen 9 9950X (Windows, a shared machine, best of three), in
multiples of real time; `cargo run --release --example vorbis_bench -- <pcm.raw>
[seconds] [runs] [file.ogg …]` measures it.

| | before | now |
|---|---|---|
| encode 60 s of a 16-bit stereo album track, q2 / q5 / q8 (one thread) | 282 / 270 / 260 | 286 / 268 / 261 |
| the same, threaded | — | 1305 / 1297 / 1190 |
| decode Xiph.Org `moog.ogg` (floor 0) | 230 | 534 |
| decode `combustion.ogg` | 222 | 660 |
| decode this encoder's q5 | 954 | 943–1054 |

Floor 0 curves take each coefficient's cosine once per frame instead of
once per bin group; the IMDCT's FFT runs on split real and imaginary
arrays with each stage's twiddles contiguous, its two shortest stages
fused (vectorised; x86-64 builds baseline and AVX2, chosen at run time;
aarch64 NEON); the bit reader peeks with one 8-byte load. The encoder
decides its block sequence first (from the transient detector alone) and
then codes the blocks one `encode` call completes in parallel
(`Encoder::set_threads`; 1 keeps the caller's thread only).

**Same output everywhere.** No fused multiply-add, and every rewritten
kernel performs its operations in the original order, so decoded PCM and
encoded packets are the same to the bit on every CPU, code path
(`force-scalar` compiles the run-time selection out; CI tests both) and
thread count — and the same as before this work, for all 26 decodable
Xiph.Org vectors and the encoder at q2, q5 and q8.

## How it is checked

- **Xiph.Org test vectors** (`tests/vectors.rs`): the 29 streams at
  <https://people.xiph.org/~xiphmont/test-vectors/vorbis/>, fetched by
  `tools/fetch-vectors.sh` and checked against the SHA-256 manifest in
  `tests/vectors/SHA256SUMS` (they are not committed). **No reference PCM is
  published for them**, and no other decoder is used here, so the check is
  what the specification fixes on its own: every one decodes; the decoded
  length equals the length the granule positions declare, end trimming and
  chained links included (28 of 29 exactly; `unused-mode-test.ogg` codes 34
  packets with a mode its setup does not define, which are discarded); the
  output is finite and at plausible levels; and the four floor-0 (LSP)
  vectors, encodings of one source, agree with each other at 11.6–15.5 dB
  waveform SNR (an unrelated vector: −0.3 dB).
- **Round trips** (`tests/roundtrip.rs`) through this crate's decoder in
  strict mode: every packet decodes, the headers parse back to the bytes
  written, the first page is the 58-byte identification page, the last page
  carries the input length and the decoded length equals it. Figures below.
- **Specification unit tests**: the bit-packing examples of 2.1.6–2.1.7;
  the Huffman example of 3.2.1 and the over/under-specified and single-entry
  rules; lookup types 1 and 2 with and without `sequence_p`; `float32_unpack`;
  the floor 1 neighbour, `render_point` and `render_line` definitions, an
  exhaustive check that every floor 1 value can be coded from every
  prediction, a hand-worked curve; floor 0's `p + q` against the power
  response of the LPC polynomial built from the same LSPs, and the whole
  curve formula; residues 0, 1 and 2 from hand-packed bitstreams, including
  interleave, cascades, "do not decode" and truncation; the inverse and
  forward MDCT against their O(N²) definitions (N = 64 to 8192), lapped
  reconstruction; the window formula and the Princen–Bradley condition for
  every pair of lapping block sizes; square polar coupling both ways,
  exhaustively over ±40; the Ogg CRC check value and packets across pages.
- **Malformed input** (`tests/fuzz.rs`): damaged files, random and damaged
  packets, damaged headers — errors, never a panic.

Round trips measured 2026-10-02 (`cargo test --release --test roundtrip --
--nocapture`). SNR is the worst channel's, over the whole signal; the
signal is harmonic tones with vibrato, low-passed noise at −24 dB, clicks
and a stretch of silence (the "tonal" rows leave out the noise). kb/s
counts the whole Ogg file.

| rate | channels | quality | SNR dB | kb/s |
|---|---|---|---|---|
| 44100 | 2 | −1 | 5.0 | 66.5 |
| 44100 | 2 | 0 | 5.8 | 79.3 |
| 44100 | 2 | 2 | 8.8 | 109.9 |
| 44100 | 2 | 4 | 12.1 | 141.3 |
| 44100 | 2 | 6 | 15.2 | 172.7 |
| 44100 | 2 | 8 | 17.6 | 202.9 |
| 44100 | 2 | 10 | 20.4 | 233.8 |
| 44100 tonal | 2 | −1 / 5 / 10 | 6.2 / 14.6 / 21.5 | 40.4 / 66.6 / 84.2 |
| 44100 tonal | 1 | 2 / 5 / 10 | 10.0 / 14.5 / 21.0 | 27.2 / 33.4 / 41.5 |
| 22050 | 1 | 3 | 10.0 | 50.9 |
| 32000 | 1 | 6 | 17.6 | 79.0 |
| 48000 | 1 | 0 | 7.1 | 43.9 |
| 22050 | 2 | 5 | 13.7 | 135.4 |
| 32000 | 2 | 2 | 10.4 | 112.8 |
| 48000 | 2 | 8 | 18.0 | 215.4 |
| 44100 | 5.1 | 5 | 12.1 | 437.3 |
| 48000 | 5.1 | 2 | 7.3 | 315.2 |
| 22050 | 5.1 | 8 | 15.0 | 447.9 |
| 44100 | 3.0 / 4.0 / 5.0 / 6.1 | 4 | 13.3 / 13.3 / 13.3 / 10.4 | 211.6 / 290.6 / 354.7 / 454.1 |
| 48000 | 7.1 | 4 | 10.3 | 546.1 |

Waveform SNR is a coarse measure of a perceptual codec: noise is coded only
to its own masking threshold, so a noisy signal comes back with a low SNR
at any rate. The encoder is a working psychoacoustic encoder, not a tuned
one; it has not been listened against other Vorbis encoders.

## Provenance and licensing

Written from the Vorbis I specification (Xiph.Org, the 2015 revision with
its errata) and RFC 3533; **no Vorbis or Ogg implementation's source was
read** — not libvorbis, libogg, Tremor, stb_vorbis, lewton, FFmpeg's or any
other — and no other implementation is used in the tests.
[docs/PROVENANCE.md](docs/PROVENANCE.md) records the sources, the places
where the specification needed reading, and how each was settled.

## Using it

```rust
// Ogg Vorbis file to PCM.
let decoded = vorbis::decode_ogg(&bytes)?;
// decoded.samples: one Vec<f32> per channel; decoded.sample_rate

// Matroska / WebM: the CodecPrivate, then each block's packet.
let mut dec = vorbis::Decoder::from_xiph_lacing(&codec_private)?;
for packet in packets {
    let planar = dec.decode(packet)?; // empty for the first packet
}

// Encoding: interleaved f32 in Vorbis channel order.
let config = vorbis::EncoderConfig { sample_rate: 48_000, channels: 2, quality: 5.0, comments: vec![] };
let file = vorbis::encode_ogg(&config, &pcm)?;
// ...or packet by packet, for another container:
let mut enc = vorbis::Encoder::new(config)?;
let mut packets = enc.encode(&pcm)?;
packets.extend(enc.finish()?);
let codec_private = enc.codec_private();
```

## License

Open Encoding Attribution License v1.0 — a source-available (not OSI open-source)
license, royalty-free, with a commercial-attribution requirement. See
[LICENSE.md](LICENSE.md) and [NOTICE](NOTICE).
