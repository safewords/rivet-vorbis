//! Vorbis I, both ways, written from the Vorbis I specification, with the
//! Ogg encapsulation of RFC 3533.
//!
//! - [`Decoder`]: one logical stream's audio packets to PCM — floor 0 and
//!   floor 1, residues 0, 1 and 2, any number of submaps and coupling
//!   steps, both block sizes and their window shapes. Container-agnostic:
//!   build it from the three header packets (or Matroska's Xiph-laced
//!   `CodecPrivate`) and feed it packets.
//! - [`OggReader`] / [`decode_ogg`]: Ogg Vorbis files, with the granule
//!   positions applied (leading samples before zero and trailing samples
//!   past the last page's position trimmed), chained links followed, other
//!   multiplexed streams skipped.
//! - [`Encoder`] / [`OggWriter`] / [`encode_ogg`]: a psychoacoustic VBR
//!   encoder — floor 1, residue 2, square polar coupling, 256/2048-sample
//!   block switching — and its Ogg pages.
//! - [`ogg`]: pages and packets for any codec.
//! - [`codebook`], [`floor0`], [`floor1`], [`residue`], [`header`],
//!   [`tables`]: the specification's building blocks, public for tools and
//!   tests.
//!
//! PCM is `f32` at full scale ±1.0 (the decoder does not clip), one vector
//! per channel or interleaved, in Vorbis channel order (specification
//! 4.3.9): mono; L R; L C R; FL FR RL RR; FL C FR RL RR; FL C FR RL RR LFE;
//! FL C FR SL SR RC LFE; FL C FR SL SR RL RR LFE.
//!
//! ```
//! // Encode a second of a 440 Hz tone, then decode it again.
//! let rate = 48_000;
//! let pcm: Vec<f32> = (0..rate * 2)
//!     .map(|i| ((i / 2) as f32 * 440.0 * std::f32::consts::TAU / rate as f32).sin() * 0.5)
//!     .collect();
//! let config = vorbis::EncoderConfig { sample_rate: rate as u32, channels: 2, quality: 4.0, comments: vec![] };
//! let file = vorbis::encode_ogg(&config, &pcm)?;
//! let decoded = vorbis::decode_ogg(&file)?;
//! assert_eq!(decoded.samples.len(), 2);
//! assert_eq!(decoded.samples[0].len(), rate);
//!
//! // Packet by packet, as from Matroska: the headers, then each packet.
//! let mut encoder = vorbis::Encoder::new(config)?;
//! let mut packets = encoder.encode(&pcm)?;
//! packets.extend(encoder.finish()?);
//! let mut decoder = vorbis::Decoder::from_xiph_lacing(&encoder.codec_private())?;
//! let mut samples = 0;
//! for p in &packets {
//!     samples += decoder.decode(&p.data)?[0].len();
//! }
//! // Without the container's end trimming the last block comes out whole.
//! assert!(samples >= rate);
//! # Ok::<(), vorbis::Error>(())
//! ```

mod bits;
pub mod codebook;
mod decode;
pub mod encode;
mod error;
pub mod floor0;
pub mod floor1;
pub mod header;
mod mdct;
pub mod ogg;
pub mod residue;
mod stream;
pub mod tables;
mod window;

pub use decode::{Decoder, decouple, interleave};
pub use encode::{EncodedPacket, Encoder, EncoderConfig, OggWriter, couple, encode_ogg};
pub use error::{Error, Result};
pub use header::{Comments, Identification, Setup, split_xiph_lacing, xiph_lacing};
pub use stream::{Block, Decoded, OggReader, decode_ogg, decode_ogg_strict};
