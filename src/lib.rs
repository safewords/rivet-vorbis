//! Vorbis I, both ways, written from the Vorbis I specification, with the
//! Ogg encapsulation of RFC 3533.

mod bits;
pub mod codebook;
mod decode;
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
pub use error::{Error, Result};
pub use header::{Comments, Identification, Setup, split_xiph_lacing, xiph_lacing};
pub use stream::{Block, Decoded, OggReader, decode_ogg, decode_ogg_strict};
