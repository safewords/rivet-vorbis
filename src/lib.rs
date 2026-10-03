//! Vorbis I, both ways, written from the Vorbis I specification, with the
//! Ogg encapsulation of RFC 3533.

mod bits;
pub mod codebook;
mod error;
mod mdct;
mod window;

pub use error::{Error, Result};
