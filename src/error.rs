//! The crate's one error type.

/// What went wrong. Every malformed input comes back as one of these; the
/// decoder never panics on bytes it is given.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The bitstream breaks the syntax or the semantics of the Vorbis I
    /// specification (or of RFC 3533 for the Ogg layer): a field out of
    /// range, a Huffman tree that is over- or under-specified, a header that
    /// ends early, a page whose CRC does not match.
    #[error("invalid Vorbis data: {0}")]
    Invalid(String),
    /// Valid Vorbis this crate does not implement, named.
    #[error("unsupported Vorbis feature: {0}")]
    Unsupported(String),
    /// A configuration the caller asked for that cannot be coded: a channel
    /// count, sample rate or quality outside what the encoder supports.
    #[error("invalid Vorbis encoder configuration: {0}")]
    Config(String),
    /// Reading from or writing to the caller's stream failed.
    #[error("I/O error: {0}")]
    Io(String),
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e.to_string())
    }
}

pub(crate) fn invalid(msg: impl Into<String>) -> Error {
    Error::Invalid(msg.into())
}

/// `std::result::Result` with this crate's [`Error`].
pub type Result<T> = std::result::Result<T, Error>;
