# Provenance

Where every part of this crate came from. The short version: the code is
this repository's own, written from the Vorbis I specification and RFC 3533;
the one normative table was copied from the specification as data; no other
Vorbis or Ogg implementation was read, run or compared against.

## Clean-room rules

- **No Vorbis or Ogg implementation's source was opened or read**, in
  writing either half or the tests: not libvorbis, libogg, Tremor,
  stb_vorbis, lewton, FFmpeg's (libavcodec `vorbis*`, libavformat `ogg*`),
  symphonia or any other. No implementation was searched for.
- **No other implementation is used in the tests**, not even as a black
  box: no ffmpeg, no oggdec, no C. The decoder is checked against the
  specification's own definitions (unit tests), against what the test
  vectors' own granule positions declare, and against this crate's encoder;
  the encoder against this crate's decoder.
- rivet's use of `lewton` (`crates/codec/src/audio/decode/vorbis.rs` in the
  rivet repository) was read only for its call shape — headers from a
  Xiph-laced buffer, then one packet at a time to planar `f32` — so that
  `Decoder::from_xiph_lacing` and `Decoder::decode` can replace it.

## Sources

**The Vorbis I specification**, Xiph.Org Foundation, the LaTeX edition of
2015 with its errata (<https://xiph.org/vorbis/doc/Vorbis_I_spec.html>,
fetched 2026-10-02), used throughout:

- Section 2 (bit packing): `bits.rs`; its coding and decoding examples are
  unit tests.
- Section 3 (codebooks): `codebook.rs` — the packed format, the codeword
  assignment rule ("the lowest valued unused binary Huffman codeword
  possible", which may lie left of a shorter codeword already assigned),
  over- and under-specified trees, the 2015-02-26 errata on single-entry
  codebooks, the VQ lookup of types 1 and 2.
- Section 4 (headers and audio packet decode): `header.rs`, `decode.rs`,
  `window.rs`.
- Section 5 (comments): `header.rs`.
- Section 6 (floor 0), with the 2015-02-27 errata on the Bark formula:
  `floor0.rs`.
- Section 7 (floor 1): `floor1.rs`.
- Section 8 (residues): `residue.rs`.
- Section 9 (`ilog`, `float32_unpack`, `lookup1_values`, `low_neighbor`,
  `high_neighbor`, `render_point`, `render_line`): where each is used.
- Section 10 (`floor1_inverse_dB_table`): `tables.rs`, the 256 values as
  the specification prints them.
- Appendix A (Ogg encapsulation, granule positions, start and end
  trimming): `stream.rs`, and the encoder's `OggWriter`.

**RFC 3533**, The Ogg Encapsulation Format Version 0: `ogg.rs` — the page
header, lacing, the CRC (polynomial 0x04c11db7, initial value 0, no
reflection, over the page with its CRC field zeroed; checked against the
standard check value of that CRC for "123456789").

**The MDCT** (4.3.7 points to Sporer, Brandenburg and Edler): the transform
is computed through the textbook identity between an MDCT and a DCT-IV of
half the size, and a DCT-IV through a complex FFT of a quarter of the size;
the derivation is this crate's and is checked against the O(N²) definitions.

**The encoder's psychoacoustics** are textbook material, not from any
encoder: the Bark scale and Terhardt's threshold in quiet, masker-dependent
offsets of 14.5 + z dB (tonal) and 5.5 dB (noise), a two-slope spreading
function. Its codebooks are designed at start-up by a Huffman construction
over Laplacian probabilities; the scale of each was chosen by measuring this
encoder's own bit rate.

**The test vectors**: the Xiph.Org Vorbis test vectors at
<https://people.xiph.org/~xiphmont/test-vectors/vorbis/>, fetched
2026-10-02, SHA-256 in `tests/vectors/SHA256SUMS`. No reference decoded
output is published with them (none was found as data), so none is compared.

## Places the specification needed reading

- **Floor 0, the running `last` (6.2.2).** Step 11 says "continue at step
  6", and step 6 sets `last` to zero — which would make steps 8 and 9 (add
  `last` to each scalar, then keep the vector's last scalar) do nothing.
  The decoder carries `last` from one vector to the next (continuing at step
  7), the only reading in which those steps mean anything, and the one under
  which the LSP angles accumulate across vectors. The test vectors settle
  it: decoded this way the four floor 0 vectors, encodings of one source,
  agree with each other at 11.6–15.5 dB waveform SNR; with `last` reset per
  vector all four decode to non-finite values.
- **End of packet during floor decode.** 6.2.2 and 7.2.3 say a floor that
  runs out of bits is "unused"; 4.3.2 says an end of packet during floor
  decode zeroes every channel. The decoder does what 4.3.2 says (the
  packet-level rule), which implies the per-floor one.
- **A codebook with no used entries.** Strictly an underspecified tree, but
  `one-entry-codebook-test.ogg` carries one in a book its audio never reads.
  The decoder accepts such a book and reports end-of-packet if it is ever
  read from; every other incomplete tree is refused.
- **An audio packet naming a mode the setup lacks.** 4.3.1 reads the mode
  number but says nothing of one past the end. `unused-mode-test.ogg` has 34
  such packets; the decoder skips them (an error in strict mode), so that
  stream decodes 24 959 samples shorter than its granule positions say.
- **Window flags.** The window shape follows the long block's two flag
  bits, as 4.3.1 specifies, not the neighbours' actual sizes; for streams
  whose flags agree with their neighbours (as 1.3.2 has a conforming encoder write them) the two
  are the same.
- **Positive and negative starting granule positions** (appendix A.2): a
  first audio page whose granule position is less than the samples its
  packets return trims the difference from the front; one greater places the
  stream later and trims nothing. On the end-of-stream page the granule
  position trims the end. A stream whose only audio page is also its last is
  trimmed at the end only.
- **Inverse MDCT scale.** The specification does not state one. The decoder
  uses the unnormalised sum of 4.3.7's definition; the test vectors decode
  at plausible levels (peaks 0.18–1.49 of full scale, RMS −35 to −12 dBFS)
  under it.
