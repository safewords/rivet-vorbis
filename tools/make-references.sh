#!/usr/bin/env bash
# Regenerate the reference decodes the opt-in test in tests/reference.rs
# compares against. Xiph.Org's oggdec (vorbis-tools, from Ubuntu 24.04) runs
# in Docker strictly as a black-box binary: Ogg Vorbis in, raw 16-bit
# little-endian signed PCM out. Its source is not read and no other
# implementation is involved.
#
# Inputs: the Xiph.Org test vectors (tools/fetch-vectors.sh) and a few files
# from this crate's encoder, written by the test itself. Outputs (not
# committed): tests/references/*.raw and the encoder .ogg files. Committed:
# tests/references/SHA256SUMS, the hashes of every input and output, and
# tests/references/oggdec-version.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
bash tools/fetch-vectors.sh >/dev/null

# The encoder's outputs, deterministic for a given encoder.
VORBIS_WRITE_ENCODER_FIXTURES=1 cargo test --release --test reference write_encoder_fixtures -- --exact >/dev/null

mkdir -p tests/references
cp tests/vectors/*.ogg tests/references/ 2>/dev/null || true

mount="$root/tests/references"
if command -v cygpath >/dev/null 2>&1; then mount="$(cygpath -m "$mount")"; fi
MSYS_NO_PATHCONV=1 docker run --rm -v "$mount:/w" ubuntu:24.04 bash -c '
  set -e
  apt-get update -qq >/dev/null
  DEBIAN_FRONTEND=noninteractive apt-get install -y -qq vorbis-tools >/dev/null
  cd /w
  dpkg-query -W -f="vorbis-tools \${Version} (ubuntu:24.04)\n" vorbis-tools > oggdec-version
  for f in *.ogg; do
    oggdec -Q -R -b 16 -e 0 -s 1 -o "${f%.ogg}.raw" "$f" || echo "oggdec failed: $f" >&2
  done
'

cd tests/references
# Inputs and outputs, so a stale reference is detected.
sha256sum *.ogg *.raw | sed 's/ \*/  /' > SHA256SUMS
echo "references written to tests/references/ ($(ls *.raw | wc -l) decodes)"
