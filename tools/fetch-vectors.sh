#!/usr/bin/env bash
# Fetch the Xiph.Org Vorbis test vectors into tests/vectors/ and check them
# against tests/vectors/SHA256SUMS. They are not committed: their licence is
# not stated, and the suite only needs them on disk.
set -euo pipefail
cd "$(dirname "$0")/../tests/vectors"
base=https://people.xiph.org/~xiphmont/test-vectors/vorbis
while read -r sum name; do
  [ -f "$name" ] || curl -sSf -o "$name" "$base/$name"
done < SHA256SUMS
sha256sum -c SHA256SUMS
