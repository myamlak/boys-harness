#!/bin/bash
# Build crlibm 1.0beta5 from the pinned public mirror, then the sinh witness
# driver.  Run inside WSL Ubuntu; paths are derived from this script's
# location so the deposited directory is relocatable.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT/third_party/crlibm-1.0beta5"

echo "=== toolchain ==="
gcc --version | head -1
uname -m
echo "=== configure ==="
cd "$SRC"
if [ -f Makefile ]; then
  make distclean || true
fi
autoreconf -i
CFLAGS="-O2" ./configure --enable-sse2 --disable-mpfr --disable-double-extended
echo "=== build scs_lib ==="
make -C scs_lib
echo "=== build libcrlibm.a ==="
make libcrlibm.a
echo "=== build driver ==="
gcc -O2 -std=c99 -msse2 -mfpmath=sse -I. \
    -o "$ROOT/raw/crlibm_sinh_witness" \
    "$ROOT/src/crlibm_sinh_witness.c" libcrlibm.a -lm
echo "=== built $ROOT/raw/crlibm_sinh_witness ==="
