#!/usr/bin/env bash
set -euo pipefail

: "${BOOST_VERSION:?BOOST_VERSION is required}"
: "${BOOST_SHA256:?BOOST_SHA256 is required}"

python -m pip install --upgrade pip
python -m pip install "numpy>=1.26" "scipy>=1.13"
brew install ccache hdf5 openmpi

archive="boost_1_${BOOST_VERSION}_0.tar.gz"
if [[ ! -f "$archive" ]]; then
  curl --fail --location --retry 3 --retry-delay 5 \
    "https://archives.boost.io/release/1.${BOOST_VERSION}.0/source/${archive}" \
    --output "$archive"
fi
printf '%s  %s\n' "$BOOST_SHA256" "$archive" | shasum -a 256 --check
tar -xzf "$archive"
