#!/usr/bin/env bash
set -euo pipefail

: "${COMP_PACK:?COMP_PACK is required}"
: "${C_COMPILER:?C_COMPILER is required}"
: "${C_VERSION:?C_VERSION is required}"
: "${BOOST_VERSION:?BOOST_VERSION is required}"
: "${BOOST_SHA256:?BOOST_SHA256 is required}"

printf 'Acquire::Retries "2";\nAcquire::http::Timeout "15";\nAcquire::https::Timeout "15";\n' \
  | sudo tee /etc/apt/apt.conf.d/99ci-timeouts > /dev/null

# Hosted runners prefer a regional Azure mirror that can stall under a large
# matrix. Use the canonical Ubuntu mirror list and put a hard wall-clock bound
# around every index refresh.
if [[ -f /etc/apt/apt-mirrors.txt ]]; then
  if [[ "$(dpkg --print-architecture)" == "arm64" ]]; then
    ubuntu_mirror=https://ports.ubuntu.com/ubuntu-ports
  else
    ubuntu_mirror=https://archive.ubuntu.com/ubuntu
  fi
  printf '%s\n' "$ubuntu_mirror" | sudo tee /etc/apt/apt-mirrors.txt > /dev/null
fi

sudo timeout 300 apt-get update
sudo apt-get install -y --no-install-recommends \
  ca-certificates ccache curl libblas-dev libhdf5-serial-dev liblapack-dev libopenmpi-dev

python -m pip install --upgrade pip
python -m pip install "numpy>=1.26" "scipy>=1.13"

if [[ "$C_COMPILER" == "clang" ]] && (( C_VERSION >= 19 )); then
  key_file="${RUNNER_TEMP:-/tmp}/apt.llvm.org.asc"
  curl --fail --location --retry 3 --retry-delay 5 \
    https://apt.llvm.org/llvm-snapshot.gpg.key --output "$key_file"
  printf '%s  %s\n' \
    '8b2a587ffd672c4687e7581dad4b2f6c1bb2ad6b480cd9771ba2ff48e0b8c75d' \
    "$key_file" | sha256sum --check --strict -
  sudo install -m 0644 "$key_file" /usr/share/keyrings/apt.llvm.org.asc

  # /etc/os-release is standard on every supported Ubuntu runner.
  # shellcheck disable=SC1091
  source /etc/os-release
  repo_suite="llvm-toolchain-${VERSION_CODENAME}-${C_VERSION}"
  printf 'deb [signed-by=/usr/share/keyrings/apt.llvm.org.asc] https://apt.llvm.org/%s/ %s main\n' \
    "$VERSION_CODENAME" "$repo_suite" \
    | sudo tee /etc/apt/sources.list.d/apt.llvm.org.list > /dev/null
  sudo timeout 300 apt-get update
elif [[ "$C_COMPILER" == "gcc" ]] && (( C_VERSION >= 15 )); then
  sudo apt-get install -y --no-install-recommends software-properties-common
  sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
  sudo timeout 300 apt-get update
fi

read -r -a compiler_packages <<< "$COMP_PACK"
sudo apt-get install -y --no-install-recommends "${compiler_packages[@]}"

archive="boost_1_${BOOST_VERSION}_0.tar.gz"
if [[ ! -f "$archive" ]]; then
  curl --fail --location --retry 3 --retry-delay 5 \
    "https://archives.boost.io/release/1.${BOOST_VERSION}.0/source/${archive}" \
    --output "$archive"
fi
printf '%s  %s\n' "$BOOST_SHA256" "$archive" | sha256sum --check --strict -
tar -xzf "$archive"
