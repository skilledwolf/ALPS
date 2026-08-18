#!/usr/bin/env bash
set -euo pipefail

: "${PLATFORM:?PLATFORM is required}"
: "${C_COMPILER:?C_COMPILER is required}"
: "${CXX_COMPILER:?CXX_COMPILER is required}"
: "${CXX_STANDARD:?CXX_STANDARD is required}"
: "${BOOST_VERSION:?BOOST_VERSION is required}"
: "${PERMISSIVE:?PERMISSIVE is required}"

workspace="${GITHUB_WORKSPACE:-$PWD}"
cc="$C_COMPILER"
cxx="$CXX_COMPILER"
if [[ "$PLATFORM" == "linux" ]]; then
  : "${C_VERSION:?C_VERSION is required for Linux}"
  cc="${C_COMPILER}-${C_VERSION}"
  cxx="${CXX_COMPILER}-${C_VERSION}"
fi

permissive_option=OFF
if [[ "$PERMISSIVE" == "true" ]]; then
  permissive_option=ON
fi

sanitizer_flags=""
if [[ -n "${SANITIZE:-}" ]]; then
  sanitizer_flags="-fsanitize=${SANITIZE} -fno-sanitize-recover=all -fno-omit-frame-pointer"
fi

cmake_args=(
  -S "$workspace"
  -B build
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_C_COMPILER="$cc"
  -DCMAKE_CXX_COMPILER="$cxx"
  -DCMAKE_C_COMPILER_LAUNCHER=ccache
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
  -DCMAKE_CXX_STANDARD="$CXX_STANDARD"
  -DALPS_ENABLE_LEGACY_PERMISSIVE="$permissive_option"
  -DBoost_SRC_DIR="$workspace/boost_1_${BOOST_VERSION}_0"
  -DPython_ROOT_DIR="$(python -c 'import sys; print(sys.prefix)')"
)

if [[ -n "$sanitizer_flags" ]]; then
  cmake_args+=(
    -DCMAKE_CXX_FLAGS="$sanitizer_flags"
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=${SANITIZE}"
    -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=${SANITIZE}"
    -DCMAKE_MODULE_LINKER_FLAGS="-fsanitize=${SANITIZE}"
  )
fi

cmake "${cmake_args[@]}"

if [[ "$PLATFORM" == "macos" ]]; then
  parallelism="$(sysctl -n hw.ncpu)"
else
  parallelism="$(nproc)"
fi
cmake --build build --parallel "$parallelism"
ccache --show-stats
