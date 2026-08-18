#!/usr/bin/env bash
set -euo pipefail

: "${PLATFORM:?PLATFORM is required}"
: "${C_COMPILER:?C_COMPILER is required}"

if [[ "$PLATFORM" == "linux" ]] && [[ "${SANITIZE:-}" == *address* ]]; then
  : "${C_VERSION:?C_VERSION is required for Linux sanitizers}"
  asan_runtime="$("${C_COMPILER}-${C_VERSION}" -print-file-name=libasan.so)"
  export LD_PRELOAD="${asan_runtime}${LD_PRELOAD:+:${LD_PRELOAD}}"
  export ASAN_OPTIONS=detect_leaks=0
fi

export UBSAN_OPTIONS=print_stacktrace=1
ctest --test-dir build --output-on-failure --output-junit junit.xml
