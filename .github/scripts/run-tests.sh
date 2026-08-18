#!/usr/bin/env bash
set -euo pipefail

: "${PLATFORM:?PLATFORM is required}"
: "${C_COMPILER:?C_COMPILER is required}"
: "${CXX_COMPILER:?CXX_COMPILER is required}"

if [[ "$PLATFORM" == "linux" ]] && [[ "${SANITIZE:-}" == *address* ]]; then
  : "${C_VERSION:?C_VERSION is required for Linux sanitizers}"
  asan_runtime="$("${C_COMPILER}-${C_VERSION}" -print-file-name=libasan.so)"
  cxx_runtime="$("${CXX_COMPILER}-${C_VERSION}" -print-file-name=libstdc++.so.6)"
  # Python is a C executable. Load libstdc++ before it dlopens the instrumented
  # extension modules so ASan can install its C++ exception interceptors.
  export LD_PRELOAD="${asan_runtime}:${cxx_runtime}${LD_PRELOAD:+:${LD_PRELOAD}}"
  # matrix_unit_tests deliberately probes allocation failure with max_size().
  export ASAN_OPTIONS=allocator_may_return_null=1:detect_leaks=0
fi

export UBSAN_OPTIONS=print_stacktrace=1
ctest --test-dir build --output-on-failure --output-junit junit.xml
