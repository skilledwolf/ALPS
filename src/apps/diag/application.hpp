// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "lattice_model.hpp"
#include <complex>
#include <ostream>

namespace diag {
// Bloch states need complex arithmetic; COMPLEX overrides the choice.
inline bool complex_matrix(alps::Parameters const& p) {
    return p.value_or_default("COMPLEX", bool(p.value_or_default("TRANSLATION_SYMMETRY", true))
                                         || p.defined("TOTAL_MOMENTUM"));
}

// Diagonalizes each run with Matrix<double> or Matrix<std::complex<double>>.
template <template <class> class Matrix>
struct application {
    static void print_copyright(std::ostream& out) { Matrix<double>::print_copyright(out); }

    template <class F>
    static void with_task(alps::run_configuration const& run, F const& f) {
        const auto p = lattice_model::parameters(run);
        if (complex_matrix(p)) { Matrix<std::complex<double>> matrix(p); f(matrix); }
        else { Matrix<double> matrix(p); f(matrix); }
    }
};
} // namespace diag
