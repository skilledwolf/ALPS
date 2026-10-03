// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Numeric lock for FourierTransformer::backward_ft's off-diagonal
// "nothing happening in this gf" short-circuit (fouriertransform.C).
//
// The guard zeroes a (flavor, s1, s2) channel whose entire high-frequency
// tail is absent, i.e. c1 == c2 == c3 == 0. The diagonal channels always
// have c1 == 1 (set by the base ctor), so this branch is exercised only on
// the OFF-diagonal (s1 != s2) multi-site path — which no production
// transformer reaches today (every concrete transformer forces n_site == 1
// and the generate_* factories throw on SITES != 1). This test constructs a
// bare two-site FourierTransformer directly to pin the intended behaviour
// ahead of any future cluster (n_site > 1) support.
//
// Regression intent: the guard's third clause used to be a bare truthiness
// test on c3 (firing when c3 != 0) instead of `c3 == 0`, the opposite of its
// two siblings and of the comment. With that typo an all-zero-tail
// off-diagonal channel took the full-transform branch and leaked the bogus
// off-diagonal input into G_tau instead of being zeroed. This test fails on
// the pre-fix code and passes once the guard reads c1==0 && c2==0 && c3==0.

#include <boost/math/constants/constants.hpp>
#include "fouriertransform.h"
#include "U_matrix.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>

namespace {
bool typed_input_contract() {
  alps::params parameters;
  parameters["BETA"] = 10.;
  parameters["MU"] = 0.;
  parameters["FLAVORS"] = 2;
  boost::shared_ptr<FourierTransformer> transformer;
  auto invalid = parameters;
  invalid["FLAVORS"] = 2.5;
  try {
    FourierTransformer::generate_transformer(invalid, transformer);
    std::puts("FAIL: fractional flavor count was truncated");
    return false;
  } catch (const alps::params_ns::exception::type_mismatch&) {}
  invalid = parameters;
  invalid["EPS_0"] = "0.25";
  try {
    FourierTransformer::generate_transformer(invalid, transformer);
    std::puts("FAIL: string band moment was coerced to a scientific number");
    return false;
  } catch (const alps::params_ns::exception::type_mismatch&) {}
  parameters["FLAVORS"] = 4;
  parameters["U"] = 4.;
  parameters["J"] = 1.;
  U_matrix interaction(parameters);
  if (interaction(0, 1) != 4. || interaction(0, 2) != 1. || interaction(0, 3) != 2. ||
      interaction(1, 3) != 1. || interaction.mu_shift() != 3.5) {
    std::puts("FAIL: native Hubbard/Hund interaction convention changed");
    return false;
  }
  return true;
}
}

int main() {
  if (!typed_input_contract()) return 1;
  const double   beta     = 10.0;
  const unsigned N_omega  = 200;   // Matsubara frequencies
  const unsigned N_tau    = 400;   // imaginary-time slices
  const unsigned n_site   = 2;
  const unsigned n_flavor = 1;

  // Bare transformer: diagonal tail c1 == 1, off-diagonal tail all zero
  // (exactly what the base constructor installs).
  FourierTransformer ft(beta, static_cast<int>(n_flavor), static_cast<int>(n_site));

  matsubara_green_function_t G_omega(N_omega, n_site, n_flavor);
  for (unsigned k = 0; k < N_omega; ++k) {
    const std::complex<double> iw(0.0, (2 * k + 1) * boost::math::constants::pi<double>() / beta);
    G_omega(k, 0, 0, 0) = 1.0 / (iw - 0.3);                  // real diagonal data
    G_omega(k, 1, 1, 0) = 1.0 / (iw + 0.3);
    G_omega(k, 0, 1, 0) = std::complex<double>(0.05, -0.02); // bogus off-diagonal
    G_omega(k, 1, 0, 0) = std::complex<double>(0.05, -0.02);
  }

  itime_green_function_t G_tau(N_tau + 1, n_site, n_flavor);
  ft.backward_ft(G_tau, G_omega);

  double max_offdiag = 0.0;
  double max_diag    = 0.0;
  for (unsigned i = 0; i <= N_tau; ++i) {
    max_offdiag = std::max(max_offdiag, std::abs(G_tau(i, 0, 1, 0)));
    max_offdiag = std::max(max_offdiag, std::abs(G_tau(i, 1, 0, 0)));
    max_diag    = std::max(max_diag, std::abs(G_tau(i, 0, 0, 0)));
  }

  bool ok = true;
  if (max_offdiag > 1e-12) {
    std::printf("FAIL: off-diagonal G_tau not zeroed (max=%g); the "
                "c1==c2==c3==0 'nothing happening' guard did not fire.\n",
                max_offdiag);
    ok = false;
  }
  if (max_diag < 1e-6) {
    std::printf("FAIL: diagonal G_tau is trivially zero (max=%g); the "
                "transform did not run.\n",
                max_diag);
    ok = false;
  }
  if (ok) {
    std::printf("OK: off-diagonal channel zeroed (max=%g), diagonal "
                "non-trivial (max=%g)\n",
                max_offdiag, max_diag);
  }
  return ok ? 0 : 1;
}
