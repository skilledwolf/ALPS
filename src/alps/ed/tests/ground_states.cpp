// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ground states by prepare_product + eigenstates: the charge sectors
// partition the full space and its ground state, Lanczos agrees with the
// dense solver where the space is large enough, and the exact
// Majumdar-Ghosh and recorded L=16 Heisenberg values are reproduced. The
// spectra of the model factories themselves are checked against the XML
// catalog in src/alps/model/tests/catalog_model_parity.cpp.

#include <alps/ed/product_problem.hpp>
#include <alps/models/builtin.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;

namespace {

// Ground state of the L=16 periodic Heisenberg chain in Sz=0 (dimension
// 12870), measured once with alps::ed::dense_eigensystem on the same
// sector (131 s; the Lanczos path below agrees to 2e-15 in 0.1 s). The
// Bethe-ansatz finite-size value is -7.1422963606.
constexpr double kHeisenberg16E0 = -7.142296360616779;

alps::ed::EigenstateOptions options(int count, alps::ed::EigenstateMethod method) {
    alps::ed::EigenstateOptions o;
    o.count = count;
    o.method = method;
    return o;
}

alps::ed::ProductOptions sector(std::map<std::string, double> charges) {
    alps::ed::ProductOptions o;
    o.sector = std::move(charges);
    return o;
}

} // namespace

TEST(GroundStates, SectorsPartitionTheSpectrum) {
    // The lowest sector ground state over every charge sector is the full
    // space's; empty sectors return no eigenvalues. The even
    // antiferromagnetic chains hold their singlet ground state in Sz=0.
    struct Case {
        mdl::Hamiltonian model;
        std::vector<std::map<std::string, double>> sectors;
        std::map<std::string, double> ground;
    };
    std::vector<Case> cases;
    auto spin_sectors = [](int L) {
        std::vector<std::map<std::string, double>> out;
        for (int twice = -L; twice <= L; twice += 2) out.push_back({{"Sz", twice / 2.}});
        return out;
    };
    cases.push_back({mdl::heisenberg(lat::chain(6), 1.0), spin_sectors(6), {{"Sz", 0.}}});
    cases.push_back({mdl::xxz(lat::chain(6), 1.0, 0.5), spin_sectors(6), {{"Sz", 0.}}});
    cases.push_back({mdl::j1j2(lat::nnn_chain(6), 1.0, 0.5), spin_sectors(6), {{"Sz", 0.}}});
    {
        Case hubbard{mdl::hubbard(lat::chain(4), 1.0, 2.0, 0.3), {}, {}};
        for (int n = 0; n <= 8; ++n)
            for (int twice = -4; twice <= 4; ++twice) hubbard.sectors.push_back({{"N", double(n)}, {"Sz", twice / 2.}});
        cases.push_back(hubbard);
    }
    {
        Case spinless{mdl::spinless_fermion_model(lat::chain(6), 1.0, 0.2, 0.5), {}, {}};
        for (int n = 0; n <= 6; ++n) spinless.sectors.push_back({{"N", double(n)}});
        cases.push_back(spinless);
    }
    for (std::size_t k = 0; k < cases.size(); ++k) {
        SCOPED_TRACE("case " + std::to_string(k));
        auto const& c = cases[k];
        auto const full = alps::ed::prepare_product(c.model);
        double const e0 = full.eigenstates(options(1, alps::ed::EigenstateMethod::dense)).eigenvalues.at(0);
        double lowest = std::numeric_limits<double>::infinity();
        std::size_t states = 0;
        for (auto const& charges : c.sectors) {
            auto const problem = alps::ed::prepare_product(c.model, sector(charges));
            states += problem.dimension();
            auto const result = problem.eigenstates(options(1, alps::ed::EigenstateMethod::dense));
            if (problem.dimension() == 0) { EXPECT_TRUE(result.eigenvalues.empty()); continue; }
            ASSERT_EQ(result.eigenvalues.size(), 1);
            lowest = std::min(lowest, result.eigenvalues[0]);
        }
        EXPECT_EQ(states, full.dimension());
        EXPECT_NEAR(lowest, e0, 1e-10);
        if (!c.ground.empty()) {
            auto const result = alps::ed::prepare_product(c.model, sector(c.ground))
                .eigenstates(options(1, alps::ed::EigenstateMethod::dense));
            EXPECT_NEAR(result.eigenvalues.at(0), e0, 1e-10);
        }
    }
}

TEST(GroundStates, LanczosAgreesWithDenseInTheSingletSector) {
    auto const problem = alps::ed::prepare_product(mdl::heisenberg(lat::chain(12)), sector({{"Sz", 0.}}));
    ASSERT_EQ(problem.dimension(), 924);
    auto const dense = problem.eigenstates(options(4, alps::ed::EigenstateMethod::dense));
    auto const lanczos = problem.eigenstates(options(4, alps::ed::EigenstateMethod::lanczos));
    ASSERT_EQ(lanczos.eigenvalues.size(), 4);
    EXPECT_EQ(lanczos.method, alps::ed::EigenstateMethod::lanczos);
    for (std::size_t i = 0; i < 4; ++i) EXPECT_NEAR(lanczos.eigenvalues[i], dense.eigenvalues[i], 1e-8) << i;
    for (double residual : lanczos.residual_norms) EXPECT_LT(residual, 1e-7);
}

TEST(GroundStates, MajumdarGhoshPointIsExact) {
    // J2 = J1/2: the dimer product states are exact ground states, E0 = -3L/8.
    auto const problem = alps::ed::prepare_product(mdl::j1j2(lat::nnn_chain(6), 1.0, 0.5));
    EXPECT_NEAR(problem.eigenstates(options(1, alps::ed::EigenstateMethod::dense)).eigenvalues.at(0),
                -3.0 / 8.0 * 6, 1e-10);
}

TEST(GroundStates, HeisenbergChain16GroundStateByLanczos) {
    // L=16 periodic Heisenberg chain, Sz=0 sector (dimension 12870): the
    // automatic method selects Lanczos. The reference is the dense ground
    // state of the same sector, computed once (see the PR text); it agrees
    // with the Bethe-ansatz finite-size value.
    auto const problem = alps::ed::prepare_product(mdl::heisenberg(lat::chain(16)), sector({{"Sz", 0.}}));
    ASSERT_EQ(problem.dimension(), 12870);
    auto const result = problem.eigenstates(options(1, alps::ed::EigenstateMethod::automatic));
    EXPECT_EQ(result.method, alps::ed::EigenstateMethod::lanczos);
    ASSERT_EQ(result.eigenvalues.size(), 1);
    EXPECT_NEAR(result.eigenvalues[0], kHeisenberg16E0, 1e-8);
}
