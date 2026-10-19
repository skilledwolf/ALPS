// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Terms on three and more sites (as <NSITETERM>s compile, ALPSim/ALPS PR 223)
// against ALPS::ed: spin and fermion products spread over the chain, so that
// their factors fall into one, two, three or four of the engine's blocks,
// their measurement, the warm-up check's product energies, and the
// Jordan-Wigner string and Trotter restrictions.
#include "support.hpp"
#include <gtest/gtest.h>

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

using namespace classic_test;

namespace {

mdl::TermSum on(mdl::SiteBasisPtr const& site, char const* name, int i) { return mdl::TermSum::located(site, name, i); }

// An open S=1/2 Heisenberg chain with a ring exchange on sites 0, 2, 5, 7, a
// three-spin term on 1, 4, 6 and a four-spin term on 0, 3, 4, 7.
mdl::Hamiltonian spin_chain() {
    auto model = mdl::heisenberg(lat::chain(8, 1.0, lat::Boundary::open), 1.0);
    auto const s = mdl::spin(0.5);
    model.add(on(s, "Sp", 0) * on(s, "Sm", 2) * on(s, "Sp", 5) * on(s, "Sm", 7) * std::complex<double>(0.6), true);
    model.add(on(s, "Sz", 1) * on(s, "Sz", 4) * on(s, "Sz", 6) * std::complex<double>(0.9));
    model.add(on(s, "Sz", 0) * on(s, "Sz", 3) * on(s, "Sz", 4) * on(s, "Sz", 7) * std::complex<double>(-1.3));
    return model;
}

// An open Hubbard chain with a four-site pair hop on 0, 2, 3, 5, a hop from
// 1 to 5 weighted by the down density on 4 (the string passes an even factor),
// and a pair hop between sites 0 and 3.
mdl::Hamiltonian fermion_chain() {
    auto model = mdl::hubbard(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 3.0);
    auto const c = mdl::spinful_fermion();
    model.add(on(c, "c_dag_up", 0) * on(c, "c_dag_dn", 2) * on(c, "c_dn", 5) * on(c, "c_up", 3)
                  * std::complex<double>(0.7), true);
    model.add(on(c, "n_dn", 4) * on(c, "c_dag_up", 1) * on(c, "c_up", 5) * std::complex<double>(-0.45), true);
    model.add(on(c, "c_dag_up", 3) * on(c, "c_dag_dn", 3) * on(c, "c_dn", 0) * on(c, "c_up", 0)
                  * std::complex<double>(0.3), true);
    return model;
}

void expect_values(cd::Measurement const& measured, std::vector<double> const& reference) {
    ASSERT_EQ(measured.values.size(), 1u);
    ASSERT_EQ(measured.values[0].size(), reference.size());
    for (std::size_t k = 0; k < reference.size(); ++k) EXPECT_NEAR(measured.values[0][k], reference[k], 1e-8) << k;
}

} // namespace

TEST(ClassicDmrgNSite, SpinTermsOnThreeAndFourSitesMatchExactDiagonalization) {
    auto const model = spin_chain();
    Sector const sector{{"Sz", 0.}};
    auto const result = cd::solve(model, exact_options(64, sector, 2));
    auto const reference = ed_energies(model, sector, 2);
    ASSERT_EQ(result.energies.size(), 2u);
    for (std::size_t k = 0; k < 2; ++k) EXPECT_NEAR(result.energies[k], reference[k], 1e-8) << k;
}

TEST(ClassicDmrgNSite, FermionProductsOnThreeAndFourSitesMatchExactDiagonalization) {
    auto const model = fermion_chain();
    for (Sector const& sector : {Sector{{"N", 4.}, {"Sz", 0.}}, Sector{{"N", 5.}, {"Sz", 0.5}}}) {
        auto const result = cd::solve(model, exact_options(128, sector));
        EXPECT_NEAR(result.energies.at(0), ed_energies(model, sector).at(0), 1e-8);
    }
}

TEST(ClassicDmrgNSite, NSiteObservablesMatchExactDiagonalization) {
    auto const model = fermion_chain();
    auto const c = mdl::spinful_fermion();
    Sector const sector{{"N", 4.}, {"Sz", 0.}};
    mdl::TermSum pair_hop = on(c, "c_dag_up", 0) * on(c, "c_dag_dn", 2) * on(c, "c_dn", 5) * on(c, "c_up", 3);
    pair_hop += pair_hop.dag();
    auto const ring = model.operator_decomposition(pair_hop).op;
    auto const density_hop = model.operator_decomposition(on(c, "n_dn", 4) * on(c, "c_dag_up", 1) * on(c, "c_up", 5)
                                                          + on(c, "n_dn", 4) * on(c, "c_dag_up", 5) * on(c, "c_up", 1)).op;
    cd::Measurements requests;
    requests.expectation["pair hop"] = ring;
    requests.terms["density hop"] = density_hop;
    auto const result = cd::solve(model, exact_options(128, sector), requests);
    auto const states = ed_states(model, sector);
    expect_values(result.measurements.at("pair hop"), {ed_value(states, ring)});
    expect_values(result.measurements.at("density hop"), ed_term_values(states, density_hop));
    EXPECT_EQ(result.measurements.at("density hop").sites.at(0), (std::vector<int>{1, 4, 5}));
}

TEST(ClassicDmrgNSite, WarmupCheckCountsTermsOnManySites) {
    // Diagonal only: staggered fields fix the Neel state, and the four-spin
    // term shifts its energy by -1.3/16. The product state the check finds
    // must carry that shift to equal the exact ground state.
    auto model = mdl::heisenberg(lat::chain(8, 1.0, lat::Boundary::open), 0.0);
    auto const s = mdl::spin(0.5);
    for (int i = 0; i < 8; ++i) model.add(on(s, "Sz", i) * std::complex<double>(i % 2 ? 1.5 : -1.5));
    model.add(on(s, "Sz", 0) * on(s, "Sz", 3) * on(s, "Sz", 4) * on(s, "Sz", 7) * std::complex<double>(-1.3));
    Sector const sector{{"Sz", 0.}};
    auto const result = cd::solve(model, exact_options(16, sector));
    auto const exact = ed_energies(model, sector).at(0);
    EXPECT_NEAR(result.energies.at(0), exact, 1e-10);
    ASSERT_TRUE(result.warmup_check);
    EXPECT_NEAR(result.warmup_check->product_energy, exact, 1e-12);
    EXPECT_FALSE(result.warmup_check->misplaced);
}

TEST(ClassicDmrgNSite, StringsOtherThanTheStrictConventionsAreRejected) {
    auto const model = fermion_chain();
    auto op = mdl::compile_hamiltonian(model).op;
    auto const basis = mdl::physical_basis(model);
    bool changed = false;
    for (auto& term : op.terms)
        if (term.factors.size() == 4 && !term.jordan_wigner.empty()) {
            term.jordan_wigner.pop_back();
            changed = true;
            break;
        }
    ASSERT_TRUE(changed);
    EXPECT_THROW(cd::solve(op, basis, exact_options(16, {{"N", 4.}, {"Sz", 0.}})), std::invalid_argument);
}

TEST(ClassicDmrgNSite, TrotterGatesRefuseTermsOnMoreThanTwoSites) {
    auto const model = spin_chain();
    cd::Evolution evolution;
    evolution.time_step = 0.05;
    evolution.steps = 2;
    evolution.states = 16;
    evolution.method = cd::EvolutionMethod::trotter;
    auto const op = mdl::compile_hamiltonian(model).op;
    auto const basis = mdl::physical_basis(model);
    try {
        cd::evolve(op, op, basis, evolution, exact_options(16, {{"Sz", 0.}}));
        FAIL() << "a Trotter evolution of a four-site term was accepted";
    } catch (std::invalid_argument const& error) {
        EXPECT_NE(std::string(error.what()).find("sites"), std::string::npos) << error.what();
    }
}
