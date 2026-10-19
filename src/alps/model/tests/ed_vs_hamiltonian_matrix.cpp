// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The catalog importer plus ALPS::ed against the legacy hamiltonian_matrix on
// the same catalog inputs and the same legacy basis enumeration: the matrices
// must agree element by element (fermionic signs included) and so must the
// complete spectra. The inputs follow the diag applications' sample
// parameter files (src/apps/test/parm*.in), sector by sector, plus n-site
// terms (<NSITETERM>, ALPSim/ALPS PR 223) on lattices and graphs whose groups
// list their sites out of order.
#include <alps/ed/eigensolver.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/model/catalog.hpp>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/numeric/matrix.hpp>
#include <alps/utility/vectorio.hpp>
#include <gtest/gtest.h>
#include "nsite_libraries.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace {
using Params = std::vector<std::pair<std::string, std::string>>;
struct Case {
    std::string name;
    Params params;
};
std::vector<Case> const& cases() {
    static std::vector<Case> const r = {
        {"spin chain, full space",
         {{"MODEL", "spin"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"J", "1"}, {"Jxy", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "Sz"}}},
        {"spin chain, Sz=0",
         {{"MODEL", "spin"}, {"LATTICE", "chain lattice"}, {"L", "6"}, {"J", "1"}, {"Jxy", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Sz_total", "0"}}},
        {"spin-1 ladder with anisotropy and field, Sz=1",
         {{"MODEL", "spin"}, {"LATTICE", "ladder"}, {"L", "3"}, {"local_S", "1"}, {"J0", "1"}, {"J1", "0.4"},
          {"D", "0.2"}, {"h", "0.1"}, {"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Sz_total", "1"}}},
        {"hardcore bosons, N=2",
         {{"MODEL", "hardcore boson"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"V", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "2"}}},
        {"boson Hubbard, N=3",
         {{"MODEL", "boson Hubbard"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"Nmax", "2"}, {"t", "1"},
          {"U", "4"}, {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "3"}}},
        {"spinless fermions, periodic, N=2",
         {{"MODEL", "spinless fermions"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"V", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "2"}}},
        {"spinless fermions, open, N=3",
         {{"MODEL", "spinless fermions"}, {"LATTICE", "open chain lattice"}, {"L", "5"}, {"t", "1"},
          {"V", "0.5"}, {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "3"}}},
        {"fermion Hubbard chain, Nup=2 Ndown=1",
         {{"MODEL", "fermion Hubbard"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"U", "4"},
          {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "2"}, {"Ndown_total", "1"}}},
        {"fermion Hubbard square, Nup=2 Ndown=2",
         {{"MODEL", "fermion Hubbard"}, {"LATTICE", "square lattice"}, {"L", "2"}, {"W", "2"}, {"t", "1"},
          {"U", "4"}, {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "2"}, {"Ndown_total", "2"}}},
        {"Ian's fermion Hubbard, N=2 Sz=0",
         {{"MODEL", "Ian's fermion Hubbard"}, {"LATTICE", "chain lattice"}, {"L", "6"}, {"t", "1"}, {"U", "2"},
          {"CONSERVED_QUANTUMNUMBERS", "N,Sz"}, {"N_total", "2"}, {"Sz_total", "0"}}},
        {"t-J chain, N=3 Sz=1/2",
         {{"MODEL", "t-J"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"J", "0.4"},
          {"CONSERVED_QUANTUMNUMBERS", "N,Sz"}, {"N_total", "3"}, {"Sz_total", "0.5"}}},
    };
    return r;
}

std::vector<double> sorted_spectrum(std::vector<double> matrix, std::size_t dimension) {
    return alps::ed::dense_eigensystem(std::move(matrix), dimension, false).eigenvalues;
}

// The catalog importer's operator through ALPS::ed against the legacy
// hamiltonian_matrix on the legacy basis, element by element, then the spectra.
void expect_legacy_matrix(alps::Parameters const& p) {
    alps::hamiltonian_matrix<alps::numeric::matrix<double>> legacy(p);
    auto const& reference = legacy.matrix();
    auto const dimension = legacy.dimension();
    ASSERT_GT(dimension, 1u);
    auto const& basis = legacy.states_vector();
    ASSERT_EQ(basis.size(), dimension);
    std::size_t const sites = basis[0].size();
    std::vector<std::int64_t> states;
    for (auto const& state : basis)
        for (auto digit : state) states.push_back(digit);

    alps::models::CatalogModel<> model(p);
    ASSERT_EQ(model.basis().size(), sites);
    auto const op = model.compile_hamiltonian();
    auto const matrix = alps::ed::build_csr_from_local<std::complex<double>>(
        op, states, dimension, model.local_parity());
    std::vector<double> dense(dimension * dimension), expected(dimension * dimension);
    for (std::size_t row = 0; row < dimension; ++row)
        for (auto i = matrix.offsets[row]; i < matrix.offsets[row + 1]; ++i) {
            auto const value = matrix.values[static_cast<std::size_t>(i)];
            EXPECT_EQ(value.imag(), 0.);
            dense[static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(i)]) * dimension + row] =
                value.real();
        }
    double scale = 1.;
    for (std::size_t row = 0; row < dimension; ++row)
        for (std::size_t column = 0; column < dimension; ++column) {
            expected[column * dimension + row] = reference(row, column);
            scale = std::max(scale, std::abs(reference(row, column)));
        }
    double difference = 0., asymmetry = 0.;
    for (std::size_t row = 0; row < dimension; ++row)
        for (std::size_t column = 0; column < dimension; ++column) {
            difference = std::max(difference,
                std::abs(dense[column * dimension + row] - expected[column * dimension + row]));
            asymmetry = std::max(asymmetry,
                std::abs(expected[column * dimension + row] - expected[row * dimension + column]));
        }
    EXPECT_LE(asymmetry, 1e-12 * scale);
    EXPECT_LE(difference, 1e-12 * scale);

    auto const spectrum = sorted_spectrum(dense, dimension), reference_spectrum = sorted_spectrum(expected, dimension);
    ASSERT_EQ(spectrum.size(), reference_spectrum.size());
    for (std::size_t i = 0; i < spectrum.size(); ++i)
        EXPECT_NEAR(spectrum[i], reference_spectrum[i], 1e-10 * scale);
}

} // namespace

TEST(EdVsHamiltonianMatrix, CompiledCatalogOperatorsReproduceTheLegacyMatrices) {
    for (auto const& item : cases()) {
        SCOPED_TRACE(item.name);
        alps::Parameters p;
        for (auto const& [key, value] : item.params) p[key] = value;
        expect_legacy_matrix(p);
    }
}

// n-site terms: four-, three- and two-site fermion products whose roles are
// out of site order (and wrap around a periodic ladder), two fermion
// operators on one role, a constant, spins and truncated bosons.
TEST(EdVsHamiltonianMatrix, CompiledNSiteTermsReproduceTheLegacyMatrices) {
    nsite_test::Libraries const libraries;
    Params const files{{"LATTICE_LIBRARY", (libraries.dir / "lattices.xml").string()},
                       {"MODEL_LIBRARY", (libraries.dir / "models.xml").string()}};
    std::vector<Case> const nsite_cases = {
        {"pair hopping, open ladder, Nup=2 Ndown=2",
         {{"MODEL", "pair hopping"}, {"LATTICE", "open ladder"}, {"L", "3"},
          {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "2"}, {"Ndown_total", "2"}}},
        {"pair hopping, periodic ladder, Nup=3 Ndown=2",
         {{"MODEL", "pair hopping"}, {"LATTICE", "periodic ladder"}, {"L", "3"},
          {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "3"}, {"Ndown_total", "2"}}},
        {"three-site fermions, scrambled graph, Nup=2 Ndown=2",
         {{"MODEL", "three-site fermions"}, {"GRAPH", "scrambled"},
          {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "2"}, {"Ndown_total", "2"}}},
        {"three-site fermions, scrambled graph, full space",
         {{"MODEL", "three-site fermions"}, {"GRAPH", "scrambled"}}},
        {"ring exchange, periodic ladder, Sz=0",
         {{"MODEL", "ring exchange"}, {"LATTICE", "periodic ladder"}, {"L", "4"},
          {"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Sz_total", "0"}}},
        {"boson ring hop, periodic ladder, N=3",
         {{"MODEL", "boson ring hop"}, {"LATTICE", "periodic ladder"}, {"L", "2"},
          {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "3"}}},
    };
    for (auto const& item : nsite_cases) {
        SCOPED_TRACE(item.name);
        alps::Parameters p;
        for (auto const& [key, value] : files) p[key] = value;
        for (auto const& [key, value] : item.params) p[key] = value;
        expect_legacy_matrix(p);
    }
}

// Momentum sectors of a periodic ladder with n-site terms: ALPS::ed on the
// catalog's translations against the legacy Bloch matrices, momentum by
// momentum, as sparsediag loops TOTAL_MOMENTUM over translation_momenta().
TEST(EdVsHamiltonianMatrix, NSiteTermsInEveryMomentumSector) {
    nsite_test::Libraries const libraries;
    alps::Parameters p;
    p["LATTICE_LIBRARY"] = (libraries.dir / "lattices.xml").string();
    p["MODEL_LIBRARY"] = (libraries.dir / "models.xml").string();
    p["MODEL"] = "pair hopping";
    p["LATTICE"] = "periodic ladder";
    p["L"] = 4;
    p["CONSERVED_QUANTUMNUMBERS"] = "Nup,Ndown";
    p["Nup_total"] = 2;
    p["Ndown_total"] = 2;
    alps::models::CatalogModel<> model(p);
    std::vector<alps::ed::SpectrumTranslation> translations;
    for (auto const& t : model.translations()) translations.push_back({t.sites, t.period});
    ASSERT_EQ(translations.size(), 1u);
    auto const op = model.compile_hamiltonian();
    auto const spaces = alps::ed::projected_spaces(model.physical_basis());
    alps::graph_helper<> graph(p);
    std::size_t total = 0;
    for (auto const& k : graph.translation_momenta()) {
        SCOPED_TRACE("momentum " + alps::write_vector(k));
        alps::ed::ModelSpectrumOptions options;
        options.sector = model.sector();
        options.momentum.assign(k.begin(), k.end());
        auto const spectrum = alps::ed::model_spectrum(op, spaces, translations, options);
        std::vector<double> levels;
        for (auto const& block : spectrum.blocks) levels.insert(levels.end(), block.eigenvalues.begin(), block.eigenvalues.end());
        std::sort(levels.begin(), levels.end());

        alps::Parameters q(p);
        q["TOTAL_MOMENTUM"] = alps::write_vector(k);
        alps::hamiltonian_matrix<alps::numeric::matrix<std::complex<double>>> legacy(q);
        auto const dimension = legacy.dimension();
        ASSERT_EQ(levels.size(), dimension);
        total += dimension;
        if (!dimension) continue;
        std::vector<std::complex<double>> dense(dimension * dimension);
        for (std::size_t column = 0; column < dimension; ++column)
            for (std::size_t row = 0; row < dimension; ++row) dense[column * dimension + row] = legacy.matrix()(row, column);
        auto const reference = alps::ed::dense_eigensystem(std::move(dense), dimension, false).eigenvalues;
        for (std::size_t i = 0; i < dimension; ++i) EXPECT_NEAR(levels[i], reference[i], 1e-10);
    }
    EXPECT_EQ(total, 28u * 28u);  // the momenta together span the sector: C(8,2)^2
}
