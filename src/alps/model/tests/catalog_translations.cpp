// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// CatalogModel::translations(): the momentum generators the diag
// applications take from the legacy lattice. The permutations are checked
// against the kit lattices' translation_generators() on rectangles, cylinders
// and multi-vertex cells, because a reversed or swapped generator leaves every
// momentum-resolved level of a real Hamiltonian unchanged. The spectra of
// catalog Heisenberg chains and squares resolved by those generators are
// checked against the kit Hamiltonian's model_spectrum, against
// hamiltonian_matrix driven as sparsediag drives it (DiagMatrix::dostep in
// src/apps/diag/diag.h: TOTAL_MOMENTUM looped over translation_momenta()),
// and against levels recorded from sparsediag runs.
#include <alps/ed/eigensolver.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/model/catalog.hpp>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/models/builtin.hpp>
#include <alps/numeric/matrix.hpp>
#include <alps/utility/vectorio.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace lat = alps::lattice;
namespace mdl = alps::models;
using Params = std::vector<std::pair<std::string, std::string>>;
double const two_pi = 2. * std::acos(-1.);

alps::Parameters parameters(Params const& values) {
    alps::Parameters p;
    for (auto const& [key, value] : values) p[key] = value;
    return p;
}

std::vector<alps::ed::SpectrumTranslation> spectrum_translations(mdl::CatalogModel<> const& model) {
    std::vector<alps::ed::SpectrumTranslation> result;
    for (auto const& t : model.translations()) result.push_back({t.sites, t.period});
    return result;
}

// The kit lattice's generators, periods and axes, which the legacy expansion
// of its LATTICEGRAPH (and of the catalog graph it names) must reproduce.
void expect_kit_generators(std::vector<mdl::CatalogTranslation> const& translations, lat::Lattice const& lattice) {
    auto const directions = lattice.periodic_directions();
    auto const generators = lattice.translation_generators();
    ASSERT_EQ(translations.size(), directions.size());
    for (std::size_t i = 0; i < translations.size(); ++i) {
        SCOPED_TRACE("direction " + std::to_string(directions[i]));
        EXPECT_EQ(translations[i].direction, directions[i]);
        EXPECT_EQ(translations[i].period, static_cast<std::size_t>(lattice.extent()[static_cast<std::size_t>(directions[i])]));
        EXPECT_EQ(translations[i].sites, generators[i]);
    }
}

// The grid index j of k = 2 pi j / L per direction.
std::vector<long> grid_index(std::vector<double> const& momentum, std::vector<std::size_t> const& periods) {
    std::vector<long> index;
    for (std::size_t d = 0; d < momentum.size(); ++d)
        index.push_back(std::lround(momentum[d] * static_cast<double>(periods[d]) / two_pi)
                        % static_cast<long>(periods[d]));
    return index;
}

// Blocks keyed by momentum grid index and charges, eigenvalues ascending.
using BlockKey = std::pair<std::vector<long>, std::vector<double>>;
std::map<BlockKey, std::vector<double>> blocks_by_key(alps::ed::ModelSpectrum const& spectrum,
                                                      std::vector<std::size_t> const& periods) {
    std::map<BlockKey, std::vector<double>> result;
    for (auto const& block : spectrum.blocks) {
        BlockKey key{grid_index(block.momentum, periods), block.quantum_numbers};
        EXPECT_TRUE(result.emplace(key, block.eigenvalues).second);
    }
    return result;
}

// The sparsediag loop over one fixed charge sector: TOTAL_MOMENTUM set to
// every legacy translation momentum in turn, the Bloch-basis matrix of each
// diagonalized densely (sparsediag runs Lanczos on the same matrix).
std::map<std::vector<long>, std::vector<double>> legacy_momentum_levels(alps::Parameters p,
                                                                        std::vector<std::size_t> const& periods) {
    alps::graph_helper<> graph(p);
    std::map<std::vector<long>, std::vector<double>> result;
    for (auto const& k : graph.translation_momenta()) {
        p["TOTAL_MOMENTUM"] = alps::write_vector(k);
        alps::hamiltonian_matrix<alps::numeric::matrix<std::complex<double>>> legacy(p);
        auto const dimension = legacy.dimension();
        std::vector<double> momentum(k.begin(), k.end());
        auto& levels = result[grid_index(momentum, periods)];
        if (!dimension) continue;
        auto const& matrix = legacy.matrix();
        std::vector<std::complex<double>> dense(dimension * dimension);
        for (std::size_t column = 0; column < dimension; ++column)
            for (std::size_t row = 0; row < dimension; ++row) dense[column * dimension + row] = matrix(row, column);
        levels = alps::ed::dense_eigensystem(std::move(dense), dimension, false).eigenvalues;
    }
    return result;
}

// One catalog Heisenberg model in one Sz sector, every momentum: ALPS::ed on
// the catalog's own translations against the kit Hamiltonian and against the
// legacy Bloch matrices; leaves the catalog blocks by momentum index.
void check_heisenberg_momenta(Params const& values, lat::Lattice const& kit_lattice,
                              std::map<std::vector<long>, std::vector<double>>& by_momentum) {
    auto const p = parameters(values);
    mdl::CatalogModel<> model(p);
    auto const translations = spectrum_translations(model);
    std::vector<std::size_t> periods;
    for (auto const& t : translations) periods.push_back(t.period);

    alps::ed::ModelSpectrumOptions options;
    options.sector = model.sector();
    options.all_momenta = true;
    auto const catalog = alps::ed::model_spectrum(model.compile_hamiltonian(),
        alps::ed::projected_spaces(model.physical_basis()), translations, options);
    auto const kit = alps::ed::model_spectrum(mdl::heisenberg(kit_lattice), options);
    auto const catalog_blocks = blocks_by_key(catalog, periods), kit_blocks = blocks_by_key(kit, periods);
    std::size_t states = 1;
    for (auto const period : periods) states *= period;
    EXPECT_EQ(catalog_blocks.size(), states);  // every momentum of the one Sz sector
    EXPECT_EQ(catalog_blocks.size(), kit_blocks.size());

    auto const legacy = legacy_momentum_levels(p, periods);
    EXPECT_EQ(legacy.size(), catalog_blocks.size());
    for (auto const& [key, levels] : catalog_blocks) {
        SCOPED_TRACE("momentum index " + alps::write_vector(key.first));
        auto const reference = kit_blocks.find(key);
        ASSERT_NE(reference, kit_blocks.end());
        ASSERT_EQ(levels.size(), reference->second.size());
        for (std::size_t i = 0; i < levels.size(); ++i) EXPECT_NEAR(levels[i], reference->second[i], 1e-10);
        auto const old = legacy.find(key.first);
        ASSERT_NE(old, legacy.end());
        ASSERT_EQ(levels.size(), old->second.size());
        for (std::size_t i = 0; i < levels.size(); ++i) EXPECT_NEAR(levels[i], old->second[i], 1e-10);
        by_momentum[key.first] = levels;
    }
}

// Levels recorded from a sparsediag run, keyed by momentum grid index: the
// two lowest distinct eigenvalues per TOTAL_MOMENTUM (NUMBER_EIGENVALUES=2;
// the IETL Lanczos reports distinct values), read from the run's
// spectrum/sectors/*/energies. The runs used the parameters of the test
// below, built from this tree's src/apps/diag/sparsediag on 2026-10-18.
void expect_sparsediag_levels(std::map<std::vector<long>, std::vector<double>> const& levels,
                              std::map<std::vector<long>, std::vector<double>> const& recorded) {
    ASSERT_EQ(levels.size(), recorded.size());
    for (auto const& [index, values] : recorded) {
        SCOPED_TRACE("momentum index " + alps::write_vector(index));
        auto const found = levels.find(index);
        ASSERT_NE(found, levels.end());
        std::vector<double> distinct;
        for (double level : found->second)
            if (distinct.empty() || level - distinct.back() > 1e-8) distinct.push_back(level);
        ASSERT_GE(distinct.size(), values.size());
        for (std::size_t i = 0; i < values.size(); ++i) EXPECT_NEAR(distinct[i], values[i], 1e-9);
    }
}
} // namespace

TEST(CatalogTranslations, ARingTranslatesEverySiteOneStepForward) {
    mdl::CatalogModel<> model(parameters({{"MODEL", "spin"}, {"LATTICE", "chain lattice"}, {"L", "7"}}));
    auto const translations = model.translations();
    ASSERT_EQ(translations.size(), 1u);
    EXPECT_EQ(translations[0].direction, 0);
    EXPECT_EQ(translations[0].period, 7u);
    for (std::size_t site = 0; site < 7; ++site) EXPECT_EQ(translations[0].sites[site], (site + 1) % 7);
    expect_kit_generators(translations, lat::chain(7));
}

TEST(CatalogTranslations, OpenLatticesAndPlainGraphsHaveNone) {
    EXPECT_TRUE(mdl::CatalogModel<>(parameters({{"MODEL", "spin"}, {"LATTICE", "open chain lattice"}, {"L", "6"}}))
                    .translations().empty());
    EXPECT_TRUE(mdl::CatalogModel<>(parameters({{"MODEL", "spin"}, {"LATTICE", "open square lattice"}, {"L", "3"}}))
                    .translations().empty());
    // A <GRAPH> carries no lattice at all.
    EXPECT_TRUE(mdl::CatalogModel<>(parameters({{"MODEL", "spin"}, {"GRAPH", "triangle"}})).translations().empty());
}

TEST(CatalogTranslations, CatalogGraphsMatchTheKitGenerators) {
    // Periodic along the first axis only: one generator, labelled axis 0.
    auto const strip = mdl::CatalogModel<>(parameters({{"MODEL", "spin"}, {"LATTICE", "square lattice 3x3"}})).translations();
    ASSERT_EQ(strip.size(), 1u);
    EXPECT_EQ(strip[0].direction, 0);
    EXPECT_EQ(strip[0].period, 3u);
    for (std::size_t x = 0; x < 3; ++x)
        for (std::size_t y = 0; y < 3; ++y) EXPECT_EQ(strip[0].sites[x * 3 + y], ((x + 1) % 3) * 3 + y);

    expect_kit_generators(
        mdl::CatalogModel<>(parameters({{"MODEL", "spin"}, {"LATTICE", "square lattice"}, {"L", "4"}, {"W", "3"}}))
            .translations(),
        lat::square(4, 3));
    expect_kit_generators(
        mdl::CatalogModel<>(parameters({{"MODEL", "spin"}, {"LATTICE", "ladder"}, {"L", "5"}})).translations(),
        lat::ladder(5));
}

TEST(CatalogTranslations, KitLatticesKeepTheirGeneratorsThroughTheLegacyReader) {
    auto const spin = parameters({{"MODEL", "spin"}, {"J", "1"}});
    // square(3, 1) is periodic with extent 1 along the second axis: an identity generator of period 1.
    std::vector<lat::Lattice> lattices = {lat::square(4, 3), lat::honeycomb(3, 2), lat::kagome(2, 3),
                                          lat::simple_cubic(2, 3, 2), lat::ladder(4), lat::square(3, 1)};
    auto const base = lat::square(4, 3);
    // A cylinder, periodic along the second axis only.
    lattices.push_back(lat::Lattice::build(base.unit_cell(), base.basis_vectors(), base.extent(),
                                           {lat::Boundary::open, lat::Boundary::periodic}));
    for (auto const& lattice : lattices) {
        SCOPED_TRACE(lattice.num_sites());
        mdl::CatalogModel<> model(lattice, spin);
        ASSERT_EQ(model.graph().num_sites(), static_cast<std::size_t>(lattice.num_sites()));
        expect_kit_generators(model.translations(), lattice);
    }
}

TEST(CatalogTranslations, HeisenbergChainMomentaMatchTheKitAndSparsediag) {
    std::map<std::vector<long>, std::vector<double>> levels;
    check_heisenberg_momenta({{"MODEL", "spin"}, {"LATTICE", "chain lattice"}, {"L", "10"}, {"J", "1"},
                              {"local_S", "0.5"}, {"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Sz_total", "0"}},
                             lat::chain(10), levels);
    expect_sparsediag_levels(levels, {
        {{0}, {-4.092207346739, -3.770597435408}},
        {{1}, {-3.246164916733, -2.886286924510}},
        {{2}, {-2.780490100408, -2.514299479211}},
        {{3}, {-2.891236100351, -2.767634884836}},
        {{4}, {-3.543279374313, -2.320334627141}},
        {{5}, {-4.515446354492, -2.975931869054}},
        {{6}, {-3.543279374313, -2.320334627141}},
        {{7}, {-2.891236100351, -2.767634884836}},
        {{8}, {-2.780490100408, -2.514299479211}},
        {{9}, {-3.246164916733, -2.886286924510}},
    });
}

TEST(CatalogTranslations, HeisenbergSquareMomentaMatchTheKitAndSparsediag) {
    std::map<std::vector<long>, std::vector<double>> levels;
    check_heisenberg_momenta({{"MODEL", "spin"}, {"LATTICE", "square lattice"}, {"L", "4"}, {"W", "4"}, {"J", "1"},
                              {"local_S", "0.5"}, {"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Sz_total", "0"}},
                             lat::square(4, 4), levels);
    expect_sparsediag_levels(levels, {
        {{0, 0}, {-11.228483208429, -9.517687983937}},
        {{0, 1}, {-8.794378908138, -7.877905391372}},
        {{0, 2}, {-8.518283596238, -8.252826957308}},
        {{0, 3}, {-8.794378908137, -7.877905391372}},
        {{1, 0}, {-8.794378908138, -7.877905391372}},
        {{1, 1}, {-8.518283596238, -7.812148008952}},
        {{1, 2}, {-8.886442407806, -7.877905391372}},
        {{1, 3}, {-8.518283596238, -7.812148008952}},
        {{2, 0}, {-8.518283596238, -8.252826957308}},
        {{2, 1}, {-8.886442407806, -7.877905391372}},
        {{2, 2}, {-10.649884872663, -7.836537573480}},
        {{2, 3}, {-8.886442407806, -7.877905391372}},
        {{3, 0}, {-8.794378908137, -7.877905391372}},
        {{3, 1}, {-8.518283596238, -7.812148008952}},
        {{3, 2}, {-8.886442407806, -7.877905391372}},
        {{3, 3}, {-8.518283596238, -7.812148008952}},
    });
}
