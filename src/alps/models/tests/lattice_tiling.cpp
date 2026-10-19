// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The tiling contract of the alps::lattice kit (lattice.hpp, README
// "Conventions") on every named factory and two custom cells. The cases and
// the custom cells are those of skilledwolf/alpsx (private repository),
// branch agent/v3-execution-and-api-audit at f8ee6e98b,
// src/model/tests/lattice/lattice_kit_parity.cpp, which compared them with
// frozen fixtures from alpsx's Python lattice module; here the contract is
// checked directly instead:
//
//   - site = cell * num_vertices + vertex, cells last dimension fastest, at
//     sum_d (n_d + f_v[d]) a_d with the vertex's type;
//   - one bond per (cell, unit-cell edge) whose endpoints exist, self-bonds
//     dropped, `wraps` iff the periodic modulo moved an endpoint cell, and
//     `vector` the in-cluster coordinate difference: the tiling is redone
//     here from unit_cell(), basis_vectors(), extent() and boundary();
//   - momenta(): k_d = 2 pi j / L_d over the periodic directions;
//   - the unit cells: the coordination of the all-periodic cases and, on
//     the isotropic ones, one length per bond type.
//
// The legacy catalog checks the catalog-named lattices' sites, types,
// coordinates and bonds independently (src/alps/model/tests/lattice_xml_bridge.cpp).
#include <alps/models/lattice.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace lat = alps::lattice;

namespace {

lat::Lattice custom_bilayer(int n) {
    lat::UnitCell uc(2);
    uc.add_vertex({0.0, 0.0}, 0);
    uc.add_vertex({0.0, 0.0}, 1);
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    uc.add_edge(1, 1, 1, {}, {1, 0});
    uc.add_edge(1, 1, 1, {}, {0, 1});
    uc.add_edge(0, 1, 2);
    return uc.tile({n, n});
}

lat::Lattice custom_cylinder(int l, int w) {
    lat::UnitCell uc(2);
    uc.add_vertex({0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    return uc.tile({l, w}, {}, {lat::Boundary::periodic, lat::Boundary::open});
}

struct Case {
    std::function<lat::Lattice()> factory;
    int coordination; // bonds at every site; -1 where a direction is open
    bool isotropic;   // one length per bond type (the anisotropic Bravais
                      // lattices give every axis bond type 0)
};

std::map<std::string, Case> const& registry() {
    using B = lat::Boundary;
    static std::map<std::string, Case> const r = {
        {"chain_L6_p", {[] { return lat::chain(6); }, 2, true}},
        {"chain_L6_o", {[] { return lat::chain(6, 1.0, B::open); }, -1, true}},
        {"chain_L2_p", {[] { return lat::chain(2); }, 2, true}},
        {"chain_L1_p", {[] { return lat::chain(1); }, 0, true}},
        {"nnn_chain_L6_p", {[] { return lat::nnn_chain(6); }, 4, true}},
        {"nnn_chain_L5_o", {[] { return lat::nnn_chain(5, 1.0, B::open); }, -1, true}},
        {"square_3x4_p", {[] { return lat::square(3, 4); }, 4, true}},
        {"square_3x4_o", {[] { return lat::square(3, 4, 1.0, B::open); }, -1, true}},
        {"frustrated_square_3x3_p", {[] { return lat::frustrated_square(3); }, 8, true}},
        {"triangular_3x3_p", {[] { return lat::triangular(3); }, 6, true}},
        {"triangular_3x3_o", {[] { return lat::triangular(3, 0, 1.0, B::open); }, -1, true}},
        {"honeycomb_2x3_p", {[] { return lat::honeycomb(2, 3); }, 3, true}},
        {"honeycomb_2x2_o", {[] { return lat::honeycomb(2, 2, 1.0, B::open); }, -1, true}},
        {"kagome_2x2_p", {[] { return lat::kagome(2); }, 4, true}},
        {"simple_cubic_2x3x2_p", {[] { return lat::simple_cubic(2, 3, 2); }, 6, true}},
        {"simple_cubic_2x2x2_o", {[] { return lat::simple_cubic(2, 0, 0, 1.0, B::open); }, -1, true}},
        {"ladder_L4_p", {[] { return lat::ladder(4); }, 3, true}},
        {"ladder_L4_o", {[] { return lat::ladder(4, 1.0, B::open); }, -1, true}},
        {"rectangular_3x2_b05_p", {[] { return lat::rectangular(3, 2, 1.0, 0.5); }, 4, false}},
        {"fcc_2x2x2_p", {[] { return lat::face_centered_cubic(2); }, 12, true}},
        {"bcc_2x2x2_p", {[] { return lat::body_centered_cubic(2); }, 8, true}},
        {"tetragonal_2x2x2_c15_p", {[] { return lat::tetragonal(2, 0, 0, 1.0, 1.5); }, 6, false}},
        {"orthorhombic_2x2x2_abc_p", {[] { return lat::orthorhombic(2, 0, 0, 1.0, 1.1, 1.2); }, 6, false}},
        {"hexagonal_2x2x2_p", {[] { return lat::hexagonal(2); }, 8, true}},
        {"custom_bilayer_3x3_p", {[] { return custom_bilayer(3); }, 5, true}},
        {"custom_cylinder_4x3", {[] { return custom_cylinder(4, 3); }, -1, true}},
    };
    return r;
}

constexpr double kTol = 1e-12;

std::vector<int> decode(int cell, std::vector<int> const& extent) {
    std::vector<int> n(extent.size());
    for (std::size_t d = extent.size(); d-- > 0;) {
        n[d] = cell % extent[d];
        cell /= extent[d];
    }
    return n;
}

int encode(std::vector<int> const& n, std::vector<int> const& extent) {
    int cell = 0;
    for (std::size_t d = 0; d < extent.size(); ++d) cell = cell * extent[d] + n[d];
    return cell;
}

std::vector<double> position(lat::Lattice const& lattice, std::vector<int> const& n, int vertex) {
    auto const& basis = lattice.basis_vectors();
    auto const& frac = lattice.unit_cell().vertex_coords().at(static_cast<std::size_t>(vertex));
    std::vector<double> r(basis.at(0).size(), 0.0);
    for (std::size_t d = 0; d < basis.size(); ++d)
        for (std::size_t w = 0; w < r.size(); ++w) r[w] += (n[d] + frac[d]) * basis[d][w];
    return r;
}

double distance(std::vector<double> const& a, std::vector<double> const& b) {
    double sum = 0;
    for (std::size_t w = 0; w < a.size(); ++w) sum += (a[w] - b[w]) * (a[w] - b[w]);
    return std::sqrt(sum);
}

using BondKey = std::tuple<int, int, int, bool>; // source, target, type, wraps

// The bonds the contract prescribes, and the length of each bond type's
// edges in the infinite lattice.
std::vector<BondKey> tile(lat::Lattice const& lattice, std::map<int, std::vector<double>>& lengths) {
    auto const& cell = lattice.unit_cell();
    auto const& extent = lattice.extent();
    int const nv = cell.num_vertices();
    int cells = 1;
    for (int e : extent) cells *= e;
    std::vector<BondKey> bonds;
    for (int id = 0; id < cells; ++id) {
        auto const n = decode(id, extent);
        for (auto const& edge : cell.edges()) {
            auto offset = [&](std::vector<int> const& o, std::size_t d) { return o.empty() ? 0 : o[d]; };
            std::vector<int> s(n.size()), t(n.size());
            bool wraps = false, outside = false;
            for (std::size_t d = 0; d < n.size(); ++d) {
                s[d] = n[d] + offset(edge.source_offset, d);
                t[d] = n[d] + offset(edge.target_offset, d);
                for (int* x : {&s[d], &t[d]}) {
                    if (*x >= 0 && *x < extent[d]) continue;
                    if (lattice.boundary()[d] == lat::Boundary::open) { outside = true; continue; }
                    wraps = true;
                    *x = ((*x % extent[d]) + extent[d]) % extent[d];
                }
            }
            if (id == 0) {
                std::vector<int> s0(n.size()), t0(n.size());
                for (std::size_t d = 0; d < n.size(); ++d) {
                    s0[d] = offset(edge.source_offset, d);
                    t0[d] = offset(edge.target_offset, d);
                }
                lengths[edge.edge_type].push_back(
                    distance(position(lattice, t0, edge.target), position(lattice, s0, edge.source)));
            }
            if (outside) continue;
            int const source = encode(s, extent) * nv + edge.source, target = encode(t, extent) * nv + edge.target;
            if (source != target) bonds.emplace_back(source, target, edge.edge_type, wraps);
        }
    }
    return bonds;
}

void check(std::string const& name, Case const& item) {
    SCOPED_TRACE(name);
    auto const lattice = item.factory();
    auto const& extent = lattice.extent();
    int const nv = lattice.unit_cell().num_vertices();

    int cells = 1;
    for (int e : extent) cells *= e;
    ASSERT_EQ(lattice.num_sites(), cells * nv);
    for (auto const& site : lattice.sites()) {
        SCOPED_TRACE("site " + std::to_string(site.index));
        int const vertex = site.index % nv;
        EXPECT_EQ(site.index, &site - lattice.sites().data());
        EXPECT_EQ(site.sublattice, vertex);
        EXPECT_EQ(site.site_type, lattice.unit_cell().vertex_types().at(static_cast<std::size_t>(vertex)));
        EXPECT_LT(distance(site.coordinate, position(lattice, decode(site.index / nv, extent), vertex)), kTol);
    }

    std::map<int, std::vector<double>> lengths;
    auto expected = tile(lattice, lengths);
    std::vector<BondKey> actual;
    std::vector<int> degree(static_cast<std::size_t>(lattice.num_sites()), 0);
    for (auto const& bond : lattice.bonds()) {
        actual.emplace_back(bond.source, bond.target, bond.bond_type, bond.wraps);
        ++degree.at(static_cast<std::size_t>(bond.source));
        ++degree.at(static_cast<std::size_t>(bond.target));
        auto const& rs = lattice.sites().at(static_cast<std::size_t>(bond.source)).coordinate;
        auto const& rt = lattice.sites().at(static_cast<std::size_t>(bond.target)).coordinate;
        ASSERT_EQ(bond.vector.size(), rs.size());
        for (std::size_t w = 0; w < rs.size(); ++w)
            EXPECT_NEAR(bond.vector[w], rt[w] - rs[w], kTol) << "bond " << bond.source << "-" << bond.target;
    }
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
    EXPECT_EQ(actual, expected);

    if (item.coordination >= 0) {
        for (std::size_t i = 0; i < degree.size(); ++i) EXPECT_EQ(degree[i], item.coordination) << "site " << i;
    }
    if (item.isotropic) {
        for (auto const& [type, values] : lengths)
            for (double length : values) EXPECT_NEAR(length, values.front(), kTol) << "bond type " << type;
    }

    auto const directions = lattice.periodic_directions();
    auto const momenta = lattice.momenta();
    if (directions.empty()) {
        ASSERT_EQ(momenta.size(), 1u);
        EXPECT_TRUE(momenta.front().empty());
        return;
    }
    std::size_t count = 1;
    for (int d : directions) count *= static_cast<std::size_t>(extent[static_cast<std::size_t>(d)]);
    ASSERT_EQ(momenta.size(), count);
    ASSERT_EQ(lattice.num_momenta(), static_cast<int>(count));
    std::vector<std::vector<int>> indices;
    for (auto const& k : momenta) {
        ASSERT_EQ(k.size(), directions.size());
        std::vector<int> j;
        for (std::size_t c = 0; c < k.size(); ++c) {
            double const L = extent[static_cast<std::size_t>(directions[c])];
            double const x = k[c] * L / (2 * M_PI);
            EXPECT_NEAR(x, std::round(x), 1e-12);
            j.push_back(static_cast<int>(std::lround(x)));
            EXPECT_GE(j.back(), 0);
            EXPECT_LT(j.back(), L);
        }
        indices.push_back(j);
    }
    std::sort(indices.begin(), indices.end());
    EXPECT_EQ(std::adjacent_find(indices.begin(), indices.end()), indices.end()) << "repeated momentum";
}

} // namespace

TEST(LatticeTiling, FactoriesFollowTheContract) {
    for (auto const& [name, item] : registry()) check(name, item);
    EXPECT_EQ(registry().size(), 26u);
}
