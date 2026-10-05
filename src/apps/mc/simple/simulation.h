// Copyright (C) 1997–2015 Synge Todo; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include "../moment_difference.hpp"

#include <alps/mcbase.hpp>
#include <alps/lattice.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/alea/checkpoint.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/alea/transform.hpp>
#include <alps/alea/transformer.hpp>
#include <alps/hdf5/stdarray.hpp>
#include <boost/filesystem/operations.hpp>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>

namespace simplemc {

// One local-update engine for H = -sum_b J_b s_i.s_j - H sum_i s_axis.
// The field axis is X for Ising/XY and Z for Heisenberg.
class simulation : public alps::mcbase, private alps::graph_helper<> {
    using spin = std::array<double, 3>;
    static constexpr double pi = 3.14159265358979323846;
public:
    simulation(alps::params const& p, std::size_t bins = 128, std::size_t chain = 0)
        : mcbase(p, chain), graph_helper<>(graph_parameters(p)), bins_(bins), chain_(chain),
          model_(p["ALGORITHM"].as<std::string>()),
          dimensions_(model_ == "ising" ? 1 : model_ == "xy" ? 2 : 3),
          beta_(inverse_temperature(p)), field_(finite(p.value_or("H", 0.), "H")),
          thermalization_(p["THERMALIZATION"].as<uint64_t>()),
          production_(p["SWEEPS"].as<uint64_t>()), spins_(num_sites()) {
        if (model_ != "ising" && model_ != "xy" && model_ != "heisenberg")
            throw std::invalid_argument("simplemc requires ising, xy or heisenberg");
        if (!num_sites() || !production_ || thermalization_ > UINT64_MAX - production_)
            throw std::invalid_argument("simplemc requires sites and a representable positive sweep count");
        double energy_bound = num_sites() * std::abs(field_);
        for (auto [it, end] = bonds(); it != end; ++it) {
            auto type = bond_type(*it);
            auto j = finite(p.value_or("J" + std::to_string(type), p.value_or("J", 1.)), "J");
            coupling_.emplace(type, j);
            energy_bound += std::abs(j);
        }
        // E^2 errors square batch sums, so bound their fourth moments too.
        if (!(energy_bound <= std::sqrt(std::sqrt(std::numeric_limits<double>::max())) /
                                 (4 * std::sqrt(double(production_)))))
            throw std::invalid_argument("simplemc energy moments would overflow");
        for (auto& s : spins_) s = proposal();
        for (auto const& name : {"Number of Sites", "Energy", "Energy Density", "Energy^2",
                                 "Magnetization Density^2", "Magnetization Density^4"})
            add_measurement(name);
        add_measurement(projection());
        if (dimensions_ > 1) {
            add_measurement(projection() + "^2");
            add_measurement(projection() + "^4");
        }
    }
    simulation(simulation const&) = delete;
    simulation& operator=(simulation const&) = delete;
    using mcbase::save;
    using mcbase::load;
    static alps::params const& checkpoint_parameters(alps::params const& p) { return p; }

    void update() override {
        if (sweeps_ == thermalization_ + production_)
            throw std::logic_error("simplemc run is complete");
        if (measurement_count() != (sweeps_ > thermalization_ ? sweeps_ - thermalization_ : 0))
            throw std::logic_error("measure the previous simplemc sweep before updating");
        for (std::size_t site = 0; site != spins_.size(); ++site) {
            auto proposed = dimensions_ == 1 ? spin{-spins_[site][0], 0, 0} : proposal();
            double difference = field_ * (spins_[site][axis()] - proposed[axis()]);
            for (auto [it, end] = neighbor_bonds(site); it != end; ++it) {
                auto neighbor = source(*it) == site ? target(*it) : source(*it);
                // A self-loop contributes constant -J to the Hamiltonian.
                if (neighbor != site)
                    difference += coupling_.at(bond_type(*it)) *
                        (dot(spins_[site], spins_[neighbor]) - dot(proposed, spins_[neighbor]));
            }
            bool accept = dimensions_ == 1
                ? random() < .5 * (1 + std::tanh(-.5 * beta_ * difference))
                : difference <= 0 || random() < std::exp(-beta_ * difference);
            if (accept) spins_[site] = proposed;
        }
        ++sweeps_;
    }

    void measure() override {
        if (sweeps_ <= thermalization_) return;
        if (measurement_count() != sweeps_ - thermalization_ - 1)
            throw std::logic_error("measure each simplemc sweep exactly once");
        double energy = 0;
        for (auto [it, end] = bonds(); it != end; ++it)
            energy -= coupling_.at(bond_type(*it)) * dot(spins_[source(*it)], spins_[target(*it)]);
        spin magnetization{};
        for (auto const& s : spins_) {
            energy -= field_ * s[axis()];
            for (int d = 0; d != 3; ++d) magnetization[d] += s[d] / spins_.size();
        }
        double m2 = dot(magnetization, magnetization), component = magnetization[axis()];
        record("Number of Sites", spins_.size());
        record("Energy", energy);
        record("Energy Density", energy / spins_.size());
        record("Energy^2", energy * energy);
        record("Magnetization Density^2", m2);
        record("Magnetization Density^4", m2 * m2);
        record(projection(), component);
        if (dimensions_ > 1) {
            record(projection() + "^2", component * component);
            record(projection() + "^4", std::pow(component, 4));
        }
    }

    double fraction_completed() const override { return double(measurement_count()) / production_; }
    uint64_t completed_sweeps() const { return sweeps_; }
    uint64_t measurement_count() const { return measurements.at("Energy")->count(); }
    std::size_t chain_id() const { return chain_; }

    void save(alps::hdf5::archive& ar) const override {
        if (measurement_count() != (sweeps_ > thermalization_ ? sweeps_ - thermalization_ : 0))
            throw std::logic_error("measure simplemc before checkpointing");
        mcbase::save(ar);
        ar["checkpoint/sweeps"] << sweeps_;
        ar["checkpoint/chain_id"] << uint64_t(chain_);
        ar["checkpoint/spins"] << spins_;
        ar["checkpoint/bonds"] << topology();
        ar["checkpoint/couplings"] << bond_couplings();
    }

    void load(alps::hdf5::archive& ar) override {
        alps::params restored_parameters;
        uint64_t sweeps, chain;
        std::vector<spin> spins;
        std::vector<std::array<uint64_t, 3>> graph;
        std::vector<double> couplings;
        ar["/parameters"] >> restored_parameters;
        ar["checkpoint/sweeps"] >> sweeps;
        ar["checkpoint/chain_id"] >> chain;
        ar["checkpoint/spins"] >> spins;
        ar["checkpoint/bonds"] >> graph;
        ar["checkpoint/couplings"] >> couplings;
        if (restored_parameters != parameters || chain != chain_ || graph != topology()
                || couplings != bond_couplings() || sweeps > thermalization_ + production_
                || spins.size() != spins_.size())
            throw std::invalid_argument("simplemc checkpoint does not match this run");
        for (auto const& s : spins) {
            if (!std::isfinite(dot(s, s)) || std::abs(dot(s, s) - 1.) > 1e-12
                    || (dimensions_ == 1 && std::abs(s[0]) != 1.)
                    || (dimensions_ < 3 && s[2] != 0.) || (dimensions_ == 1 && s[1] != 0.))
                throw std::invalid_argument("invalid simplemc checkpoint spin");
        }
        // The shared codec checks the complete native layout before commit.
        alps::alea::hdf5_serializer serializer(ar, "measurements");
        if (ar.list_children("measurements").size() != measurements.size())
            throw std::invalid_argument("unexpected simplemc checkpoint measurements");
        for (auto const& entry : measurements) {
            auto name = ar.encode_segment(entry.first);
            alps::alea::batch_acc<double> value;
            alps::alea::deserialize(serializer, name, value);
            if (value.size() != 1 || value.num_batches() != bins_
                    || value.current_batch_size() != value.cursor().factor()
                    || value.count() != (sweeps > thermalization_ ? sweeps - thermalization_ : 0)
                    || !value.store().batch().allFinite())
                throw std::invalid_argument("invalid simplemc checkpoint measurement state");
        }
        mcbase::load(ar);
        spins_ = std::move(spins);
        sweeps_ = sweeps;
    }

    static results_type derive(results_type const& raw, alps::params const& p) {
        auto results = raw;
        for (auto const& entry : raw) validate_result(entry.second);
        if (raw.empty() || (raw.at("Energy").store().count().array() > 0).count() < 2) return results;
        double beta = inverse_temperature(p), n = raw.at("Number of Sites").mean()(0);
        auto difference = native_mc::moment_difference(raw.at("Energy").count(), raw.at("Energy").num_batches());
        auto transform = [&](std::string const& name, std::string const& first, std::string const& second,
                             std::function<double(double, double)> function) {
            auto inputs = alps::alea::join(raw.at(first), raw.at(second));
            auto result = alps::alea::transform(alps::alea::jackknife_prop{},
                alps::alea::make_transformer<double>(std::move(function)), inputs);
            if (result.store().batch().allFinite()) {
                validate_result(result);
                results.emplace(name, std::move(result));
            }
        };
        transform("Specific Heat", "Energy", "Energy^2",
                  [=](double e, double e2) { return beta == 0 ? 0. : beta * (beta * difference(e2, e * e)) / n; });
        for (auto const& suffix : {"", " X", " Z"}) {
            auto first = std::string("Magnetization Density") + suffix + "^2", second = first.substr(0, first.size()-1) + "4";
            auto it = raw.find(second);
            // The full mean and every leave-one-out denominator must exist.
            if (it != raw.end() && it->second.store().batch().sum() > it->second.store().batch().maxCoeff())
                transform(std::string("Binder Ratio of Magnetization") + suffix, first, second,
                          [](double m2, double m4) { return m2 * m2 / m4; });
        }
        return results;
    }

    void snapshot(std::filesystem::path const& filename) const {
        spin center{};
        for (std::size_t site = 0; site != spins_.size(); ++site) {
            if (coordinate(site).size() > 3)
                throw std::invalid_argument("VTK snapshots require at most three spatial dimensions");
            for (std::size_t d = 0; d != coordinate(site).size(); ++d)
                center[d] += coordinate(site)[d] / spins_.size();
        }
        auto target = std::filesystem::absolute(filename), temporary = target.parent_path();
        do { temporary = target.parent_path() / boost::filesystem::unique_path(".simplemc-%%%%-%%%%-%%%%").string(); }
        while (!std::filesystem::create_directory(temporary));
        try {
            std::ofstream out(temporary / "snapshot.vtk");
            out.exceptions(std::ios::failbit | std::ios::badbit);
            out << std::setprecision(17) << "# vtk DataFile Version 3.0\nsimplemc snapshot\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
                << spins_.size() << " double\n";
            for (std::size_t site = 0; site != spins_.size(); ++site) {
                auto const& point = coordinate(site);
                for (int d = 0; d != 3; ++d) out << (d < point.size() ? point[d] - center[d] : 0.) << (d == 2 ? '\n' : ' ');
            }
            out << "CELLS 0 0\nCELL_TYPES 0\nPOINT_DATA " << spins_.size() << '\n';
            if (dimensions_ == 1) {
                out << "SCALARS spins double\nLOOKUP_TABLE default\n";
                for (auto const& s : spins_) out << s[0] << '\n';
            } else {
                out << "VECTORS spins double\n";
                for (auto const& s : spins_) out << s[0] << ' ' << s[1] << ' ' << s[2] << '\n';
                out << "SCALARS angle double\nLOOKUP_TABLE default\n";
                for (auto const& s : spins_) {
                    double angle = std::atan2(s[1], s[0]) / (2 * pi);
                    out << (angle < 0 ? angle + 1 : angle) << '\n';
                }
            }
            out.close();
            // An atomic link publishes the complete file and refuses collisions.
            std::filesystem::create_hard_link(temporary / "snapshot.vtk", target);
        } catch (...) { std::filesystem::remove_all(temporary); throw; }
        std::filesystem::remove_all(temporary);
    }

private:
    static alps::Parameters graph_parameters(alps::params const& p) {
        auto seed = p.value_or<uint32_t>("DISORDER_SEED", uint32_t(p.value_or("SEED", 42)));
        alps::Disorder::seed(seed);
        return alps::make_deprecated_parameters(p);
    }
    static void validate_result(alps::alea::batch_result<double> const& result) {
        if (!result.valid() || !result.store().batch().allFinite()
                || (result.count() && !result.mean().allFinite())
                || ((result.store().count().array() > 0).count() > 1 && !result.stderror().allFinite()))
            throw std::overflow_error("simplemc result moments or errors are not representable");
    }
    static double finite(double value, char const* name) {
        if (!std::isfinite(value)) throw std::invalid_argument(std::string("nonfinite simplemc ") + name);
        return value;
    }
    static double inverse_temperature(alps::params const& p) {
        if (!p.exists("T")) return 0;
        auto temperature = finite(p["T"].as<double>(), "T");
        if (!(temperature > 0) || !std::isfinite(1. / temperature))
            throw std::invalid_argument("simplemc T must be positive with finite inverse");
        return 1. / temperature;
    }
    static double dot(spin const& a, spin const& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
    int axis() const { return dimensions_ == 3 ? 2 : 0; }
    std::string projection() const {
        return std::string("Magnetization Density") + (dimensions_ == 1 ? "" : dimensions_ == 2 ? " X" : " Z");
    }
    spin proposal() {
        if (dimensions_ == 1) return {random() < .5 ? -1. : 1., 0, 0};
        double angle = 2 * pi * random(), z = dimensions_ == 3 ? 2 * random() - 1 : 0;
        double radius = std::sqrt(std::max(0., 1 - z*z));
        return {radius * std::cos(angle), radius * std::sin(angle), z};
    }
    void add_measurement(std::string const& name) {
        measurements.emplace(name, std::make_shared<alps::alea::batch_acc<double>>(1, bins_));
    }
    void record(std::string const& name, double value) {
        *measurements.at(name) << alps::alea::make_adapter(value);
    }
    std::vector<std::array<uint64_t, 3>> topology() const {
        std::vector<std::array<uint64_t, 3>> result;
        for (auto [it, end] = bonds(); it != end; ++it)
            result.push_back({uint64_t(source(*it)), uint64_t(target(*it)), uint64_t(bond_type(*it))});
        return result;
    }
    std::vector<double> bond_couplings() const {
        std::vector<double> result;
        for (auto [it, end] = bonds(); it != end; ++it) result.push_back(coupling_.at(bond_type(*it)));
        return result;
    }
    std::size_t bins_, chain_;
    std::string model_;
    int dimensions_;
    double beta_, field_;
    uint64_t thermalization_, production_, sweeps_ = 0;
    std::map<alps::type_type, double> coupling_;
    std::vector<spin> spins_;
};
}
