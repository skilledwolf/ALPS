// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "../spinsim.h"
#include <boost/filesystem/operations.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace {
using vector = std::vector<double>;
using result = alps::alea::batch_result<double>;
constexpr double pi = 3.14159265358979323846;
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F function, char const* message) {
    bool rejected = false;
    try { function(); } catch (std::exception const&) { rejected = true; }
    require(rejected, message);
}
void agree(result const& measured, double expected, char const* message) {
    if (!(std::abs(measured.mean()(0) - expected) <= 6 * measured.stderror()(0) + 1e-3)) {
        std::cerr << message << ": " << measured.mean()(0) << " +/- " << measured.stderror()(0)
                  << ", expected " << expected << '\n';
        throw std::runtime_error(message);
    }
}
void step(spinmc::simulation& simulation, uint64_t count) {
    while (count--) { simulation.update(); simulation.measure(); }
}
void finish(spinmc::simulation& simulation) {
    for (uint64_t i = 0; simulation.fraction_completed() < 1; ++i) {
        require(i < 2000000, "simulation did not finish");
        step(simulation, 1);
    }
}
alps::params parameters(std::filesystem::path const& library, std::string const& model,
                        std::string const& graph = "pair", std::string const& update = "local") {
    alps::params p;
    p["LATTICE_LIBRARY"] = library.string(); p["GRAPH"] = graph; p["MODEL"] = model;
    p["SEED"] = 731; p["DISORDER_SEED"] = 42; p["THERMALIZATION"] = 300;
    p["SWEEPS"] = 50000; p["T"] = 1.1; p["UPDATE"] = update; p["J"] = vector{.7};
    return p;
}
auto sample(alps::params const& p) {
    spinmc::simulation simulation(p, 64);
    finish(simulation);
    require(simulation.measurement_count() == p["SWEEPS"].as<uint64_t>(), "incorrect production sample count");
    return simulation.collect_results();
}
vector diagonal(vector const& values) {
    vector matrix(values.size() * values.size());
    for (size_t i = 0; i != values.size(); ++i) matrix[i * values.size() + i] = values[i];
    return matrix;
}
vector packed(vector const& matrix, size_t dimensions) {
    vector values;
    for (size_t i = 0; i != dimensions; ++i)
        for (size_t j = 0; j <= i; ++j) values.push_back(matrix[i * dimensions + j]);
    return values;
}
double bilinear(vector const& x, vector const& matrix, vector const& y) {
    double value = 0;
    for (size_t i = 0; i != x.size(); ++i)
        for (size_t j = 0; j != y.size(); ++j) value += x[i] * matrix[i * y.size() + j] * y[j];
    return value;
}
struct point { vector spin; double weight; };
std::vector<point> sphere(size_t dimensions) {
    std::vector<point> points;
    if (dimensions == 2) {
        for (size_t i = 0; i != 96; ++i) {
            double angle = 2 * pi * (i + .5) / 96;
            points.push_back({{std::cos(angle), std::sin(angle)}, 1. / 96});
        }
    } else {
        auto const& nodes = boost::math::quadrature::gauss<double, 12>::abscissa();
        auto const& weights = boost::math::quadrature::gauss<double, 12>::weights();
        for (size_t z = 0; z != nodes.size(); ++z)
            for (double sign : {-1., 1.})
                for (size_t i = 0; i != 24; ++i) {
                    double angle = 2 * pi * (i + .5) / 24, height = sign * nodes[z];
                    double radius = std::sqrt(1 - height * height);
                    points.push_back({{radius * std::cos(angle), radius * std::sin(angle), height}, weights[z] / 48});
                }
    }
    return points;
}
// Numerical integration uses the Hamiltonian directly, independently of the
// engine's matrix decoder, local deltas and cluster estimators.
struct moments { double partition = 0, energy = 0, energy2 = 0, magnetization2 = 0, staggered2 = 0; };
moments continuous(size_t dimensions, vector const& j, vector const& d, vector const& h,
                   double s0, double s1, double temperature, bool loop = false) {
    auto points = sphere(dimensions);
    moments sum;
    for (auto const& a : points)
        for (size_t i = 0; i != (loop ? 1 : points.size()); ++i) {
            auto const& b = loop ? a : points[i];
            double energy = -s0 * (loop ? s0 : s1) * bilinear(a.spin, j, b.spin)
                          - s0 * s0 * bilinear(a.spin, d, a.spin);
            if (!loop) energy -= s1 * s1 * bilinear(b.spin, d, b.spin);
            double m2 = 0, staggered = 0;
            for (size_t i = 0; i != dimensions; ++i) {
                energy -= h[i] * (s0 * a.spin[i] + (loop ? 0 : s1 * b.spin[i]));
                m2 += std::pow(s0 * a.spin[i] + (loop ? 0 : s1 * b.spin[i]), 2) / (loop ? 1 : 4);
                staggered += std::pow(s0 * a.spin[i] - (loop ? 0 : s1 * b.spin[i]), 2) / (loop ? 1 : 4);
            }
            double weight = a.weight * (loop ? 1 : b.weight) * std::exp(-energy / temperature);
            sum.partition += weight; sum.energy += weight * energy; sum.energy2 += weight * energy * energy;
            sum.magnetization2 += weight * m2; sum.staggered2 += weight * staggered;
        }
    return sum;
}
void physical_moments(alps::params const& p, moments const& expected) {
    auto measured = sample(p);
    agree(measured.at("Energy"), expected.energy / expected.partition, "matrix Hamiltonian energy is incorrect");
    agree(measured.at("Energy^2"), expected.energy2 / expected.partition, "matrix Hamiltonian second moment is incorrect");
    agree(measured.at("Magnetization^2"), expected.magnetization2 / expected.partition, "spin normalization is incorrect");
}
void equal(spinmc::simulation const& first, spinmc::simulation const& second) {
    require(first.completed_sweeps() == second.completed_sweeps()
            && first.measurement_count() == second.measurement_count()
            && first.fraction_completed() == second.fraction_completed(), "restart changed counters");
    auto a = first.collect_results(), b = second.collect_results();
    require(a.size() == b.size(), "restart changed measurement names");
    for (auto const& entry : a) {
        auto const& other = b.at(entry.first);
        require(entry.second.store().batch() == other.store().batch()
                && entry.second.store().count() == other.store().count(), "restart changed raw native bins");
    }
}
vector saved_spins(std::string const& filename) {
    alps::hdf5::archive archive(filename);
    auto path = "/simulation/realizations/0/clones/0/checkpoint/spins";
    auto shape = archive.extent(path);
    size_t count = 1;
    for (auto dimension : shape) count *= dimension;
    vector spins(count);
    archive.read(path, spins.data(), shape);
    return spins;
}
} // namespace

int main() {
    auto directory = std::filesystem::temp_directory_path() / boost::filesystem::unique_path("spinmc-native-%%%%-%%%%").string();
    std::filesystem::create_directory(directory);
    struct cleanup {
        std::filesystem::path directory;
        ~cleanup() { std::error_code ignored; std::filesystem::remove_all(directory, ignored); }
    } remove{directory};
    auto library = directory / "graphs.xml";
    {
        std::ofstream xml(library);
        xml << R"(<LATTICES>
<GRAPH name="isolated" vertices="1" dimension="1"><VERTEX id="1"><COORDINATE>0</COORDINATE></VERTEX></GRAPH>
<GRAPH name="self" vertices="1" dimension="1"><VERTEX id="1"><COORDINATE>0</COORDINATE></VERTEX><EDGE source="1" target="1" type="0"/></GRAPH>
<GRAPH name="pair" vertices="2" dimension="1"><VERTEX id="1" type="0"><COORDINATE>0</COORDINATE></VERTEX><VERTEX id="2" type="1"><COORDINATE>1</COORDINATE></VERTEX><EDGE source="1" target="2" type="0"/></GRAPH>
<GRAPH name="parallel" vertices="2" dimension="1"><VERTEX id="1" type="0"><COORDINATE>0</COORDINATE></VERTEX><VERTEX id="2" type="1"><COORDINATE>1</COORDINATE></VERTEX><EDGE source="1" target="2" type="0"/><EDGE source="1" target="2" type="1"/></GRAPH>
<GRAPH name="mixed" vertices="3" dimension="1"><VERTEX id="1" type="0"><COORDINATE>0</COORDINATE></VERTEX><VERTEX id="2" type="1"><COORDINATE>1</COORDINATE></VERTEX><VERTEX id="3" type="0"><COORDINATE>2</COORDINATE></VERTEX><EDGE source="1" target="2" type="0"/><EDGE source="2" target="3" type="1"/><EDGE source="1" target="3" type="2"/><EDGE source="3" target="3" type="3"/></GRAPH>
</LATTICES>)";
    }

    // All live continuous matrix representations: scalar, diagonal, symmetric
    // packed and full asymmetric. Full J and scalar D exercise mixed decoding.
    for (size_t dimensions : {2, 3}) {
        auto p = parameters(library, dimensions == 2 ? "XY" : "Heisenberg");
        p["S0"] = .8; p["S1"] = 1.2;
        vector h(dimensions, .17); h[0] = .21; h[1] = -.31; p["h"] = h;
        vector diag(dimensions, .2); diag[0] = .65; diag[1] = -.15;
        vector anisotropy(dimensions, .04); anisotropy[0] = .12; anisotropy[1] = -.08;
        auto d = diagonal(anisotropy); d[1] = d[dimensions] = .07;
        auto symmetric = diagonal(diag); symmetric[1] = symmetric[dimensions] = .19;
        auto asymmetric = symmetric; asymmetric[1] = .35; asymmetric[dimensions] = -.22;
        auto isotropic = diagonal(vector(dimensions, .45));
        for (int shape : {0, 1, 2, 3}) {
            auto j = shape == 0 ? isotropic : shape == 1 ? diagonal(diag) : shape == 2 ? symmetric : asymmetric;
            p["J"] = shape == 0 ? vector{.45} : shape == 1 ? diag : shape == 2 ? packed(j, dimensions) : j;
            auto onsite = shape == 0 ? diagonal(vector(dimensions, .13)) : shape == 1 ? diagonal(anisotropy) : d;
            p["D"] = shape == 0 ? vector{.13} : shape == 1 ? anisotropy : shape == 2 ? packed(d, dimensions) : d;
            physical_moments(p, continuous(dimensions, j, onsite, h, .8, 1.2, 1.1));
        }
        p["GRAPH"] = "parallel"; p["J0"] = asymmetric; p["J1"] = .23; p["D"] = .13;
        auto combined = asymmetric;
        for (size_t i = 0; i != dimensions; ++i) combined[i * dimensions + i] += .23;
        physical_moments(p, continuous(dimensions, combined, diagonal(vector(dimensions, .13)), h, .8, 1.2, 1.1));
        p["GRAPH"] = "self"; p["J0"] = asymmetric;
        physical_moments(p, continuous(dimensions, asymmetric, diagonal(vector(dimensions, .13)), h, .8, 0, 1.1, true));
    }

    // Enumeration independently checks mixed signs/types, quantum normalization,
    // onsite constants, fields, self-loops and derived thermodynamic identities.
    for (bool quantum : {false, true}) {
        auto p = parameters(library, "Ising", "mixed");
        p["S0"] = .7; p["S1"] = 1.3; p["g"] = 1.2; p["h"] = vector{-.35};
        p["D0"] = .17; p["D1"] = -.09;
        p["J0"] = .7; p["J1"] = -.4; p["J2"] = .3; p["J3"] = 1.1;
        if (quantum) p["CONVENTION"] = "quantum";
        // Quantum-field combinations are intentionally rejected in the public contract.
        if (quantum) p["h"] = vector{0.};
        vector length{quantum ? std::sqrt(.7 * 1.7) : .7, quantum ? std::sqrt(1.3 * 2.3) : 1.3,
                      quantum ? std::sqrt(.7 * 1.7) : .7};
        double z = 0, e = 0, e2 = 0, abs_m = 0, m2 = 0, m4 = 0, em2 = 0, em4 = 0;
        for (int state = 0; state != 8; ++state) {
            double a = state & 1 ? length[0] : -length[0], b = state & 2 ? length[1] : -length[1],
                   c = state & 4 ? length[2] : -length[2];
            double energy = (quantum ? 1 : -1) * (.7*a*b - .4*b*c + .3*a*c + 1.1*c*c)
                          - .17*.7*.7*2 + .09*1.3*1.3 + (quantum ? 0 : 1.2*.35*(a+b+c));
            double m = (a+b+c)/3, weight = std::exp(-energy/1.1);
            z += weight; e += weight*energy; e2 += weight*energy*energy; abs_m += weight*std::abs(m);
            m2 += weight*m*m; m4 += weight*std::pow(m,4); em2 += weight*energy*m*m; em4 += weight*energy*std::pow(m,4);
        }
        auto measured = spinmc::derive(sample(p), p);
        agree(measured.at("Energy"), e/z, "Ising Hamiltonian disagrees with enumeration");
        agree(measured.at("Magnetization^2"), m2/z, "Ising spin factors are incorrect");
        agree(measured.at("Specific Heat"), (e2/z - std::pow(e/z,2))/(3*1.1*1.1), "heat capacity is incorrect");
        agree(measured.at("Binder Cumulant"), m4*z/(m2*m2), "Binder U4 is incorrect");
        agree(measured.at("Binder Cumulant U2"), m2*z/(abs_m*abs_m), "Binder U2 is incorrect");
        agree(measured.at("Magnetization^2 slope"), (em2/z - e*m2/(z*z))/(1.1*1.1), "temperature derivative is incorrect");
        agree(measured.at("Magnetization^4 slope"), (em4/z - e*m4/(z*z))/(1.1*1.1), "fourth-moment derivative is incorrect");
    }

    for (int q : {3, 4, 10})
        for (double j : {.7, -.7}) {
            auto p = parameters(library, "Potts", "pair", j > 0 ? "cluster" : "auto");
            p["q"] = q; p["J"] = j;
            double z = 0, energy = 0;
            for (int a = 0; a != q; ++a) for (int b = 0; b != q; ++b) {
                double e = a == b ? -j : 0, weight = std::exp(-e/1.1);
                z += weight; energy += e*weight;
            }
            spinmc::simulation run(p, 64);
            require(run.effective_update() == (j > 0 ? "cluster" : "local"), "Potts auto update is invalid");
            finish(run); agree(run.collect_results().at("Energy"), energy/z, "Potts Boltzmann distribution is incorrect");
            if (j < 0) {
                p["UPDATE"] = "cluster";
                rejects([&] { spinmc::simulation invalid(p); }, "negative Potts couplings accepted cluster updates");
            }
        }

    // Isolated spins in fields expose ignored-field default dispatch. O(4) is
    // checked against its sphere integral rather than another sampler.
    for (auto const& model : {"Ising", "XY", "Heisenberg", "O(4)", "Potts"}) {
        auto p = parameters(library, model, "isolated", "auto");
        p["h"] = vector{-.9}; p["D"] = 0.; p["q"] = 3;
        spinmc::simulation run(p, 64);
        require(run.effective_update() == "local", "auto cluster ignored a nonzero field");
        finish(run);
        double x = .9/1.1, magnetization = std::string(model) == "Ising" ? std::tanh(x)
            : std::string(model) == "XY" ? boost::math::cyl_bessel_i(1,x)/boost::math::cyl_bessel_i(0,x)
            : std::string(model) == "Heisenberg" ? 1/std::tanh(x)-1/x
            : std::string(model) == "O(4)" ? boost::math::cyl_bessel_i(2,x)/boost::math::cyl_bessel_i(1,x)
            : -std::exp(-x)/(std::exp(-x)+2);
        agree(run.collect_results().at("Energy"), -.9*magnetization, "isolated-spin field energy is incorrect");
        if (std::string(model) != "Potts")
            agree(run.collect_results().at("Magnetization along Field"), magnetization, "negative-field projection is incorrect");
        p["UPDATE"] = "cluster";
        rejects([&] { spinmc::simulation invalid(p); }, "explicit clusters accepted a field");
    }

    for (auto const& model : {"Ising", "XY", "Heisenberg", "O(4)"})
        for (double j : {.7, -.7}) {
            auto p = parameters(library, model, "pair", "cluster"); p["J"] = j;
            auto measured = sample(p);
            double x = .7/1.1;
            size_t dimensions = std::string(model) == "Ising" ? 1 : std::string(model) == "XY" ? 2 : std::string(model) == "Heisenberg" ? 3 : 4;
            double correlation = dimensions == 1 ? std::tanh(x) : dimensions == 2 ? boost::math::cyl_bessel_i(1,x)/boost::math::cyl_bessel_i(0,x)
                : dimensions == 3 ? 1/std::tanh(x)-1/x : boost::math::cyl_bessel_i(2,x)/boost::math::cyl_bessel_i(1,x);
            agree(measured.at("Energy"), -.7*correlation, "cluster distribution is incorrect");
            agree(measured.at("Magnetization^2"), (1+(j > 0 ? correlation : -correlation))/2, "cluster magnetization is incorrect");
            auto prefix = j > 0 ? "Improved " : "Improved Staggered ";
            agree(measured.at(std::string(prefix)+"Magnetization^2"), (1+correlation)/2, "cluster improved estimator is incorrect");
            agree(measured.at(std::string(prefix)+"Susceptibility"), (1+correlation)/(1.1*dimensions), "cluster susceptibility normalization is incorrect");
        }

    // Full continuation evidence includes the actual spins and RNG, not just
    // means. Stop during warmup, at its boundary and inside an ALEA merge.
    for (auto const& model : {"Ising", "XY", "Heisenberg", "O(4)", "Potts"})
        for (auto const& update : {"local", "cluster"})
            for (uint64_t split : {3, 17, 148}) {
                auto p = parameters(library, model, "pair", update);
                p["q"] = 4; p["THERMALIZATION"] = 17; p["SWEEPS"] = 237;
                spinmc::simulation reference(p, 16, 3), interrupted(p, 16, 3), resumed(p, 16, 3);
                finish(reference); step(interrupted, split);
                if (std::string(update) == "local" && split <= 17)
                    require(interrupted.measurement_count() == 0, "warmup leaked into production measurements");
                auto checkpoint = (directory / "state.h5").string();
                interrupted.save(checkpoint); resumed.load(checkpoint); finish(resumed);
                equal(reference, resumed);
                reference.save((directory / "reference.h5").string()); resumed.save((directory / "resumed.h5").string());
                require(saved_spins((directory / "reference.h5").string()) == saved_spins((directory / "resumed.h5").string()), "restart changed spin state");
                require(reference.get_random()() == resumed.get_random()(), "restart changed RNG continuation");
                auto counts = resumed.collect_results();
                {
                    alps::hdf5::archive archive(checkpoint, "a");
                    archive["/simulation/realizations/0/clones/0/checkpoint/chain_id"] << uint64_t(4);
                }
                rejects([&] { resumed.load(checkpoint); }, "checkpoint accepted a different chain");
                for (auto const& entry : counts)
                    require(entry.second.store().batch() == resumed.collect_results().at(entry.first).store().batch(), "failed restart changed measurements");
                interrupted.save(checkpoint);
                auto invalid = saved_spins(checkpoint);
                std::fill(invalid.begin(), invalid.end(), 99.);
                {
                    alps::hdf5::archive archive(checkpoint, "a");
                    auto path = "/simulation/realizations/0/clones/0/checkpoint/spins";
                    auto shape = archive.extent(path);
                    archive.write(path, invalid.data(), shape);
                }
                rejects([&] { resumed.load(checkpoint); }, "checkpoint accepted an invalid spin state");
                for (auto const& entry : counts)
                    require(entry.second.store().batch() == resumed.collect_results().at(entry.first).store().batch(), "invalid spin load changed measurements");
                equal(reference, resumed);
                require(reference.get_random()() == resumed.get_random()(), "failed restart changed RNG continuation");
            }
    auto zero = parameters(library, "Ising", "isolated"); zero["THERMALIZATION"] = 0; zero["SWEEPS"] = 1;
    spinmc::simulation first(zero, 8); step(first, 1);
    require(first.measurement_count() == 1 && first.fraction_completed() == 1, "zero warmup lost the first sample");
    // Always-flip Metropolis observations after N attempts trap even-N systems
    // in one parity sector when every move accepts. Infinite-temperature local
    // updates must still sample all four pair states, not just ++ and --.
    auto hot = parameters(library, "Ising"); hot.erase("T"); hot["beta"] = 0.; hot["J"] = 0.;
    auto hot_results = sample(hot);
    agree(hot_results.at("Magnetization^2"), .5, "infinite-temperature Ising observations preserve a spurious parity sector");
    agree(hot_results.at("|Magnetization|"), .5, "infinite-temperature Ising states are not equiprobable");

    // Raw moments cannot resolve a tiny variance around a nonzero energy.
    // Keep the evidence, but do not amplify rounding into a heat capacity.
    auto cold = parameters(library, "Ising", "pair", "cluster");
    cold.erase("T"); cold["beta"] = 1e100; cold["THERMALIZATION"] = 0; cold["SWEEPS"] = 64;
    auto frozen = sample(cold);
    auto frozen_derived = spinmc::derive(frozen, cold);
    require(!frozen_derived.count("Specific Heat"), "frozen energy published a cancellation-dominated heat capacity");
    require(frozen_derived.at("Energy").store().batch() == frozen.at("Energy").store().batch(),
            "conditioning rejection changed raw energy evidence");
    auto energy = frozen.at("Energy").store(), energy2 = frozen.at("Energy^2").store();
    for (Eigen::Index i = 0; i < energy.count().size(); ++i) {
        double e = .7 + (i % 2 ? 1e-8 : -1e-8);
        energy.batch()(0,i) = double(energy.count()(i)) * e;
        energy2.batch()(0,i) = double(energy2.count()(i)) * e * e;
    }
    frozen["Energy"] = result(energy); frozen["Energy^2"] = result(energy2);
    frozen["E.Magnetization^2"] = result(energy);
    frozen["E.Magnetization^4"] = result(energy);
    auto unresolved = spinmc::derive(frozen, cold);
    require(!unresolved.count("Specific Heat"), "underresolved real energy variance was clamped or published");
    require(unresolved.at("Energy").store().batch() == energy.batch()
            && unresolved.at("Energy^2").store().batch() == energy2.batch(),
            "conditioning rejection changed raw moments");
    auto hot_derived = spinmc::derive(hot_results, hot);
    require(hot_derived.at("Specific Heat").mean()(0) == 0., "beta zero heat capacity is not exactly zero");
    auto potts_field = parameters(library, "Potts", "isolated");
    potts_field["q"] = 3; potts_field["h"] = vector{.1, 0.};
    rejects([&] { spinmc::simulation invalid(potts_field); }, "Potts silently discarded a field component");

    auto warmup = parameters(library, "Ising", "pair", "cluster"); warmup["THERMALIZATION"] = 170;
    spinmc::simulation fresh(warmup, 8);
    auto checkpoint = (directory / "warmup.h5").string();
    for (uint64_t updates : {0, 1}) {
        fresh.save(checkpoint);
        auto rng = fresh.get_random();
        {
            alps::hdf5::archive archive(checkpoint, "a");
            archive["/simulation/realizations/0/clones/0/checkpoint/updates"] << updates;
            archive["/simulation/realizations/0/clones/0/checkpoint/warmup_updates"] << updates;
            archive["/simulation/realizations/0/clones/0/checkpoint/warmup_sites"] << uint64_t(updates ? 3 : 100);
        }
        rejects([&] { fresh.load(checkpoint); }, "checkpoint accepted impossible cluster warmup counters");
        require(fresh.completed_sweeps() == 0 && fresh.measurement_count() == 0
                && fresh.get_random() == rng, "failed warmup load changed continuation state");
    }
    std::cout << "spinmc independent Hamiltonian, cluster and exact restart gates passed\n";
}
