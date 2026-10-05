/* ndim_spin.hpp
 * adapted from heisenberg_lattice/heisenberg.hpp
 */

#ifndef HEISENBERG_HPP
#define HEISENBERG_HPP

#include "tinyvector/tinyvector.hpp"

#include <alps/mcbase.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/random/uniform_on_sphere_n.h>
#include <alps/lattice.h>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/array.hpp>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/array.hpp>

#include <vector>
#include <string>

/**
 * The simulation class derives from mcbase.
 */
template <int N>
class ALPS_DECL ndim_spin_sim : public alps::mcbase {

    public:
    using results_type = std::map<std::string,alps::alea::batch_result<double>>;
    results_type collect_results(result_names_type const& names={}) const {
        return collect_results_as<alps::alea::batch_result<double>>(names);
    }

        typedef tinyvector<double, N, INTRIN_OPT> spintype;
        ndim_spin_sim(parameters_type const & parms, std::size_t seed_offset = 0);

        void update() override;
        void measure() override;
        double fraction_completed() const override;

        using alps::mcbase::save;
        void save(alps::hdf5::archive & ar) const override;

        using alps::mcbase::load;
        void load(alps::hdf5::archive & ar) override;

        // convenience function to get a random spin (uniformly distancesributed direction)
        const spintype random_spin();

    private:
        
        alps::graph_helper<> lattice;
        int num_sites;
        std::vector<spintype> spins;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        alps::uniform_on_sphere_n<N, double, spintype > random_spin_gen;
        std::vector<double> distances;
};

/* implementation */

template<int N>
ndim_spin_sim<N>::ndim_spin_sim(parameters_type const & parms, std::size_t seed_offset)
    : alps::mcbase(parms, seed_offset)
    , lattice(alps::make_deprecated_parameters(parms))
    , num_sites(lattice.num_sites())
    , spins(num_sites)
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , random_spin_gen()
{
    // initialize spins with random values
    for(int i = 0; i < num_sites; ++i) {
        spins[i] = random_spin();
    }

    measurements.emplace("Energy", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
    measurements.emplace("Magnetization", std::make_shared<alps::alea::batch_acc<double>>(N, 64));
    measurements.emplace("Magnetization^2", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
    measurements.emplace("Magnetization^4", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
    measurements.emplace("Correlations", std::make_shared<alps::alea::batch_acc<double>>(num_sites, 64));
    measurements.emplace("Distances", std::make_shared<alps::alea::batch_acc<double>>(num_sites, 64));

    std::vector<double> ref = lattice.coordinate(0);
    std::vector<double> a;
    for (int i = 0; i < num_sites; ++i) {
        double d = 0;
        a = lattice.coordinate(i);
        for (std::size_t j = 0; j < a.size(); ++j) a[j] -= ref[j];
        for (int j = 0; j < a.size(); ++j) {
            d += a[j] * a[j];
        }
        distances.push_back(std::sqrt(d));
    }
    *measurement("Distances") << alps::alea::make_adapter(distances);
}

template<int N>
void ndim_spin_sim<N>::update() {
    for (int j = 0; j < num_sites; ++j) {
        using std::exp;
        int i = int(double(num_sites) * random());

        // get a random site
        alps::graph_helper<>::site_descriptor site_i = lattice.site(i);

        // sum over all it's neighbors for delta_H later
        typename ndim_spin_sim<N>::spintype nn_sum;
        nn_sum.initialize(0);
        alps::graph_helper<>::neighbor_iterator nn_it, nn_end;
        for(boost::tie(nn_it, nn_end) = lattice.neighbors(site_i); nn_it != nn_end; ++nn_it) {
            nn_sum += spins[*nn_it];
        }

        // generate a new random spin and decide if we keep it
        typename ndim_spin_sim<N>::spintype new_spin = random_spin();
        double delta_H = dot(new_spin - spins[i], nn_sum);
        double p = exp(beta * delta_H);
        if ( p >= 1. || random() < p )
            spins[i] = new_spin;
    }
}

template <int N>
void ndim_spin_sim<N>::measure() {
    sweeps++;
    if (sweeps > thermalization_sweeps) {
        typename ndim_spin_sim<N>::spintype magnetization;
        magnetization.initialize(0);
        double energy = 0;
        std::vector<double> correlations(num_sites, 0);
        // To measure magnetization, magnetic susceptibility and correlationselations we sum over all sites.
        for (int i = 0; i < num_sites; ++i) {
            magnetization += spins[i];
            correlations[i] = dot(spins[0], spins[i]);
        }
        // To measure the Energy we sum only over neighbored sites (bonds in terms of the lattice).
        alps::graph_helper<>::bond_iterator bond_it, bond_end;
        for(boost::tie(bond_it, bond_end) = lattice.bonds(); bond_it != bond_end; ++bond_it) {
            energy += - dot(spins[lattice.source(*bond_it)], spins[lattice.target(*bond_it)]);
        }
        energy /= num_sites;                // $\frac{1}{V} \sum_{\text{i,j nn}}{\sigma_i \sigma_j}$
        magnetization /= num_sites;         // $\frac{1}{V} \sum_{i}{\sigma_i}$
        double magnetization2 = dot(magnetization, magnetization);

        // store the measurements
        *measurement("Energy") << alps::alea::make_adapter(energy);
        *measurement("Magnetization") << alps::alea::make_adapter(spintype::vector(magnetization));
        *measurement("Magnetization^2") << alps::alea::make_adapter(magnetization2);
        *measurement("Magnetization^4") << alps::alea::make_adapter(magnetization2 * magnetization2);
        *measurement("Correlations") << alps::alea::make_adapter(correlations);
    }
}

template <int N>
double ndim_spin_sim<N>::fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

template <int N>
void ndim_spin_sim<N>::save(alps::hdf5::archive & ar) const {
    mcbase::save(ar);
    ar["checkpoint/sweeps"] << sweeps;
    ar["checkpoint/spins"] << spins;
}

template <int N>
void ndim_spin_sim<N>::load(alps::hdf5::archive & ar) {
    parameters_type restored_parameters;
    ar["/parameters"] >> restored_parameters;
    if (restored_parameters["L"].template as<int>() != parameters["L"].template as<int>())
        throw std::invalid_argument("checkpoint geometry differs from the simulation");
    for (auto const* key : {"LATTICE", "LATTICE_LIBRARY"})
        if (restored_parameters.template value_or<std::string>(key, "") != parameters.template value_or<std::string>(key, ""))
            throw std::invalid_argument("checkpoint lattice differs from the simulation");
    int restored_thermalization = restored_parameters["THERMALIZATION"].template as<int>();
    int restored_total = restored_parameters["SWEEPS"].template as<int>();
    double restored_beta = 1. / restored_parameters["T"].template as<double>();
    int restored_sweeps;
    std::vector<spintype> restored_spins;
    ar["checkpoint/sweeps"] >> restored_sweeps;
    ar["checkpoint/spins"] >> restored_spins;
    if (!std::isfinite(restored_beta) || restored_beta <= 0 || restored_thermalization < 0
            || restored_total <= 0 || restored_sweeps < 0
            || std::int64_t(restored_sweeps) > std::int64_t(restored_thermalization) + restored_total
            || restored_spins.size() != std::size_t(num_sites))
        throw std::invalid_argument("invalid Heisenberg checkpoint progress or shape");
    for (auto const& spin : restored_spins)
        if (!std::isfinite(dot(spin, spin)) || std::abs(dot(spin, spin) - 1.) > 1e-12)
            throw std::invalid_argument("invalid Heisenberg checkpoint spin");
    mcbase::load(ar);
    thermalization_sweeps = restored_thermalization;
    total_sweeps = restored_total;
    beta = restored_beta;
    sweeps = restored_sweeps;
    spins = std::move(restored_spins);
}

template <int N>
const typename ndim_spin_sim<N>::spintype ndim_spin_sim<N>::random_spin() {
    return random_spin_gen(random.engine());
}

#endif 
