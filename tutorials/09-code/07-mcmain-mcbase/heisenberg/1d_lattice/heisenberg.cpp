/* heisenberg.cpp
 * adapted from alps/tutorials/09-code/07-mcmain-mcbase/heisenberg.cpp
 * 1d grid, 3d spins
 */

#include "heisenberg.hpp"
#include "helper.hpp"


heisenberg_sim::heisenberg_sim(parameters_type const & parms, std::size_t seed_offset)
    : alps::mcbase(parms, seed_offset)
    , length(parameters["L"])
    , sweeps(0)
    , thermalization_sweeps(int(parameters["THERMALIZATION"]))
    , total_sweeps(int(parameters["SWEEPS"]))
    , beta(1. / double(parameters["T"]))
    , random_spin_gen()
    , spins(length)
{
    /* initialize spins */
    for(int i = 0; i < length; ++i)
    {
        spins[i] = random_spin();
    }
    measurements.emplace("Energy", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
    measurements.emplace("Magnetization", std::make_shared<alps::alea::batch_acc<double>>(3, 64));
    measurements.emplace("Magnetization^2", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
    measurements.emplace("Magnetization^4", std::make_shared<alps::alea::batch_acc<double>>(1, 64));
    measurements.emplace("Correlations", std::make_shared<alps::alea::batch_acc<double>>(length, 64));
}

void heisenberg_sim::update() {
    for (int j = 0; j < length; ++j) {
        using std::exp;
        int i = int(double(length) * random());
        int right = ( i + 1 < length ? i + 1 : 0 );
        int left = ( i - 1 < 0 ? length - 1 : i - 1 );
        /* generate a new random spin for update and accept or not */
        spintype new_spin = random_spin();
        double delta_H = dot(new_spin - spins[i], spins[right] + spins[left]);
        double p = exp(beta * delta_H);
        if ( p >= 1. || random() < p )
            spins[i] = new_spin;
    }
}

void heisenberg_sim::measure() {
    sweeps++;
    if (sweeps > thermalization_sweeps) {
        spintype tmag = {{0, 0, 0}};
        double ten = 0;
        //~ double sign = 1;
        std::vector<double> corr(length);
        for (int i = 0; i < length; ++i) {
            tmag += spins[i];
            //~ sign *= spins[i];
            ten += dot(-spins[i], spins[ i + 1 < length ? i + 1 : 0 ]);
            for (int d = 0; d < length; ++d)
                corr[d] += dot(spins[i], spins[( i + d ) % length ]);
        }
        for (auto& correlation : corr) correlation /= length;
        ten /= length;
        tmag /= length;
        *measurement("Energy") << alps::alea::make_adapter(ten);
        *measurement("Magnetization") << alps::alea::make_adapter(vector_from_spintype(tmag));
        *measurement("Magnetization^2") << alps::alea::make_adapter(dot(tmag, tmag));
        *measurement("Magnetization^4") << alps::alea::make_adapter(dot(tmag, tmag) * dot(tmag, tmag));
        *measurement("Correlations") << alps::alea::make_adapter(corr);
    }
}

double heisenberg_sim::fraction_completed() const {
    return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
}

void heisenberg_sim::save(alps::hdf5::archive & ar) const {
    mcbase::save(ar);
    ar["checkpoint/sweeps"] << sweeps;
    ar["checkpoint/spins"] << spins;
}

void heisenberg_sim::load(alps::hdf5::archive & ar) {
    parameters_type restored_parameters;
    ar["/parameters"] >> restored_parameters;
    if (restored_parameters["L"].as<int>() != parameters["L"].as<int>())
        throw std::invalid_argument("checkpoint geometry differs from the simulation");
    int restored_thermalization = restored_parameters["THERMALIZATION"].as<int>();
    int restored_total = restored_parameters["SWEEPS"].as<int>();
    double restored_beta = 1. / restored_parameters["T"].as<double>();
    int restored_sweeps;
    std::vector<spintype> restored_spins;
    ar["checkpoint/sweeps"] >> restored_sweeps;
    ar["checkpoint/spins"] >> restored_spins;
    if (!std::isfinite(restored_beta) || restored_beta <= 0 || restored_thermalization < 0
            || restored_total <= 0 || restored_sweeps < 0
            || std::int64_t(restored_sweeps) > std::int64_t(restored_thermalization) + restored_total
            || restored_spins.size() != std::size_t(length))
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

const spintype heisenberg_sim::random_spin() {
    return spin_from_vector(random_spin_gen(random.engine()));
}

const spintype heisenberg_sim::random_spin_accept_reject() {
    double r = 0;
    spintype spin = { {0, 0, 0} };
    /* get a random vector from uniformly distributed numbers in [-1,1) */
    do 
    {
        for(int i = 0; i < 3; ++i)
        {
            spin[i] = 2. * random() - 1.;
        }
        r = abs(spin);
    }
    while(r >= 1 || r == 0); /* reject and replace if not in unit sphere */
    /* normalize to get a point on the 3d unit sphere with uniformly distributed angles*/
    for(int i = 0; i < 3; ++i)
    {
        spin[i] /= r;
    } 
    return spin;
}
