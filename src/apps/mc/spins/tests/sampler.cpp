// SPDX-License-Identifier: MIT
// Drive the spin simulation without the scheduler, with its own parameters,
// measurements and random number generator, as another driver would. From
// the same inputs it must measure exactly what the scheduler run measures.
#include "factory.h"
#include "spinsim.h"

#include <alps/osiris/comm.h>
#include <alps/random/rngfactory.h>
#include <alps/scheduler/montecarlo.h>
#include <boost/shared_ptr.hpp>
#include <gtest/gtest.h>

#include <memory>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

namespace {

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

// Count, mean and error of every observable, exactly.
std::vector<std::string> summary(const alps::ObservableSet& measurements) {
    std::vector<std::string> out;
    for (const auto& entry : measurements) {
        std::ostringstream line;
        line << std::hexfloat << entry.first;
        if (const auto* scalar = dynamic_cast<const alps::RealObservable*>(entry.second)) {
            line << " " << scalar->count();
            if (scalar->count())
                line << " " << scalar->mean() << " " << scalar->error();
        } else if (const auto* vector =
                       dynamic_cast<const alps::RealVectorObservable*>(entry.second)) {
            line << " " << vector->count();
            if (vector->count())
                for (double value : std::valarray<double>(vector->mean()))
                    line << " " << value;
        }
        out.push_back(line.str());
    }
    return out;
}

// Seeds as a scheduler worker does, before the simulation is built: the
// engine from SEED, the disorder generator from DISORDERSEED.
template <class Simulation>
std::vector<std::string> sampled(const std::string& text, int sweeps) {
    alps::Parameters p = parameters(text);
    boost::shared_ptr<alps::buffered_rng_base> engine(alps::rng_factory.create("mt19937"));
    engine->seed(static_cast<uint32_t>(p["SEED"]));
    alps::Disorder::seed(p.value_or_default("DISORDERSEED", 0));
    typename Simulation::random_type random(*engine, boost::uniform_real<>());
    alps::ObservableSet measurements;

    Simulation simulation(p, measurements, random);
    for (int i = 0; i < sweeps; ++i)
        simulation.dostep();
    EXPECT_TRUE(simulation.is_thermalized());
    return summary(measurements);
}

std::vector<std::string> scheduled(const std::string& text, int sweeps) {
    const alps::ProcessList here(1, alps::local_process());
    std::unique_ptr<alps::scheduler::Worker> worker(
        SpinFactory().make_worker(here, parameters(text), 0));
    for (int i = 0; i < sweeps; ++i)
        worker->dostep();
    return summary(dynamic_cast<const alps::scheduler::MCRun&>(*worker).get_measurements());
}

const std::string ising_chain =
    R"p(LATTICE="chain lattice";L=8;MODEL="Ising";J=1;T=1.5;SWEEPS=95;UPDATE="local";SEED=11)p";
const std::string ising_square =
    R"p(LATTICE="square lattice";L=4;MODEL="Ising";J=1;T=2.5;SWEEPS=50;THERMALIZATION=5;SEED=3)p";
const std::string two_band =
    R"p(LATTICE="2 band chain lattice";L=4;MODEL="Heisenberg";J0="random()";J1="0.5+random()";J2=-1;S0="0.5+random()";S1="1+random()";D0="random() 0.2 0.3";D1="random()";DISORDERSEED=7;T=1;SWEEPS=40;THERMALIZATION=4;SEED=13)p";
const std::string xy =
    R"p(LATTICE="chain lattice";L=6;MODEL="XY";J="1 0.5";D="0.1 0.2";CONVENTION="quantum";S=1.5;T=0.7;SWEEPS=50;THERMALIZATION=5;SEED=6)p";

TEST(SpinmcSampler, MeasuresLikeTheSchedulerRun) {
    EXPECT_EQ((sampled<SpinSim<IsingMoment, MIdMatrix<double, 1>>>(ising_chain, 30)),
              scheduled(ising_chain, 30));
    EXPECT_EQ((sampled<SpinSim<IsingMoment, MIdMatrix<double, 1>>>(ising_square, 25)),
              scheduled(ising_square, 25));
    EXPECT_EQ((sampled<SpinSim<ONMoment<3>, DiagMatrix<double, 3>>>(two_band, 20)),
              scheduled(two_band, 20));
    EXPECT_EQ((sampled<SpinSim<XYMoment, DiagMatrix<double, 2>>>(xy, 20)), scheduled(xy, 20));
}

// Current behavior, kept from the scheduler run: building the simulation
// inserts an empty UPDATE into the parameters it is given.
TEST(SpinmcSampler, InsertsUpdateIntoTheGivenParameters) {
    alps::Parameters p = parameters(ising_square);
    boost::shared_ptr<alps::buffered_rng_base> engine(alps::rng_factory.create("mt19937"));
    SpinSim<IsingMoment, MIdMatrix<double, 1>>::random_type random(*engine,
                                                                   boost::uniform_real<>());
    alps::ObservableSet measurements;
    SpinSim<IsingMoment, MIdMatrix<double, 1>> simulation(p, measurements, random);
    ASSERT_TRUE(p.defined("UPDATE"));
    EXPECT_EQ(std::string(p["UPDATE"]), "");
    EXPECT_TRUE(measurements.has("Cluster size"));
}

} // namespace
