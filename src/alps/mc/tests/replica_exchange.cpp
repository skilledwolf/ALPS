// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/mc/replica_exchange.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <tuple>

namespace {
void require(bool condition,char const* message) {
    if (!condition) throw std::runtime_error(message);
}
alps::params parameters() {
    alps::params p;
    p["INVERSE_TEMPERATURE_SET"]=std::vector<double>{.25,.7,1.3};
    p["SWEEPS"]=200000; p["THERMALIZATION"]=1000; p["SEED"]=89;
    return p;
}
double log_weight(double energy,double beta) {return -beta*energy;}
using exchange=alps::mc::replica_exchange<double>;
using event=std::tuple<size_t,std::string,double>;
std::vector<event> scripted_step(exchange& state) {
    std::vector<event> events;
    state.step([&](size_t walker,size_t slot,double beta,bool sampling) {
        events.emplace_back(slot,"walker",double(walker));
        events.emplace_back(slot,"beta",beta);
        events.emplace_back(slot,"sampling",double(sampling));
    },[]{return std::vector<double>{-4.,-2.,1.};},log_weight,
    [&](size_t slot,char const* name,double value){events.emplace_back(slot,name,value);});
    return events;
}

// Two-state physical walkers with E=+-2 have exactly known Gibbs moments.
// Their state persists through both local updates and temperature exchanges;
// an incorrect exchange sign consequently changes the measured distribution.
void thermodynamics(bool randomized,bool disabled,bool zero_beta=false) {
    auto p=parameters(); p["RANDOM_EXCHANGE"]=randomized; p["NO_EXCHANGE"]=disabled;
    if (zero_beta) p["INVERSE_TEMPERATURE_SET"]=std::vector<double>{0.,.7,1.3};
    exchange state(p,0,0.,zero_beta);
    alps::random01 random(773);
    std::vector<double> energies{-2.,2.,-2.},sum(3),square(3);
    std::vector<uint64_t> counts(3);
    while (state.fraction_completed()<1) {
        state.step([&](size_t walker,size_t slot,double beta,bool sampling) {
            const double change=-2*energies[walker];
            // A full two-state heat bath would redraw the exact distribution
            // every sweep and hide exchange errors. Retain temporal dependence.
            if (random()<.2 && random()<1/(1+std::exp(beta*change))) energies[walker]=-energies[walker];
            if (sampling) {sum[slot]+=energies[walker];square[slot]+=energies[walker]*energies[walker];++counts[slot];}
        },[&]{return energies;},log_weight,[](size_t,char const*,double){});
    }
    for (size_t i=0;i<3;++i) {
        require(counts[i]==200000,"Exchange lost or duplicated production samples");
        require(std::abs(sum[i]/counts[i]+2*std::tanh(2*state.beta(i)))<.025,"Scalar exchange violates Gibbs thermodynamics");
        require(square[i]/counts[i]==4,"Exchange changed physical energy bounds");
    }
    require(state.completed_sweeps()==201000,"Exchange warmup count changed");
}

void identity_diagnostics() {
    auto p=parameters(); p["THERMALIZATION"]=0;
    exchange state(p,0,0.);
    // Equal weights accept every edge. After six sweeps walker 0 returns to
    // the hottest slot; after eight sweeps walker 2 returns there instead.
    for (size_t step=1;step<=8;++step) {
        size_t returns=0;
        state.step([](size_t,size_t,double,bool){},[]{return std::vector<double>(3,0.);},log_weight,
        [&](size_t index,char const* name,double value) {
            const std::string label(name);
            if (label=="EXMC: Acceptance Rate")
                require(index==step%2 && value==1,"Acceptance diagnostic no longer refers to its temperature edge");
            if (label=="EXMC: Inverse Round-Trip Time" && value) {
                require((step==6 && index==0) || (step==8 && index==2),"Round-trip diagnostic lost physical walker identity");
                ++returns;
            }
        });
        require(returns==size_t(step==6 || step==8),"Round-trip completion was lost or duplicated");
    }
}

void classical_zero_endpoint() {
    auto p=parameters();p["INVERSE_TEMPERATURE_SET"]=std::vector<double>{0.,.4,1.3};
    bool rejected=false;
    try {exchange quantum(p,0,0.);} catch (std::invalid_argument const&) {rejected=true;}
    require(rejected,"Quantum grids accepted a zero inverse temperature");
    alps::mc::temperature_grid grid(p,true);
    grid.optimize_rate(std::vector<double>{0.,-2.,-3.},log_weight);
    require(grid[0]==0 && grid[2]==1.3 && grid[1]>0 && grid[1]<1.3 && std::abs(grid[1]-.4)>.001,
            "Classical zero endpoint disabled rate feedback");
    p["OPTIMIZE_TEMPERATURE"]=true;p["OPTIMIZATION_TYPE"]="population";
    rejected=false;
    try {exchange population(p,0,0.,true);} catch (std::invalid_argument const&) {rejected=true;}
    require(rejected,"Population feedback accepted an infinite temperature endpoint");
}

void checkpoints(boost::filesystem::path const& path,std::string const& method,std::string const& rng) {
    auto p=parameters(); p["OPTIMIZE_TEMPERATURE"]=true; p["OPTIMIZATION_TYPE"]=method;
    p["INITIAL_BLOCK_SWEEPS"]=5; p["OPTIMIZATION_ITERATIONS"]=1;
    p["THERMALIZATION"]=3; p["RNG"]=rng;
    exchange state(p,0,0.);
    for (int step=0;step<9;++step) scripted_step(state);
    {alps::hdf5::archive ar(path,"w");state.save(ar);}
    exchange restored(p,0,0.);
    {
        alps::hdf5::archive ar(path);
        std::vector<std::vector<double>> weights; ar["exchange/weights"]>>weights;
        require(weights.size()==3 && weights[0].size()==1,"Scalar exchange checkpoint shape changed");
        require(std::any_of(weights.begin(),weights.end(),[](auto const& w){return w[0]<0;}),"Signed feedback fixture lacks negative weights");
        restored.load(ar);
    }
    for (int step=0;step<200;++step)
        require(scripted_step(state)==scripted_step(restored),"Exchange restart changed feedback, routing or RNG");

    // A rejected load must leave the complete live state, including RNG,
    // unchanged. Duplicate walker IDs cannot form a temperature permutation.
    auto unchanged=restored;
    {alps::hdf5::archive ar(path,"a");ar["exchange/walkers"]<<std::vector<size_t>{0,0,2};}
    bool failed=false;
    try {alps::hdf5::archive ar(path);restored.load(ar);} catch (std::invalid_argument const&) {failed=true;}
    require(failed,"Invalid exchange permutation was accepted");
    for (int step=0;step<20;++step)
        require(scripted_step(unchanged)==scripted_step(restored),"Rejected exchange checkpoint mutated live state");
}

void rewind_feedback(boost::filesystem::path const& path) {
    auto p=parameters();p["OPTIMIZE_TEMPERATURE"]=true;p["OPTIMIZATION_TYPE"]="population";
    p["INITIAL_BLOCK_SWEEPS"]=1;p["OPTIMIZATION_ITERATIONS"]=1;
    exchange state(p,0,0.),fresh(p,0,0.);
    {alps::hdf5::archive ar(path,"w");state.save(ar);}
    scripted_step(state); // Too few round trips: the live optimization blocks grow.
    {alps::hdf5::archive ar(path);state.load(ar);}
    for (int step=0;step<20;++step)
        require(scripted_step(state)==scripted_step(fresh),"Rewinding feedback changed the initial schedule or RNG");
}
}

int main() {
    const auto path=boost::filesystem::temp_directory_path()/boost::filesystem::unique_path("alps-exchange-%%%%-%%%%.h5");
    try {
        thermodynamics(false,false);thermodynamics(true,false);thermodynamics(false,true);
        thermodynamics(true,false,true);
        identity_diagnostics();
        classical_zero_endpoint();
        for (auto method:{"rate","population"}) for (auto rng:{"mt19937","lagged_fibonacci607"}) checkpoints(path,method,rng);
        rewind_feedback(path);
        boost::filesystem::remove(path);
        return 0;
    } catch (std::exception const& error) {
        boost::filesystem::remove(path);
        std::cerr<<error.what()<<'\n';return 1;
    }
}
