// Copyright (C) 1997-2010 Synge Todo; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
//
// Replica exchange of heat-bath Ising walkers. Each walker may itself be
// decomposed over a team of MPI processes, so a ladder runs on
// teams x processes_per_walker ranks.
#include "../classical.hpp"
#include "../multiple/kernel.hpp"

// The spatial kernel draws every proposal from its team root's RNG, so a
// walker's trajectory does not depend on how many processes share it.
class spatial_walker {
    alps::random01 random_;
    spatial_ising_kernel kernel_;
public:
    spatial_walker(alps::params const& p,size_t offset):random_(p["SEED"].as<int>()+offset,p["RNG"].as<std::string>()),
        kernel_(alps::graph_helper<>(alps::make_deprecated_parameters(p)),p["J"].as<double>(),[&]{return random_();}) {}
    static auto names() {return classical_walker<1>::names();}
    void step(double beta) {kernel_.step(beta,[&]{return random_();});}
    std::vector<double> sample() const {return kernel_.sample();}
    void synchronize() {kernel_.synchronize();}
    void validate(size_t processes) const {kernel_.validate_partition(int(processes));}
#ifdef ALPS_HAVE_MPI
    void distribute(boost::mpi::communicator const& team) {kernel_.distribute(team);}
#endif
    void save(alps::hdf5::archive& ar) const {kernel_.save(ar); ar["rng"]<<random_;}
    void load(alps::hdf5::archive& ar) {
        auto random=random_;
        ar["rng"]>>random;
        if (random.name()!=random_.name()) throw std::invalid_argument("Checkpoint RNG changed");
        kernel_.load(ar); random_=std::move(random);
    }
};

inline constexpr char team_schema[]=R"toml(
[execution.processes_per_walker]
type = "int64"
default = 1
min = 1
max = 2147483647
)toml";

// Consecutive ranks form a team; walker w belongs to team w % teams. Only a
// team's root contributes samples and stores the walker's checkpoint.
struct team_parallel : alps::mc::replica_parallel {
#ifdef ALPS_HAVE_MPI
    static constexpr auto threading=ising_group::threading;
    std::optional<boost::mpi::communicator> team;
#endif
    size_t processes=1;
    void configure(alps::run_configuration const& run) {
        replica_parallel::configure(run);
        processes=run.execution["processes_per_walker"].as<int64_t>();
        if (processes>1 && !distributed)
            throw std::invalid_argument("execution.processes_per_walker requires execution.parallel = \"replicas\"");
        if (size_t(size())%processes)
            throw std::invalid_argument("The number of MPI processes must be a multiple of execution.processes_per_walker");
    }
    template<class Simulation> auto make(alps::params const& p,size_t bins,size_t id) const {
        return std::make_unique<Simulation>(p,bins,id,*this);
    }
    size_t teams() const {return size_t(size())/processes;}
    bool owns_walker(size_t w) const {return !distributed || w%teams()==size_t(rank())/processes;}
    bool stores_walker(size_t w) const {return owns_walker(w) && size_t(rank())%processes==0;}
    void collect(std::vector<double>& values) const {
        if (size_t(rank())%processes) std::fill(values.begin(),values.end(),0.);
        replica_parallel::collect(values);
    }
    template<class Walkers> void check(Walkers const& walkers) const {
        for (auto const& walker:walkers) if (walker) walker->validate(processes);
    }
    // Every rank reaches its first update or checkpoint together, after the
    // runner's collective preflight; only then is the world split.
    template<class Walkers> void bind(Walkers& walkers) {
#ifdef ALPS_HAVE_MPI
        if (processes==1) return;
        if (!team) team=world.split(rank()/int(processes));
        for (auto& walker:walkers) if (walker) walker->distribute(*team);
#else
        (void)walkers;
#endif
    }
};

int main(int argc,char** argv) {
    return classical_main<spatial_walker,team_parallel>(argc,argv,"exchange",
        {"ising; exchange","multiple parallel ising; exchange","ising"},team_schema);
}
