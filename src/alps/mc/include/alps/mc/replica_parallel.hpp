// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "parallel.hpp"
#include <optional>

namespace alps::mc {
inline constexpr char replica_parallel_schema[]=R"toml(
[execution.parallel]
type = "string"
choices = ["chains", "replicas"]
default = "chains"
[execution.ranks_per_replica]
type = "int64"
min = 1
default = 1
)toml";

// Physical walkers have stable owners; the exchange permutation and each
// temperature's chronological statistics are replicated. Communicators are
// split lazily, after the runner has agreed on all inputs and checkpoints.
struct replica_parallel : collective {
    bool distributed=false;
    int ranks_per_replica=1;
#ifdef ALPS_HAVE_MPI
    std::optional<boost::mpi::communicator> worker;
#endif
    void configure(alps::run_configuration const& run) {
        distributed=run.execution["parallel"].as<std::string>()=="replicas";
        auto n=run.execution["ranks_per_replica"].as<uint64_t>();
        if (!n || n>uint64_t(size()) || size()%n || (!distributed && n!=1))
            throw std::invalid_argument("ranks_per_replica must divide the MPI size and requires parallel = replicas");
        ranks_per_replica=int(n);
    }
    bool owns(size_t id) const {return distributed || parallel::owns(id);}
    bool stopped(bool local) const {return distributed ? any(local) : local;}
    int teams() const {return distributed ? size()/ranks_per_replica : 1;}
    int owner(size_t walker) const {return int(walker%teams())*ranks_per_replica;}
    bool owns_walker(size_t walker) const {return !distributed || owner(walker)/ranks_per_replica==rank()/ranks_per_replica;}
    bool head() const {return !distributed || rank()%ranks_per_replica==0;}
    void initialize() {
#ifdef ALPS_HAVE_MPI
        if (distributed && ranks_per_replica>1 && !worker) worker=world.split(rank()/ranks_per_replica,rank());
#endif
    }
    template<class F> void update(F const& f) const {
        if (distributed) checked(f);else f();
    }
    void collect(std::vector<double>& values) const {
#ifdef ALPS_HAVE_MPI
        if (distributed && size()>1) {
            if (values.size()>INT_MAX) throw std::overflow_error("Replica sample exceeds MPI capacity");
            BOOST_MPI_CHECK_RESULT(MPI_Allreduce,(MPI_IN_PLACE,values.data(),int(values.size()),MPI_DOUBLE,MPI_SUM,world));
        }
#else
        (void)values;
#endif
    }
    template<class Simulation> auto make(alps::params const& p,size_t bins,size_t id) const {
        return std::make_unique<Simulation>(p,bins,id,*this);
    }
    template<class Chains> void synchronize(Chains& chains) const {
        if (distributed) for (auto& chain:chains) chain->synchronize();
        else parallel::synchronize(chains);
    }
    template<class Walkers> void synchronize_walkers(Walkers& walkers) const {
        if (!distributed) return;
        transport([&](auto& ar) {
            if (head()) for (size_t w=0;w<walkers.size();++w) if (owns_walker(w))
                ar["/walkers/"+std::to_string(w)]<<walkers[w];
        },[&](auto& ar,int source) {
            for (size_t w=0;w<walkers.size();++w) if (owner(w)==source)
                ar["/walkers/"+std::to_string(w)]>>walkers[w];
        });
    }
};
}
