// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "parallel.hpp"

namespace alps::mc {
inline constexpr char replica_parallel_schema[]=R"toml(
[execution.parallel]
type = "string"
choices = ["chains", "replicas"]
default = "chains"
)toml";

// Physical walkers have stable owners; the exchange permutation and each
// temperature's chronological statistics are replicated. Construction and
// checkpoint loading stay local until the runner's collective preflight.
struct replica_parallel : collective {
    bool distributed=false;
    void configure(alps::run_configuration const& run) {
        distributed=run.execution["parallel"].as<std::string>()=="replicas";
    }
    bool owns(size_t id) const {return distributed || parallel::owns(id);}
    bool stopped(bool local) const {return distributed ? any(local) : local;}
    int owner(size_t walker) const {return int(walker%size());}
    bool owns_walker(size_t walker) const {return !distributed || owner(walker)==rank();}
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
    template<class T> void broadcast(T& value,size_t walker) const {
#ifdef ALPS_HAVE_MPI
        if (distributed) boost::mpi::broadcast(world,value,owner(walker));
#else
        (void)value;(void)walker;
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
            for (size_t w=0;w<walkers.size();++w) if (owns_walker(w))
                ar["/walkers/"+std::to_string(w)]<<walkers[w];
        },[&](auto& ar,int source) {
            for (size_t w=0;w<walkers.size();++w) if (owner(w)==source)
                ar["/walkers/"+std::to_string(w)]>>walkers[w];
        });
    }
};
}
