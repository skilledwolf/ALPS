/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "kernel.hpp"
#include "../single/native.hpp"
#include <array>

// All ranks cooperate on each physical chain; publication remains root-only.
struct spatial_group : ising_group {
    bool owns(size_t) const {return true;}
    bool stopped(bool local) const {return any(local);}
    template<class Runs> void verify(Runs const& runs,bool validate=false) const {
#ifdef ALPS_HAVE_MPI
        if (size()==1) return;
        std::vector<std::string> local,expected;
        checked([&] {
            local.push_back(validate ? "validate" : "execute");
            for (auto const& run:runs) local.push_back(alps::format_run_configuration(run));
            if (rank()==0) expected=local;
        });
        boost::mpi::broadcast(world,expected,0);
        checked([&] {
            if (local!=expected) throw std::invalid_argument("Spatial ranks require identical run configurations");
        });
        // Matching paths alone do not guarantee matching rank-local input
        // files. Compare bytes once before any physical collective, keeping
        // the existing schemas and bounded memory even for large checkpoints.
        std::set<std::string> inputs;
        checked([&] {for (auto const& run:runs) for (auto const& [key,value]:run.input)
            for (auto const& path:alps::run_paths(value)) inputs.insert(path);});
        for (auto const& path:inputs) {
            std::ifstream file(path,std::ios::binary);
            checked([&] {if (!file) throw std::runtime_error("Cannot read spatial input: "+path);});
            std::array<char,65536> actual{},reference{};
            bool different=false;
            for (;;) {
                int count=0;
                if (rank()==0) {file.read(reference.data(),reference.size());count=int(file.gcount());}
                boost::mpi::broadcast(world,count,0);
                if (!count) break;
                boost::mpi::broadcast(world,reference.data(),count,0);
                if (rank()!=0) {
                    file.read(actual.data(),count);
                    different=different || file.gcount()!=count || !std::equal(actual.begin(),actual.begin()+count,reference.begin());
                }
            }
            different=different || file.bad() || file.peek()!=std::char_traits<char>::eof();
            checked([&] {if (different) throw std::invalid_argument("Spatial ranks require identical input files: "+path);});
        }
#else
        (void)runs;(void)validate;
#endif
    }
    template<class Simulation> auto make(alps::params const& p,size_t bins,size_t id) const {
#ifdef ALPS_HAVE_MPI
        return std::make_unique<Simulation>(p,bins,id,world);
#else
        return std::make_unique<Simulation>(p,bins,id);
#endif
    }
    template<class Chains> void synchronize(Chains& chains) const {
        for (auto& chain:chains) chain->synchronize();
    }
};

int main(int argc,char** argv) {
    return ising_main<basic_ising<spatial_ising_kernel>,spatial_group>(argc,argv,"ising_multiple");
}
