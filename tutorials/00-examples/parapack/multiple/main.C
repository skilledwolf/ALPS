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

// All ranks cooperate on each physical chain; publication remains root-only.
struct spatial_group : alps::mc::collective {
#ifdef ALPS_HAVE_MPI
    static constexpr auto threading=ising_group::threading;
#endif
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
