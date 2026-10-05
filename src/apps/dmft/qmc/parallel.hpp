// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/config.h>
#include <stdexcept>
#include <string>
#ifdef ALPS_HAVE_MPI
#include <alps/alea/mpi.hpp>
#include <climits>
#endif

namespace alps::solvers {
// An MPI-enabled library remains usable from an ordinary serial caller.
// The caller owns MPI initialization/finalization when using multiple ranks.
struct parallel {
    int rank=0, size=1;
    parallel() {
#ifdef ALPS_HAVE_MPI
        int initialized=0, finalized=0;
        MPI_Initialized(&initialized); MPI_Finalized(&finalized);
        if (initialized && !finalized) {
            MPI_Comm_rank(MPI_COMM_WORLD,&rank);
            MPI_Comm_size(MPI_COMM_WORLD,&size);
        }
#endif
    }
    void agree(std::string message) const {
#ifdef ALPS_HAVE_MPI
        if (size>1) {
            int origin=message.empty() ? size : rank;
            MPI_Allreduce(MPI_IN_PLACE,&origin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
            if (origin==size) return;
            if (message.size()>INT_MAX) message="Solver error exceeds MPI message limit";
            int length=static_cast<int>(message.size());
            MPI_Bcast(&length,1,MPI_INT,origin,MPI_COMM_WORLD);
            message.resize(length);
            MPI_Bcast(message.data(),length,MPI_CHAR,origin,MPI_COMM_WORLD);
        }
#endif
        if (!message.empty()) throw std::runtime_error(message);
    }
    template<class F> void checked(F operation) const {
        std::string error;
        try { operation(); }
        catch (std::exception const& e) { error=*e.what() ? e.what() : "Solver operation failed"; }
        catch (...) { error="Unknown solver failure"; }
        agree(std::move(error));
    }
    template<class Simulation> auto collect(Simulation const& simulation) const {
#ifdef ALPS_HAVE_MPI
        if (size>1) {
            alps::alea::mpi_reducer reducer(MPI_COMM_WORLD);
            return simulation.collect_results(&reducer);
        }
#endif
        return simulation.collect_results();
    }
};
}
