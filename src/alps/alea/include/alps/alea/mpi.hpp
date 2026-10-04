// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once

#include <alps/alea/core.hpp>
#include <mpi.h>
#include <climits>
#include <stdexcept>
#include <string>

namespace alps::alea {

// Optional transport: link MPI::MPI_CXX alongside ALPS::statistics.
// The caller initializes MPI and keeps this borrowed intracommunicator alive.
// Construction and reductions are collective, in the same order on every rank.
class mpi_reducer final : public reducer {
public:
    explicit mpi_reducer(MPI_Comm comm, int root = 0) : comm_(comm), root_(root) {
        int inter;
        checked(MPI_Comm_test_inter(comm_, &inter));
        if (inter) throw std::invalid_argument("ALEA reductions require an intracommunicator");
        checked(MPI_Comm_rank(comm_, &rank_));
        checked(MPI_Comm_size(comm_, &size_));
        auto maximum = get_max(root_), minimum = -get_max(-int64_t(root_));
        if (minimum != maximum || minimum < 0 || maximum >= size_)
            throw std::invalid_argument("ALEA reduction root must agree on every rank and belong to the communicator");
    }

    reducer_setup get_setup() const override {
        return {size_t(rank_), uint64_t(size_), rank_ == root_};
    }
    int64_t get_max(int64_t value) const override {
        checked(MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT64_T, MPI_MAX, comm_));
        return value;
    }

    using reducer::reduce;
    void reduce(view<double> data) const override { sum(data, MPI_DOUBLE); }
    void reduce(view<int32_t> data) const override { sum(data, MPI_INT32_T); }
    void reduce(view<int64_t> data) const override { sum(data, MPI_INT64_T); }
    void reduce(view<uint64_t> data) const override { sum(data, MPI_UINT64_T); }
    void commit() const override {} // Transfers finish before reduce returns.

private:
    static void checked(int error) {
        if (error != MPI_SUCCESS)
            throw std::runtime_error("ALEA MPI operation failed: " + std::to_string(error));
    }
    template<class T> void sum(view<T> data, MPI_Datatype datatype) const {
        int64_t count = data.size() <= INT_MAX && (!data.size() || data.data())
                      ? int64_t(data.size()) : -1;
        auto maximum = get_max(count), minimum = -get_max(-count);
        if (minimum < 0)
            throw std::invalid_argument("ALEA MPI buffers must be non-null and contain at most INT_MAX elements");
        if (minimum != maximum) throw size_mismatch();
        if (!count) return;
        checked(MPI_Reduce(rank_ == root_ ? MPI_IN_PLACE : data.data(),
                           rank_ == root_ ? data.data() : nullptr,
                           int(count), datatype, MPI_SUM, root_, comm_));
    }
    MPI_Comm comm_;
    int root_, rank_, size_;
};

}
