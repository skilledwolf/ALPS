// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Bespoke adaptation of local error agreement in ALPS reconciliation
// 90a36c0d4, src/alps/mc/include/alps/mc/parallel.hpp. The callback contract is
// deliberately restricted to local work; this is not MPI fault recovery.
#pragma once

#include <mpi.h>
#include <array>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>

namespace alps::detail {

class collective_error : public std::runtime_error {
public:
    collective_error(int rank, const std::string &message)
        : std::runtime_error(message), rank_(rank) {}
    int failing_rank() const noexcept { return rank_; }
private:
    int rank_;
};

// The caller initializes MPI and keeps this borrowed intracommunicator alive.
// any(), all() and checked_local() are collective and must be called in the
// same order by every rank. No signal handlers, threads or error handlers are
// installed. MPI communication failures or missing participants are outside
// this contract; the communicator's existing MPI error handler is respected.
class mpi_group {
public:
    explicit mpi_group(MPI_Comm communicator) : communicator_(communicator) {
        int inter;
        checked(MPI_Comm_test_inter(communicator_, &inter));
        if (inter) throw std::invalid_argument("MC consensus requires an intracommunicator");
        checked(MPI_Comm_rank(communicator_, &rank_));
        checked(MPI_Comm_size(communicator_, &size_));
    }

    int rank() const noexcept { return rank_; }
    int size() const noexcept { return size_; }
    bool any(bool local) const { return agree(local, MPI_MAX); }
    bool all(bool local) const { return agree(local, MPI_MIN); }

    // Execute local work, then report the lowest failing rank on every rank.
    // The callback MUST NOT call MPI, including through this group or another
    // library. Wrapping a callback that diverges within MPI cannot prevent a
    // deadlock. Local mutations are not rolled back; publication belongs in a
    // later phase after checked_local() succeeds everywhere.
    //
    // Diagnostics are bounded to 1023 bytes before transport, avoiding dynamic
    // allocation while recording a failure. Each rank that threw retains its
    // own exception nested inside the shared collective_error. Successful
    // ranks receive only the agreed diagnostic. After catching the error on
    // all ranks, callers may use the group again with matching operations.
    template<class Operation> void checked_local(Operation &&operation) const {
        std::exception_ptr local_error;
        std::array<char, 1024> diagnostic{};
        try {
            std::forward<Operation>(operation)();
        } catch (...) {
            local_error = std::current_exception();
            try {
                std::rethrow_exception(local_error);
            } catch (const std::exception &error) {
                copy_message(diagnostic, error.what());
            } catch (...) {
                copy_message(diagnostic, "non-standard exception");
            }
        }
        int failed = local_error ? rank_ : size_;
        checked(MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MIN, communicator_));
        if (failed == size_) return;
        checked(MPI_Bcast(diagnostic.data(), int(diagnostic.size()), MPI_CHAR, failed, communicator_));
        const auto message = "MC local operation failed on rank " + std::to_string(failed)
                           + ": " + diagnostic.data();
        if (local_error) {
            try { std::rethrow_exception(local_error); }
            catch (...) { std::throw_with_nested(collective_error(failed, message)); }
        }
        throw collective_error(failed, message);
    }

private:
    static void copy_message(std::array<char, 1024> &destination, const char *message) noexcept {
        if (!message || !*message) message = "local callback raised an exception";
        for (std::size_t i = 0; i + 1 < destination.size() && message[i]; ++i)
            destination[i] = message[i];
    }
    static void checked(int error) {
        if (error != MPI_SUCCESS)
            throw std::runtime_error("MC MPI operation failed: " + std::to_string(error));
    }
    bool agree(bool local, MPI_Op operation) const {
        int value = local;
        checked(MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, operation, communicator_));
        return value != 0;
    }
    MPI_Comm communicator_;
    int rank_, size_;
};

} // namespace alps::detail
