// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <complex>
#include <functional>
#include <iostream>
#include <limits>

namespace {
void require(const boost::mpi::communicator &comm, bool valid, const char *message) {
    if (!boost::mpi::all_reduce(comm, valid, std::logical_and<bool>()))
        throw std::runtime_error(message);
}
alps::params payload(int rank) {
    alps::params parameters;
    parameters["source_rank"] = rank;
    parameters["flag"] = rank % 2 == 0;
    parameters["signed"] = std::int64_t(9007199254740993LL);
    parameters["unsigned"] = std::numeric_limits<std::uint64_t>::max();
    parameters["real"] = 0.125 * (rank + 1);
    parameters["complex"] = std::complex<double>(rank + 0.5, -2.0);
    parameters["string"] = std::string("embedded\0nul", 12);
    parameters["bool_vector"] = std::vector<bool>{true, false, rank % 2 == 0};
    parameters["signed_vector"] = std::vector<std::int64_t>{
        std::numeric_limits<std::int64_t>::min(), 9007199254740993LL,
        std::numeric_limits<std::int64_t>::max()};
    parameters["unsigned_vector"] = std::vector<std::uint64_t>{
        0, 9007199254740993ULL, std::numeric_limits<std::uint64_t>::max()};
    parameters["real_vector"] = std::vector<double>{-0.25, 0.0, 0.5};
    parameters["complex_vector"] = std::vector<std::complex<double>>{{1.0, -2.0}, {0.0, 0.5}};
    parameters["string_vector"] = std::vector<std::string>{"a,b", "", "root " + std::to_string(rank)};
    parameters["empty_bool"] = std::vector<bool>{};
    parameters["empty_signed"] = std::vector<std::int64_t>{};
    parameters["empty_unsigned"] = std::vector<std::uint64_t>{};
    parameters["empty_real"] = std::vector<double>{};
    parameters["empty_complex"] = std::vector<std::complex<double>>{};
    parameters["empty_string"] = std::vector<std::string>{};
    parameters["unset"];
    return parameters;
}
}
int main(int argc, char **argv) {
    boost::mpi::environment environment(argc, argv);
    boost::mpi::communicator world;
    try {
        require(world, world.size() >= 2, "params MPI contract requires at least two ranks");
        // Both roots exercise replacement of rank-specific content, not just filling an empty map.
        for (const int root : {0, world.size() - 1}) {
            auto received = payload(world.rank());
            if (world.rank() != root)
                received["receiver_only_" + std::to_string(world.rank())] = world.rank();
            const auto before = received;
            const auto expected = payload(root);
            received.broadcast(world, root);
            require(world, received == expected, "broadcast changed type/value or kept stale receiver keys");
            require(world, before["source_rank"].as<int>() == world.rank(), "broadcast mutated an owning copy");
            require(world, received["signed"].as<std::int64_t>() == 9007199254740993LL,
                    "int64 precision was lost in transport");
            require(world, received["unsigned"].as<std::uint64_t>() ==
                           std::numeric_limits<std::uint64_t>::max(), "uint64 range was lost in transport");
            require(world, received.exists<std::vector<bool>>("empty_bool") &&
                           received.exists<std::vector<std::int64_t>>("empty_signed") &&
                           received.exists<std::vector<std::uint64_t>>("empty_unsigned") &&
                           received.exists<std::vector<double>>("empty_real") &&
                           received.exists<std::vector<std::complex<double>>>("empty_complex") &&
                           received.exists<std::vector<std::string>>("empty_string"),
                    "empty vector type changed in transport");
            require(world, !received.exists("unset") && received.find("unset") != received.end(),
                    "unset entry changed in transport");
            bool lossy_conversion_rejected = false;
            try { received["signed"].as<double>(); }
            catch (const alps::params_ns::exception::value_mismatch &) { lossy_conversion_rejected = true; }
            require(world, lossy_conversion_rejected, "transported integer allowed a lossy conversion");
            auto copied = received;
            copied["string_vector"] = std::vector<std::string>{"changed"};
            auto extracted = received["complex_vector"].as<std::vector<std::complex<double>>>();
            extracted[0] = {42.0, 0.0};
            require(world, received == expected && before["source_rank"].as<int>() == world.rank(),
                    "copy or extracted vector retained aliases");
        }
        auto cleared = payload(world.rank());
        if (world.rank() == 0) cleared = alps::params{};
        cleared.broadcast(world, 0);
        require(world, cleared.empty(), "empty root did not replace a nonempty receiver map");
        if (world.rank() == 0) std::cout << "Typed params MPI contracts passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "params MPI rank " << world.rank() << ": " << error.what() << '\n';
        world.abort(1);
        return 1;
    }
}
