// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Compile against the isolated ALPSCore reference SDK, never ALPS headers.
#include <alps/dictionary.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    // A command-line run: Core records its INI text, origins and help with the values.
    const char* command[] = {"core-fixture", "--L=8", "--model=heisenberg"};
    alps::params p(3, command);
    p.description("ALPSCore 2.3.3 params checkpoint fixture")
        .define<int>("L", 4, "linear lattice size")
        .define<double>("T", 2.5, "temperature")
        .define<bool>("flag", true, "a Boolean switch")
        .define<std::string>("model", "ising", "model name")
        .define<double>("unset", "declared without a value");
    p["off"] = false;
    p["negative"] = -3;
    p["unsigned"] = 4294967295u;
    p["wide"] = 1099511627776L;
    p["unsigned_wide"] = 18446744073709551615UL;
    p["single"] = 1.25F;
    p["text with space"] = std::string("a,b");
    p["a/b&c"] = 7;
    p["reals"] = std::vector<double>{1.5, -2.0};
    p["integers"] = std::vector<int>{1, 2, 3};
    p["flags"] = std::vector<bool>{true, false, true};
    p["words"] = std::vector<std::string>{"x", "y z"};
    p["empty"] = std::vector<double>{};
    alps::params_ns::dictionary d;
    d["count"] = 3;
    d["ratio"] = 0.5;
    alps::hdf5::archive archive(argv[1], "w");
    archive["/parameters"] << p;
    archive["/dictionary"] << d;
}
