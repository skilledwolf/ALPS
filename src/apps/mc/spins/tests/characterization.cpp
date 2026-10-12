// SPDX-License-Identifier: MIT
// Characterize how spinmc interprets its parameters. Workers come from the
// factory the scheduler uses, but each case runs a fixed number of sweeps
// instead of the scheduler's wall-clock completion checks, so the results do
// not depend on timing. The checks pin behavior that refactoring must keep:
// defaults and precedence, update selection, the parameter set saved with a
// run, and the first error reported for invalid input.
//
// Set ALPS_SPINMC_DUMP_DIR to write each case's complete worker state
// (parameters, RNG state and measurements) after its sweeps. Comparing the
// dumps of two builds shows whether a change alters the Monte Carlo
// trajectory; the dumps are platform-specific and are not committed.
#include "factory.h"

#include <alps/hdf5/archive.hpp>
#include <alps/osiris/comm.h>
#include <alps/parameter.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

#include <cstdlib>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

std::unique_ptr<alps::scheduler::Worker> make_worker(const alps::Parameters& p) {
    const alps::ProcessList where(1, alps::local_process());
    return std::unique_ptr<alps::scheduler::Worker>(SpinFactory().make_worker(where, p, 0));
}

std::string construction_error(const std::string& text) {
    try {
        make_worker(parameters(text));
    } catch (const std::exception& error) {
        return error.what();
    }
    return "";
}

std::set<std::string> names(const alps::Parameters& p) {
    std::set<std::string> out;
    for (const auto& entry : p)
        out.insert(entry.key());
    return out;
}

struct Case {
    const char* name;
    const char* input;
    int sweeps;                          // dostep() calls before the checks
    bool cluster;                        // whether cluster updates were selected
    double work_done;                    // negative: depends on cluster sizes
    std::set<std::string> saved;         // parameter names saved with the run
};

// Local updates thermalize after exactly THERMALIZATION sweeps, so their
// completed fraction is (sweeps - THERMALIZATION) / SWEEPS. The default
// thermalization is SWEEPS/10 in integer arithmetic: 95/10 = 9.
const std::vector<Case> cases = {
    {"IsingChainLocalDefaultThermalization",
     R"p(LATTICE="chain lattice";L=8;MODEL="Ising";J=1;T=1.5;SWEEPS=95;UPDATE="local";SEED=11)p",
     30, false, (30. - 9.) / 95.,
     {"LATTICE", "L", "MODEL", "J", "T", "SWEEPS", "UPDATE", "SEED"}},
    {"IsingSquareAutomaticCluster",
     R"p(LATTICE="square lattice";L=4;MODEL="Ising";J=1;T=2.5;SWEEPS=50;THERMALIZATION=5;SEED=3)p",
     25, true, -1.,
     {"LATTICE", "L", "MODEL", "J", "T", "SWEEPS", "THERMALIZATION", "SEED", "UPDATE"}},
    {"IsingTriangularFrustratedAutomaticLocal",
     R"p(LATTICE="triangular lattice";L=6;MODEL="Ising";J=1;T=2;SWEEPS=40;THERMALIZATION=4;SEED=4)p",
     20, false, (20. - 4.) / 40.,
     {"LATTICE", "L", "MODEL", "J", "T", "SWEEPS", "THERMALIZATION", "SEED", "UPDATE"}},
    {"IsingTriangularUnfrustratedAutomaticCluster",
     R"p(LATTICE="triangular lattice";L=6;MODEL="Ising";J=-1;T=1;SWEEPS=40;THERMALIZATION=4;SEED=4)p",
     20, true, -1.,
     {"LATTICE", "L", "MODEL", "J", "T", "SWEEPS", "THERMALIZATION", "SEED", "UPDATE"}},
    {"HeisenbergLadderBetaPerBondCluster",
     R"p(LATTICE="ladder";L=4;MODEL="Heisenberg";J0=1;J1=0.5;beta=2;T=100;SWEEPS=40;THERMALIZATION=4;UPDATE="cluster";SEED=5)p",
     20, true, -1.,
     {"LATTICE", "L", "MODEL", "J0", "J1", "beta", "T", "SWEEPS", "THERMALIZATION", "UPDATE",
      "SEED"}},
    {"HeisenbergFieldExpressionsLocal",
     R"p(LATTICE="chain lattice";L=6;MODEL="Heisenberg";J=1;T=0.8;h="0.2*J";g="2*J";S="J/2";SWEEPS=60;THERMALIZATION=6;UPDATE="local";SEED=9)p",
     25, false, (25. - 6.) / 60.,
     {"LATTICE", "L", "MODEL", "J", "T", "h", "g", "S", "SWEEPS", "THERMALIZATION", "UPDATE",
      "SEED"}},
    {"XYQuantumAnisotropicAutomaticLocal",
     R"p(LATTICE="chain lattice";L=6;MODEL="XY";J="1 0.5";D="0.1 0.2";CONVENTION="quantum";S=1.5;T=0.7;SWEEPS=50;THERMALIZATION=5;SEED=6)p",
     20, false, (20. - 5.) / 50.,
     {"LATTICE", "L", "MODEL", "J", "D", "CONVENTION", "S", "T", "SWEEPS", "THERMALIZATION", "SEED",
      "UPDATE"}},
    {"O4ChainAutomaticCluster",
     R"p(LATTICE="chain lattice";L=6;MODEL="O(4)";J=1;T=1;SWEEPS=30;THERMALIZATION=3;SEED=2)p",
     15, true, -1.,
     {"LATTICE", "L", "MODEL", "J", "T", "SWEEPS", "THERMALIZATION", "SEED", "UPDATE"}},
};

class Characterization : public ::testing::TestWithParam<Case> {};

TEST_P(Characterization, InterpretsParametersAndSavesTheRun) {
    const Case& c = GetParam();
    auto worker = make_worker(parameters(c.input));
    ASSERT_TRUE(worker);
    for (int i = 0; i < c.sweeps; ++i)
        worker->dostep();

    const auto& run = dynamic_cast<const alps::scheduler::MCRun&>(*worker);
    EXPECT_EQ(run.get_measurements().has("Cluster size"), c.cluster);
    if (c.work_done >= 0.)
        EXPECT_DOUBLE_EQ(worker->work_done(), c.work_done);

    alps::testing::TemporaryDirectory directory;
    const auto file = (directory.path() / "run.h5").string();
    {
        alps::hdf5::archive archive(file, "w");
        worker->save(archive);
    }
    alps::hdf5::archive archive(file);
    alps::Parameters saved;
    archive["/parameters"] >> saved;
    EXPECT_EQ(names(saved), c.saved);

    if (const char* dump = std::getenv("ALPS_SPINMC_DUMP_DIR")) {
        alps::hdf5::archive out(std::string(dump) + "/" + c.name + ".h5", "w");
        worker->save(out);
    }
}

INSTANTIATE_TEST_SUITE_P(Spinmc, Characterization, ::testing::ValuesIn(cases),
                         [](const ::testing::TestParamInfo<Case>& info) {
                             return std::string(info.param.name);
                         });

// An absent UPDATE is inserted while the worker is constructed and saved
// with the run. Its empty value round-trips through HDF5 as the number 0.
TEST(SpinmcCharacterization, SavesAnInsertedEmptyUpdate) {
    auto worker = make_worker(parameters(
        R"p(LATTICE="chain lattice";L=4;MODEL="Ising";J=1;T=1;SWEEPS=10;SEED=1)p"));
    alps::testing::TemporaryDirectory directory;
    const auto file = (directory.path() / "run.h5").string();
    {
        alps::hdf5::archive archive(file, "w");
        worker->save(archive);
    }
    alps::hdf5::archive archive(file);
    alps::Parameters saved;
    archive["/parameters"] >> saved;
    ASSERT_TRUE(saved.defined("UPDATE"));
    EXPECT_EQ(std::string(saved["UPDATE"]), "0");
}

const std::string chain = R"p(LATTICE="chain lattice";L=4;MODEL="Heisenberg";J=1;SEED=1;)p";

TEST(SpinmcCharacterization, RequiresTemperatureBeforeSweeps) {
    EXPECT_EQ(construction_error(chain + "SWEEPS=10"), "parameter T not defined");
    EXPECT_EQ(construction_error(chain), "parameter T not defined");
    EXPECT_EQ(construction_error(chain + "T=1"), "parameter SWEEPS not defined");
}

TEST(SpinmcCharacterization, RejectsInconsistentFieldSettings) {
    EXPECT_EQ(construction_error(chain + R"p(T=1;SWEEPS=10;h=0.5;CONVENTION="quantum")p"),
              "Cannot use quantum convention with nonzero magnetic field .");
    EXPECT_EQ(construction_error(chain + R"p(T=1;SWEEPS=10;h=0.5;UPDATE="cluster")p"),
              "Illegal update type cluster with magnetic field not set to zero.");
}

TEST(SpinmcCharacterization, ReportsCouplingErrorsBeforeUnknownUpdates) {
    EXPECT_EQ(construction_error(chain + R"p(T=1;SWEEPS=10;UPDATE="fast")p"),
              "Illegal update type fast");
    EXPECT_EQ(construction_error(chain + R"p(T=1;SWEEPS=10;D="1 2 3 4 5 6 7 8 9")p"),
              "Invalid matrix for D");
    EXPECT_EQ(construction_error(chain + R"p(T=1;SWEEPS=10;D="1 2 3 4 5 6 7 8 9";UPDATE="fast")p"),
              "Invalid matrix for D");
}

} // namespace
