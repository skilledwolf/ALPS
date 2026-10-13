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
#include <boost/core/demangle.hpp>
#include <gtest/gtest.h>

#include <cstdlib>
#include <algorithm>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
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
    // Random couplings draw from the disorder generator, so these dumps change
    // if the couplings are resolved in a different order: site types before
    // bond types, S before D within a site type, and bonds in edge order.
    {"HeisenbergLadderRandomCouplingsCluster",
     R"p(LATTICE="ladder";L=4;MODEL="Heisenberg";J0="random()";J1="0.5*random()";S="0.5+random()";DISORDERSEED=7;T=1;SWEEPS=40;THERMALIZATION=4;SEED=12)p",
     20, true, -1.,
     {"LATTICE", "L", "MODEL", "J0", "J1", "S", "DISORDERSEED", "T", "SWEEPS", "THERMALIZATION",
      "SEED", "UPDATE"}},
    {"HeisenbergTwoBandRandomSiteTermsLocal",
     R"p(LATTICE="2 band chain lattice";L=4;MODEL="Heisenberg";J0="random()";J1="0.5+random()";J2=-1;S0="0.5+random()";S1="1+random()";D0="random() 0.2 0.3";D1="random()";DISORDERSEED=7;T=1;SWEEPS=40;THERMALIZATION=4;SEED=13)p",
     20, false, (20. - 4.) / 40.,
     {"LATTICE", "L", "MODEL", "J0", "J1", "J2", "S0", "S1", "D0", "D1", "DISORDERSEED", "T",
      "SWEEPS", "THERMALIZATION", "SEED", "UPDATE"}},
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

// Site types are resolved before bonds, site type 0 before site type 1, and
// S before D within a site type.
TEST(SpinmcCharacterization, ReportsSiteTermErrorsInTypeOrderBeforeBondErrors) {
    const std::string two_band =
        R"p(LATTICE="2 band chain lattice";L=4;MODEL="Heisenberg";T=1;SWEEPS=10;SEED=1;)p";
    const std::string bad_d0 = R"p(D0="1 2 3 4 5 6 7 8 9";)p";
    const std::string bad_d1 = R"p(D1="1 2 3 4 5 6 7 8 9";)p";
    const std::string bad_s1 = R"p(S1="1+";)p";
    const std::string bad_j0 = R"p(J0="1 2";)p";
    const std::string bad_expression = "Illegal term in expression";
    EXPECT_EQ(construction_error(two_band + "J=1;" + bad_s1), bad_expression);
    EXPECT_EQ(construction_error(two_band + "J=1;" + bad_d1 + bad_s1), bad_expression);
    EXPECT_EQ(construction_error(two_band + "J=1;" + bad_d0 + bad_s1), "Invalid matrix for D");
    EXPECT_EQ(construction_error(two_band + "J=1;" + bad_d1), "Invalid matrix for D");
    EXPECT_EQ(construction_error(two_band + bad_j0 + R"p(D0="1 2 3";)p"), "invalid input");
    EXPECT_EQ(construction_error(two_band + bad_j0 + bad_d1), "Invalid matrix for D");
}

// The simulation class the factory builds, without whitespace. The moment
// headers cannot be included in a second translation unit (ising.h defines a
// static member), so compare type names. XYMoment names ONMoment<2>.
std::string built(const std::string& text) {
    const auto worker = make_worker(parameters(text));
    const alps::scheduler::Worker& built_worker = *worker;
    std::string name = boost::core::demangle(typeid(built_worker).name());
    name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
    return name;
}

const std::string heisenberg =
    R"p(LATTICE="chain lattice";L=4;MODEL="Heisenberg";T=1;SWEEPS=10;SEED=1;)p";
const std::string xy = R"p(LATTICE="chain lattice";L=4;MODEL="XY";T=1;SWEEPS=10;SEED=1;)p";
const std::string ising = R"p(LATTICE="chain lattice";L=4;MODEL="Ising";T=1;SWEEPS=10;SEED=1;)p";

// The factory chooses the matrix type from the largest number of
// space-separated values among J, J', J0..J3, D, D' and D0..D3.
TEST(SpinmcCharacterization, ChoosesTheMatrixTypeFromTheLongestCouplingText) {
    EXPECT_EQ(built(heisenberg + "J=1"), "SpinSim<ONMoment<3>,MIdMatrix<double,3>>");
    EXPECT_EQ(built(heisenberg + R"p(J="1 2 3")p"), "SpinSim<ONMoment<3>,DiagMatrix<double,3>>");
    EXPECT_EQ(built(heisenberg + R"p(J="1 0.5 1 0 0 1")p"),
              "SpinSim<ONMoment<3>,SquareMatrix<double,3>>");
    EXPECT_EQ(built(heisenberg + R"p(J="1 0 0 0 1 0 0 0 1")p"),
              "SpinSim<ONMoment<3>,SquareMatrix<double,3>>");
    EXPECT_EQ(built(heisenberg + R"p(J=1;D="0.1 0.2 0.3")p"),
              "SpinSim<ONMoment<3>,DiagMatrix<double,3>>");

    EXPECT_EQ(built(xy + "J=1"), "SpinSim<ONMoment<2>,MIdMatrix<double,2>>");
    EXPECT_EQ(built(xy + R"p(J="1 2")p"), "SpinSim<ONMoment<2>,DiagMatrix<double,2>>");
    EXPECT_EQ(built(xy + R"p(J="1 0.5 1")p"), "SpinSim<ONMoment<2>,SquareMatrix<double,2>>");
    EXPECT_EQ(built(xy + R"p(J="1 0 0 1")p"), "SpinSim<ONMoment<2>,SquareMatrix<double,2>>");

    EXPECT_EQ(built(ising + "J=1"), "SpinSim<IsingMoment,MIdMatrix<double,1>>");
    EXPECT_EQ(built(R"p(LATTICE="chain lattice";L=4;MODEL="O(4)";J=1;T=1;SWEEPS=10;SEED=1)p"),
              "SpinSim<ONMoment<4>,MIdMatrix<double,4>>");
    EXPECT_EQ(built(R"p(LATTICE="chain lattice";L=4;MODEL="Potts";q=3;J=1;T=1;SWEEPS=10;SEED=1)p"),
              "SpinSim<PottsMoment<3u>,MIdMatrix<double,2>>");
}

// J' and D' take part in the choice although no bond or site reads them, and
// J4 is read for bond type 4 but takes no part in the choice.
TEST(SpinmcCharacterization, CountsPrimedCouplingsButNotHigherBondTypes) {
    EXPECT_EQ(built(heisenberg + R"p(J=1;J'="1 2 3")p"),
              "SpinSim<ONMoment<3>,DiagMatrix<double,3>>");
    EXPECT_EQ(built(heisenberg + R"p(J=1;D'="1 2 3")p"),
              "SpinSim<ONMoment<3>,DiagMatrix<double,3>>");
    EXPECT_EQ(built(heisenberg + R"p(J=1;J4="1 2 3")p"),
              "SpinSim<ONMoment<3>,MIdMatrix<double,3>>");
    // Values are split at spaces only, so an expression with spaces counts
    // as several values: three here, which selects the diagonal type.
    EXPECT_EQ(construction_error(heisenberg + R"p(J="2 * 0.5")p"), "Illegal term in expression");
}

// The model library names the coupling of bond type 1 J'. A ladder's rungs
// have bond type 1, so J' without J1 is rejected instead of giving the rungs
// J. A chain has no such bonds and still accepts an unused J'.
TEST(SpinmcCharacterization, RejectsPrimedCouplingsThatWouldBeIgnored) {
    const std::string ladder =
        R"p(LATTICE="ladder";L=4;MODEL="Heisenberg";T=1;SWEEPS=10;SEED=1;)p";
    EXPECT_EQ(construction_error(ladder + "J=1;J'=0.5"),
              R"(J' is not read: bond type 1 takes J1, falling back to J. Set J1, for example J1="J'", or remove J')");
    EXPECT_EQ(built(ladder + "J=1;J1=0.5;J'=0.5"), "SpinSim<ONMoment<3>,MIdMatrix<double,3>>");
    EXPECT_EQ(built(ladder + R"p(J=1;J'=0.5;J1="J'")p"),
              "SpinSim<ONMoment<3>,MIdMatrix<double,3>>");
    EXPECT_EQ(built(heisenberg + "J=1;J'=0.5"), "SpinSim<ONMoment<3>,MIdMatrix<double,3>>");
}

// An unsupported number of values is reported with the model and every
// counted parameter that has that many values. Absent parameters are skipped.
TEST(SpinmcCharacterization, ReportsUnsupportedModelsAndMatrixSizes) {
    EXPECT_EQ(construction_error(ising + R"p(J="1 2")p"),
              R"(invalid parameters: cannot produce model Ising from 2 values in J="1 2")");
    EXPECT_EQ(construction_error(ising + R"p(J="1 2";D=0)p"),
              R"(invalid parameters: cannot produce model Ising from 2 values in J="1 2")");
    EXPECT_EQ(construction_error(heisenberg + R"p(J="1 2")p"),
              R"(invalid parameters: cannot produce model Heisenberg from 2 values in J="1 2")");
    EXPECT_EQ(construction_error(heisenberg + R"p(J=1;D="1 2")p"),
              R"(invalid parameters: cannot produce model Heisenberg from 2 values in D="1 2")");
    EXPECT_EQ(construction_error(xy + R"p(J="1 2 3 4 5")p"),
              R"(invalid parameters: cannot produce model XY from 5 values in J="1 2 3 4 5")");
    EXPECT_EQ(construction_error(heisenberg + R"p(J=1;J'="1 2")p"),
              R"(invalid parameters: cannot produce model Heisenberg from 2 values in J'="1 2")");
    EXPECT_EQ(construction_error(heisenberg + R"p(J0="1 2";D1="3 4")p"),
              "invalid parameters: cannot produce model Heisenberg from 2 values in "
              R"(J0="1 2" D1="3 4")");
    EXPECT_EQ(construction_error(ising),
              "invalid parameters: cannot produce model Ising without J, J', J0-J3, D, D' or D0-D3");
    // O(4) and Potts ignore the count. An O(4) coupling with two values
    // fails only when the coupling is read.
    EXPECT_EQ(construction_error(
                  R"p(LATTICE="chain lattice";L=4;MODEL="O(4)";J="1 2";T=1;SWEEPS=10;SEED=1)p"),
              "invalid input");
    EXPECT_EQ(construction_error(
                  R"p(LATTICE="chain lattice";L=4;MODEL="Potts";J=1;T=1;SWEEPS=10;SEED=1)p"),
              "parameter q not defined");
    EXPECT_EQ(construction_error(
                  R"p(LATTICE="chain lattice";L=4;MODEL="Potts";q=5;J=1;T=1;SWEEPS=10;SEED=1)p"),
              "5-state Potts model not implemented in file factory.C");
    EXPECT_EQ(construction_error(R"p(LATTICE="chain lattice";L=4;J=1;T=1;SWEEPS=10;SEED=1)p"),
              "parameter MODEL not defined");
    EXPECT_EQ(construction_error(
                  R"p(LATTICE="chain lattice";L=4;MODEL="Clock";J=1;T=1;SWEEPS=10;SEED=1)p"),
              "Clock model not implemented in file factory.C");
}

} // namespace
