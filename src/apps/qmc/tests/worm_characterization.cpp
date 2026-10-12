// SPDX-License-Identifier: MIT
// Characterize the existing worm/QMCRun parameter contract before extracting it.
#include "../worms/WRun.h"

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/osiris/xdrdump.h>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

namespace {

alps::Parameters worm_parameters() {
  alps::Parameters p;
  p["LATTICE"] = "square lattice";
  p["L"] = 2;
  p["MODEL"] = "boson Hubbard";
  p["Nmax"] = 2;
  p["NONLOCAL"] = 0;
  p["U"] = 1;
  p["mu"] = 0.5;
  p["t"] = 0.1;
  p["T"] = 2;
  p["THERMALIZATION"] = 2;
  p["SWEEPS"] = 8;
  p["SEED"] = 12345;
  p["DISORDERSEED"] = 17;
  return p;
}

class WormProbe : public WRun {
public:
  explicit WormProbe(const alps::Parameters& p)
      : WRun(alps::ProcessList(1, alps::Process(0)), p, 0) {}
  double inverse_temperature() const { return beta; }
};

void start_without_log(WormProbe& run) {
  // Invoke the application's start hook without scheduler timestamps in dumps.
  static_cast<alps::scheduler::Worker&>(run).start();
}

void dump_reference(const std::string& name, const WormProbe& run) {
  const char* directory = std::getenv("ALPS_QMC_DUMP_DIR");
  if (!directory) return;
  std::filesystem::create_directories(directory);
  const auto prefix = std::filesystem::path(directory) / ("worm_" + name);
  // Keep the raw solver payload for auditing. The worker wrapper also contains
  // scheduler log data, so it is deliberately excluded from these references.
  {
    alps::OXDRFileDump dump(prefix.string() + ".xdr");
    run.save(dump);
  }
  int version = 0, num_kinks = 0;
  double worms_per_update = 0, sign = 0;
  std::vector<std::uint8_t> initial_state;
  std::int32_t steps = 0;
  unsigned last_id = 0;
  std::vector<unsigned> kink_site_counts, kink_ids;
  std::vector<double> kink_times;
  std::vector<std::uint8_t> kink_states;
  {
    alps::IXDRFileDump dump(prefix.string() + ".xdr");
    dump >> version >> num_kinks >> worms_per_update >> initial_state >> steps;
    std::uint32_t site_count = 0;
    dump >> site_count;
    for (std::uint32_t site = 0; site < site_count; ++site) {
      std::uint32_t count = 0;
      dump >> count;
      kink_site_counts.push_back(count);
      for (std::uint32_t kink = 0; kink < count; ++kink) {
        double time = 0;
        std::uint8_t state = 0;
        unsigned id = 0;
        // Read primitives rather than Kink::load: time_struct wraps endpoints,
        // while the reference must retain the exact serialized double value.
        dump >> time >> state >> id;
        kink_times.push_back(time);
        kink_states.push_back(state);
        kink_ids.push_back(id);
      }
    }
    dump >> sign >> last_id;
  }
  // Existing WRun leaves last_id_ uninitialized before incrementing it.
  // Its arbitrary starting offset propagates to kink IDs; normalize the
  // shared offset modulo unsigned arithmetic, preserving every linkage.
  // remove_jump compares IDs for equality, so zeroing IDs individually
  // would lose meaningful state. Raw bytes remain available alongside HDF5.
  for (auto& id : kink_ids) id -= last_id;
  {
    alps::hdf5::archive archive(prefix.string() + ".h5", "w");
    static_cast<const alps::scheduler::MCRun&>(run).save(archive);
    archive["/characterization/beta"] << run.inverse_temperature();
    archive["/characterization/thermalized"] << run.is_thermalized();
    archive["/characterization/work_done"] << run.work_done();
    archive["/characterization/solver/version"] << version;
    archive["/characterization/solver/num_kinks"] << num_kinks;
    archive["/characterization/solver/worms_per_update"] << worms_per_update;
    archive["/characterization/solver/initial_state"] << initial_state;
    archive["/characterization/solver/steps"] << steps;
    archive["/characterization/solver/sign"] << sign;
    archive["/characterization/solver/kink_site_counts"] << kink_site_counts;
    archive["/characterization/solver/kink_times"] << kink_times;
    archive["/characterization/solver/kink_states"] << kink_states;
    archive["/characterization/solver/kink_ids_relative_to_last_id"] << kink_ids;
  }
}

} // namespace

TEST(WormCharacterization, TemperatureAliasPrecedence) {
  auto p = worm_parameters();
  p["Beta"] = 2;
  p["beta"] = 3;
  p["BETA"] = 4;
  p["T"] = 5;
  p["TEMPERATURE"] = 6;
  p["temperature"] = 8;
  const char* aliases[] = {"Beta", "beta", "BETA", "T", "TEMPERATURE", "temperature"};
  const double expected[] = {2, 3, 4, 0.2, 1. / 6., 0.125};
  for (int i = 0; i < 6; ++i) {
    SCOPED_TRACE(aliases[i]);
    WormProbe run(p);
    EXPECT_DOUBLE_EQ(run.inverse_temperature(), expected[i]);
    p.erase(aliases[i]);
  }
}

TEST(WormCharacterization, TemperatureExpressionsUseParameters) {
  for (const char* name : {"Beta", "beta", "BETA", "T", "TEMPERATURE", "temperature"}) {
    SCOPED_TRACE(name);
    auto p = worm_parameters();
    p.erase("T");
    p["U"] = 1.5;
    p[name] = "4*U";
    WormProbe run(p);
    const bool inverse = std::string(name) == "Beta" || std::string(name) == "beta" ||
                         std::string(name) == "BETA";
    EXPECT_DOUBLE_EQ(run.inverse_temperature(), inverse ? 6 : 1. / 6.);
  }
}

TEST(WormCharacterization, MissingOrUnresolvedTemperatureThrows) {
  auto p = worm_parameters();
  p.erase("T");
  EXPECT_THROW({ WormProbe run(p); }, std::exception);
  p["T"] = "4*unknown_temperature_parameter";
  EXPECT_THROW({ WormProbe run(p); }, std::exception);
}

TEST(WormCharacterization, NegativeTemperatureRejectedAndZeroIsAccepted) {
  auto p = worm_parameters();
  p["T"] = -2;
  EXPECT_THROW({ WormProbe run(p); }, std::out_of_range);
  p["Beta"] = -2;
  EXPECT_THROW({ WormProbe run(p); }, std::out_of_range);
  p["Beta"] = 0;
  WormProbe zero_beta(p);
  EXPECT_DOUBLE_EQ(zero_beta.inverse_temperature(), 0);
  p.erase("Beta");
  p["T"] = 0;
  WormProbe zero_temperature(p);
  EXPECT_TRUE(std::isinf(zero_temperature.inverse_temperature()));
}

TEST(WormCharacterization, ThermalizationRequiresNonnegativeInteger) {
  auto p = worm_parameters();
  p.erase("THERMALIZATION");
  EXPECT_THROW({ WormProbe run(p); }, std::runtime_error);
  p["THERMALIZATION"] = "2*U";
  EXPECT_THROW({ WormProbe run(p); }, boost::bad_lexical_cast);
  p["THERMALIZATION"] = 2.5;
  EXPECT_THROW({ WormProbe run(p); }, boost::bad_lexical_cast);
  p["THERMALIZATION"] = -1;
  EXPECT_THROW({ WormProbe run(p); }, std::out_of_range);
  p["THERMALIZATION"] = 0;
  WormProbe run(p);
  EXPECT_TRUE(run.is_thermalized());
}

TEST(WormCharacterization, ValidationOrderRetainsTemperatureThenEtaThenThermalization) {
  auto p = worm_parameters();
  p["T"] = -2;
  p["eta"] = -1;
  p["THERMALIZATION"] = -1;
  try {
    WormProbe run(p);
    FAIL() << "Expected a negative-temperature error";
  } catch (const std::out_of_range& error) {
    EXPECT_STREQ(error.what(), "negative inverse temperature beta is illegal");
  }
  p["T"] = 2;
  try {
    WormProbe run(p);
    FAIL() << "Expected a negative-eta error";
  } catch (const std::out_of_range& error) {
    EXPECT_STREQ(error.what(), "negative eta is illegal");
  }
  p["eta"] = 1;
  try {
    WormProbe run(p);
    FAIL() << "Expected a negative-thermalization error";
  } catch (const std::out_of_range& error) {
    EXPECT_STREQ(error.what(), "negative thermalization is illegal");
  }
}

TEST(WormCharacterization, SweepsAreRequiredOnlyAfterThermalization) {
  auto p = worm_parameters();
  p.erase("SWEEPS");
  WormProbe run(p);
  EXPECT_DOUBLE_EQ(run.work_done(), 0);
  start_without_log(run);
  run.dostep();
  EXPECT_DOUBLE_EQ(run.work_done(), 0);
  run.dostep();
  EXPECT_TRUE(run.is_thermalized());
  EXPECT_THROW(run.work_done(), std::runtime_error);
}

TEST(WormCharacterization, ProgressStartsAtThermalizationBoundary) {
  auto p = worm_parameters();
  WormProbe run(p);
  start_without_log(run);
  EXPECT_FALSE(run.is_thermalized());
  EXPECT_DOUBLE_EQ(run.work_done(), 0);
  for (int step = 1; step <= 12; ++step) {
    run.dostep();
    EXPECT_EQ(run.is_thermalized(), step >= 2);
    EXPECT_DOUBLE_EQ(run.work_done(), step < 2 ? 0 : (step - 2) / 8.);
  }
  dump_reference("progress", run);
}

TEST(WormCharacterization, ChangedParametersExceptSeedAreRejected) {
  auto p = worm_parameters();
  WormProbe run(p);
  auto changed = p;
  changed["SEED"] = 6789;
  EXPECT_NO_THROW(run.set_parameters(changed));
  for (const char* name : {"SWEEPS", "THERMALIZATION", "T"}) {
    SCOPED_TRACE(name);
    changed = p;
    changed[name] = 17;
    EXPECT_THROW(run.set_parameters(changed), std::runtime_error);
  }
  changed = p;
  changed["T"] = "4/2";
  EXPECT_NO_THROW(run.set_parameters(changed));
}

TEST(WormCharacterization, CheckpointRetainsProgressAndCannotExtendSweeps) {
  auto p = worm_parameters();
  WormProbe original(p);
  start_without_log(original);
  for (int step = 0; step < 5; ++step) original.dostep();
  alps::testing::TemporaryDirectory temporary;
  const auto checkpoint = temporary.path() / "worm.xdr";
  {
    alps::OXDRFileDump dump(checkpoint.string());
    original.save(dump);
  }
  WormProbe resumed(p);
  {
    alps::IXDRFileDump dump(checkpoint.string());
    resumed.load(dump);
  }
  EXPECT_EQ(resumed.is_thermalized(), original.is_thermalized());
  EXPECT_DOUBLE_EQ(resumed.work_done(), original.work_done());
  EXPECT_DOUBLE_EQ(resumed.inverse_temperature(), original.inverse_temperature());
  auto extension = p;
  extension["SWEEPS"] = 16;
  EXPECT_THROW(resumed.set_parameters(extension), std::runtime_error);
  start_without_log(resumed);
  for (int step = 0; step < 5; ++step) resumed.dostep();
  dump_reference("resumed", resumed);
}

TEST(WormCharacterization, ReferenceTrajectoriesForTemperatureAndThermalization) {
  for (int case_index = 0; case_index < 3; ++case_index) {
    auto p = worm_parameters();
    if (case_index == 1) {
      p["Beta"] = "2*U";
      p["T"] = 99;
    } else if (case_index == 2) {
      p["THERMALIZATION"] = 0;
      p["TEMPERATURE"] = "4*U";
      p.erase("T");
    }
    WormProbe run(p);
    start_without_log(run);
    for (int step = 0; step < 18; ++step) run.dostep();
    dump_reference("trajectory_" + std::to_string(case_index), run);
  }
}
