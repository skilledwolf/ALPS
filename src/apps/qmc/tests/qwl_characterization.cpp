// SPDX-License-Identifier: MIT
// Characterizes the existing QWL worker before extracting its settings logic.
#include <gtest/gtest.h>
#include "../qwl/qwl_sse.h"
#include "construction_draws.h"
#include <alps/osiris/xdrdump.h>
#include <alps/testing/temporary_directory.hpp>
#include <cstdlib>
#include <filesystem>
#include <tuple>

namespace {
alps::Parameters qwl_parameters() {
  alps::Parameters p;
  p["LATTICE"] = "chain lattice";
  p["MODEL"] = "spin";
  p["local_S"] = "1/2";
  p["L"] = 4;
  p["J"] = 1;
  p["SEED"] = 731;
  p["MEASURE_MAGNETIC_PROPERTIES"] = 0;
  return p;
}

alps::ProcessList local_process() { return alps::ProcessList(1); }

struct SavedHistogram {
  std::vector<double> values;
  unsigned int left = 0;
  void read(alps::IDump& dump) { dump >> values >> left; }
  bool operator==(const SavedHistogram& rhs) const {
    return values == rhs.values && left == rhs.left;
  }
  void write(alps::hdf5::archive& archive, const std::string& path) const {
    archive[path + "/values"] << values;
    archive[path + "/left"] << left;
  }
};

// Decode the established simulation checkpoint format rather than exposing
// private members. Identity operators have uninitialized bond fields in the
// original implementation; their unused bonds are normalized to zero here.
// Active operator bonds, all histogram values and the spin state are retained.
struct QwlState {
  SavedHistogram g, histo, histoup, nmeasurements, transition_prob;
  SavedHistogram uniform_structure_factor, staggered_structure_factor;
  unsigned int sweeps, cutoff, non0, logf_step, block_sweeps, block_sweeps_total;
  unsigned int norder, doing_multicanonical, upwalker, thistime, all_done;
  double logf, minimum_histogram, flatness_treshold;
  std::vector<unsigned int> spins, ops, bonds;

  void read(alps::IDump& dump, bool magnetic) {
    g.read(dump);
    histo.read(dump);
    histoup.read(dump);
    nmeasurements.read(dump);
    transition_prob.read(dump);
    if (magnetic) {
      uniform_structure_factor.read(dump);
      staggered_structure_factor.read(dump);
    }
    dump >> sweeps >> cutoff >> non0 >> logf >> logf_step;
    dump >> block_sweeps >> block_sweeps_total >> norder;
    dump >> minimum_histogram >> flatness_treshold;
    dump >> doing_multicanonical >> upwalker >> thistime >> all_done >> spins;
    ops.resize(cutoff);
    bonds.resize(cutoff);
    for (unsigned int i = 0; i < cutoff; ++i) {
      dump >> ops[i] >> bonds[i];
      if (!ops[i]) bonds[i] = 0;
    }
  }

  bool operator==(const QwlState& rhs) const {
    return std::tie(g, histo, histoup, nmeasurements, transition_prob,
                    uniform_structure_factor, staggered_structure_factor,
                    sweeps, cutoff, non0, logf_step, block_sweeps,
                    block_sweeps_total, norder, doing_multicanonical, upwalker,
                    thistime, all_done, logf, minimum_histogram,
                    flatness_treshold, spins, ops, bonds) ==
           std::tie(rhs.g, rhs.histo, rhs.histoup, rhs.nmeasurements,
                    rhs.transition_prob, rhs.uniform_structure_factor,
                    rhs.staggered_structure_factor, rhs.sweeps, rhs.cutoff,
                    rhs.non0, rhs.logf_step, rhs.block_sweeps,
                    rhs.block_sweeps_total, rhs.norder, rhs.doing_multicanonical,
                    rhs.upwalker, rhs.thistime, rhs.all_done, rhs.logf,
                    rhs.minimum_histogram, rhs.flatness_treshold,
                    rhs.spins, rhs.ops, rhs.bonds);
  }

  void write(alps::hdf5::archive& archive) const {
    g.write(archive, "/qwl/g");
    histo.write(archive, "/qwl/histo");
    histoup.write(archive, "/qwl/histoup");
    nmeasurements.write(archive, "/qwl/nmeasurements");
    transition_prob.write(archive, "/qwl/transition_prob");
    if (!uniform_structure_factor.values.empty()) {
      uniform_structure_factor.write(archive, "/qwl/uniform_structure_factor");
      staggered_structure_factor.write(archive, "/qwl/staggered_structure_factor");
    }
#define WRITE_QWL_FIELD(field) archive["/qwl/" #field] << field
    WRITE_QWL_FIELD(sweeps);
    WRITE_QWL_FIELD(cutoff);
    WRITE_QWL_FIELD(non0);
    WRITE_QWL_FIELD(logf);
    WRITE_QWL_FIELD(logf_step);
    WRITE_QWL_FIELD(block_sweeps);
    WRITE_QWL_FIELD(block_sweeps_total);
    WRITE_QWL_FIELD(norder);
    WRITE_QWL_FIELD(minimum_histogram);
    WRITE_QWL_FIELD(flatness_treshold);
    WRITE_QWL_FIELD(doing_multicanonical);
    WRITE_QWL_FIELD(upwalker);
    WRITE_QWL_FIELD(thistime);
    WRITE_QWL_FIELD(all_done);
    WRITE_QWL_FIELD(spins);
    WRITE_QWL_FIELD(ops);
    WRITE_QWL_FIELD(bonds);
#undef WRITE_QWL_FIELD
  }
};

QwlState state_of(const QWL_SSE_Simulation& run, bool magnetic = false) {
  alps::testing::TemporaryDirectory directory;
  const boost::filesystem::path path((directory.path() / "qwl.dump").string());
  {
    alps::OXDRFileDump dump(path);
    run.save(dump);
  }
  alps::IXDRFileDump dump(path);
  QwlState result;
  result.read(dump, magnetic);
  return result;
}

void maybe_dump(const QWL_SSE_Simulation& run, const std::string& name,
                bool magnetic = false) {
  const char* root = std::getenv("ALPS_QMC_DUMP_DIR");
  if (!root || !*root) return;
  std::filesystem::create_directories(root);
  const auto path = std::filesystem::path(root) / ("qwl_" + name + ".h5");
  alps::hdf5::archive archive(path.string(), "w");
  run.alps::scheduler::MCRun::save(archive);
  state_of(run, magnetic).write(archive);
}

alps::Parameters one_bin_schedule() {
  auto p = qwl_parameters();
  p["EXPANSION_ORDER_MAXIMUM"] = 0;
  p["CUTOFF"] = 4;
  p["USE_ZHOU_BHATT_METHOD"] = 0;
  p["BLOCK_SWEEPS"] = 4;
  p["INITIAL_INCREASE_FACTOR"] = std::exp(1.);
  p["NUMBER_OF_WANG_LANDAU_STEPS"] = 1;
  p["SWEEPS"] = 3;
  return p;
}
} // namespace

TEST(QwlCharacterization, DefaultsUseZhouBhattAndSixteenWeightSteps) {
  QWL_SSE_Simulation run(local_process(), qwl_parameters(), 0);
  const auto state = state_of(run);
  EXPECT_EQ(state.g.left, 0u);
  EXPECT_EQ(state.g.values.size(), 501u);
  EXPECT_EQ(state.cutoff, 500u);
  EXPECT_EQ(state.block_sweeps_total, 501u);
  EXPECT_DOUBLE_EQ(state.logf, 1.);
  EXPECT_DOUBLE_EQ(state.minimum_histogram, 1.);
  EXPECT_DOUBLE_EQ(state.flatness_treshold, 1E10);
  EXPECT_EQ(state.logf_step, 1u);
  EXPECT_TRUE(run.get_measurements().has("Coefficients 16"));
  EXPECT_FALSE(run.is_thermalized());
  EXPECT_DOUBLE_EQ(run.work_done(), 0.);
  maybe_dump(run, "default_initial");
}

TEST(QwlCharacterization, OriginalScheduleHasDifferentDefaults) {
  auto p = qwl_parameters();
  p["USE_ZHOU_BHATT_METHOD"] = 0;
  QWL_SSE_Simulation run(local_process(), p, 0);
  const auto state = state_of(run);
  EXPECT_EQ(state.block_sweeps_total, 10000u);
  EXPECT_DOUBLE_EQ(state.logf, std::log(std::exp(500 * std::log(4.) / 10000)));
  EXPECT_DOUBLE_EQ(state.minimum_histogram, 0.);
  EXPECT_DOUBLE_EQ(state.flatness_treshold, .2);
  maybe_dump(run, "original_initial");
}

TEST(QwlCharacterization, ExplicitScheduleAndOrderBounds) {
  auto p = qwl_parameters();
  p["EXPANSION_ORDER_MINIMUM"] = 2;
  p["EXPANSION_ORDER_MAXIMUM"] = 4;
  p["CUTOFF"] = 8;
  p["NUMBER_OF_WANG_LANDAU_STEPS"] = 3;
  p["INITIAL_MODIFICATION_FACTOR"] = 4.;
  p["FLATNESS_TRESHOLD"] = .125;
  p["BLOCK_SWEEPS"] = 97; // ignored by the Zhou-Bhatt schedule
  QWL_SSE_Simulation run(local_process(), p, 0);
  const auto state = state_of(run);
  EXPECT_EQ(state.g.left, 2u);
  EXPECT_EQ(state.g.values.size(), 3u);
  EXPECT_EQ(state.cutoff, 8u);
  EXPECT_EQ(state.block_sweeps_total, 3u);
  EXPECT_DOUBLE_EQ(state.logf, std::log(4.));
  EXPECT_DOUBLE_EQ(state.minimum_histogram, 1. / std::log(4.));
  EXPECT_DOUBLE_EQ(state.flatness_treshold, .125);
  EXPECT_TRUE(run.get_measurements().has("Coefficients 3"));
  // The original worker crashes while stepping at a positive minimum order;
  // this constructor characterization deliberately does not step that case.
  maybe_dump(run, "positive_minimum_initial");
}

TEST(QwlCharacterization, RejectsInvalidOrderBounds) {
  auto p = qwl_parameters();
  p["EXPANSION_ORDER_MINIMUM"] = -1;
  EXPECT_THROW(QWL_SSE_Simulation(local_process(), p, 0), std::runtime_error);
  p["EXPANSION_ORDER_MINIMUM"] = 0;
  p["EXPANSION_ORDER_MAXIMUM"] = 5;
  p["CUTOFF"] = 4;
  EXPECT_THROW(QWL_SSE_Simulation(local_process(), p, 0), std::runtime_error);
}

TEST(QwlCharacterization, ScheduleUsesPlainNumbers) {
  auto p = qwl_parameters();
  p["CUTOFF"] = "2*U";
  p["U"] = 4;
  EXPECT_THROW(QWL_SSE_Simulation(local_process(), p, 0), boost::bad_lexical_cast);
}

TEST(QwlCharacterization, TemperatureAndThermalizationAreIgnored) {
  auto p = qwl_parameters();
  p["CUTOFF"] = 12;
  QWL_SSE_Simulation reference(local_process(), p, 0);
  for (int i = 0; i < 10; ++i) reference.dostep();
  const auto expected = state_of(reference);
  for (const char* key : {"Beta", "beta", "BETA", "T", "TEMPERATURE",
                          "temperature", "INVERSE_TEMPERATURE", "THERMALIZATION"}) {
    SCOPED_TRACE(key);
    auto changed = p;
    changed[key] = "invalid temperature";
    QWL_SSE_Simulation run(local_process(), changed, 0);
    for (int i = 0; i < 10; ++i) run.dostep();
    EXPECT_TRUE(state_of(run) == expected);
    EXPECT_TRUE(run.is_thermalized());
  }
  // QWL hardcodes beta=1; changing temperature leaves updates and RNG state
  // unchanged. Save a trajectory with nonzero expansion order for A/B checks.
  maybe_dump(reference, "temperature_ignored_after_ten");
}

TEST(QwlCharacterization, ProgressStartsAfterFinalWeightStep) {
  QWL_SSE_Simulation run(local_process(), one_bin_schedule(), 0);
  for (int i = 0; i < 4; ++i) {
    run.dostep();
    EXPECT_TRUE(run.is_thermalized());
    EXPECT_DOUBLE_EQ(run.work_done(), 0.);
  }
  auto state = state_of(run);
  EXPECT_EQ(state.doing_multicanonical, 1u);
  EXPECT_EQ(state.block_sweeps_total, 3u);
  EXPECT_EQ(state.sweeps, 4u);
  maybe_dump(run, "multicanonical_start");
  run.dostep();
  EXPECT_DOUBLE_EQ(run.work_done(), 1. / 3.);
  run.dostep();
  EXPECT_DOUBLE_EQ(run.work_done(), 2. / 3.);
  run.dostep();
  EXPECT_DOUBLE_EQ(run.work_done(), 1.);
  EXPECT_EQ(state_of(run).all_done, 1u);
  maybe_dump(run, "multicanonical_done");
}

TEST(QwlCharacterization, OnlySweepsCanChange) {
  QWL_SSE_Simulation run(local_process(), one_bin_schedule(), 0);
  EXPECT_TRUE(run.change_parameter("SWEEPS", alps::StringValue(12)));
  EXPECT_EQ(state_of(run).block_sweeps_total, 12u);
  for (const char* key : {"SEED", "THERMALIZATION", "T", "CUTOFF",
                          "NUMBER_OF_WANG_LANDAU_STEPS", "BLOCK_SWEEPS"}) {
    SCOPED_TRACE(key);
    EXPECT_FALSE(run.change_parameter(key, alps::StringValue(8)));
  }
  maybe_dump(run, "changed_sweeps");
}

TEST(QwlCharacterization, CheckpointRestoresScheduleAndSpinState) {
  const auto p = one_bin_schedule();
  QWL_SSE_Simulation original(local_process(), p, 0);
  for (int i = 0; i < 5; ++i) original.dostep();
  const auto expected = state_of(original);
  alps::testing::TemporaryDirectory directory;
  const boost::filesystem::path path((directory.path() / "qwl.dump").string());
  {
    alps::OXDRFileDump dump(path);
    original.save(dump);
  }
  auto changed = p;
  changed["BLOCK_SWEEPS"] = 9;
  changed["INITIAL_INCREASE_FACTOR"] = 4.;
  changed["FLATNESS_TRESHOLD"] = .7;
  QWL_SSE_Simulation restored(local_process(), changed, 0);
  {
    alps::IXDRFileDump dump(path);
    restored.load(dump);
  }
  EXPECT_TRUE(state_of(restored) == expected);
  EXPECT_DOUBLE_EQ(restored.work_done(), original.work_done());
  maybe_dump(restored, "checkpoint_restored");
  restored.dostep();
  maybe_dump(restored, "checkpoint_continued");
}

class QwlProbe : public QWL_SSE_Simulation {
public:
  explicit QwlProbe(const alps::Parameters& p) : QWL_SSE_Simulation(local_process(), p, 0) {}
  double draw() { return random_01(); }
};

// The start configuration is drawn while the run is constructed: one number
// per site of the four-site chain.
TEST(QwlCharacterization, DrawsTheStartConfigurationDuringConstruction) {
  QwlProbe run(qwl_parameters());
  EXPECT_EQ(construction_draws(qwl_parameters(), [&] { return run.draw(); }), 4);
}
