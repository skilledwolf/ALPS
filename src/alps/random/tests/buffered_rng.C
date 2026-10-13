// SPDX-License-Identifier: MIT
// Restoring or reseeding a buffered generator must continue from the new
// engine state, not from numbers buffered before it was replaced.
#include <alps/random/buffered_rng.h>
#include <alps/osiris/comm.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/testing/temporary_directory.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

namespace {

typedef alps::buffered_rng<boost::mt19937> generator;

std::vector<uint32_t> draw(alps::buffered_rng_base& rng, int n) {
  std::vector<uint32_t> out;
  for (int i = 0; i < n; ++i)
    out.push_back(rng());
  return out;
}

std::string state_after(int draws) {
  generator rng;
  rng.seed(42u);
  draw(rng, draws);
  std::ostringstream state;
  rng.write(state);
  return state.str();
}

TEST(BufferedRng, ReadContinuesFromTheReadStateOnly) {
  const std::string state = state_after(50000);
  // Two generators with different histories read the same state.
  generator fresh, used;
  used.seed(42u);
  draw(used, 10);
  std::istringstream a(state), b(state);
  fresh.read(a);
  used.read(b);
  EXPECT_EQ(draw(fresh, 30000), draw(used, 30000));
}

TEST(BufferedRng, SeedingDiscardsBufferedNumbers) {
  generator fresh, used;
  used.seed(7u);
  draw(used, 10);
  fresh.seed(42u);
  used.seed(42u);
  EXPECT_EQ(draw(fresh, 30000), draw(used, 30000));
}

TEST(BufferedRng, WriteAllRestoresTheBufferToo) {
  generator original, restored;
  original.seed(42u);
  draw(original, 10);
  std::stringstream state;
  original.write_all(state);
  restored.seed(7u);
  draw(restored, 3);
  restored.read_all(state);
  EXPECT_EQ(draw(restored, 30000), draw(original, 30000));
}

std::vector<double> drawn;

// A run that draws an initial configuration when it is built, as
// simulations do, and one number per step.
class drawing_run : public alps::scheduler::MCRun {
public:
  drawing_run(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : alps::scheduler::MCRun(w, p, n) {
    for (int i = 0; i < 10; ++i)
      random_01();
  }
  void dostep() override { drawn.push_back(random_01()); }
  double work_done() const override { return 0.; }
  bool is_thermalized() const override { return true; }
  void save(alps::ODump&) const override {}
  void load(alps::IDump&) override {}
};

// A run resumed from a checkpoint is built with the same SEED as the
// original run. It must not continue with the numbers its own seed put in
// the buffer, which the original run already used in its first steps.
TEST(BufferedRng, ResumedRunDoesNotRepeatTheFirstSteps) {
  alps::Parameters p;
  std::istringstream input("SEED=42");
  p.parse(input);
  const alps::ProcessList here(1, alps::local_process());
  alps::testing::TemporaryDirectory directory;
  const auto checkpoint = (directory.path() / "run").string();

  drawn.clear();
  drawing_run original(here, p, 0);
  original.start_worker();
  for (int i = 0; i < 20000; ++i)
    original.dostep();
  original.save_to_file(checkpoint, checkpoint + ".h5");
  const std::vector<double> first(drawn.begin(), drawn.begin() + 1000);

  drawn.clear();
  drawing_run resumed(here, p, 0);
  resumed.load_from_file(checkpoint, checkpoint + ".h5");
  for (int i = 0; i < 1000; ++i)
    resumed.dostep();
  EXPECT_NE(drawn, first);
  int repeated = 0;
  for (std::size_t i = 0; i < first.size(); ++i)
    repeated += drawn[i] == first[i];
  EXPECT_LT(repeated, 5);
}

} // namespace
