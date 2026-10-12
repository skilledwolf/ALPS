#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
// Numerical regression tests for every FREQUENCY_GRID of the maxent
// tool. ContiParameters allocates t_array_ with nfreq_+1 knots; leaving the
// last knot unwritten corrupts the final bin width to ~(OMEGA_MIN-OMEGA_MAX),
// which blows the default-model normalisation up to +-inf and makes every
// linear-grid run return an all-NaN spectrum.
//
// The test synthesises a clean fermionic G(tau) from a known two-Gaussian
// spectral function (peaks at omega = +-1.5), runs MaxEntSimulation
// in-process on each grid, and asserts the recovered spectrum is
// finite and satisfies the sum rule int A(omega) domega ~= 1.
#include "maxent.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/ngs/params.hpp>
#include <boost/filesystem.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
bool all_finite(const std::vector<double>& v) {
  if (v.empty()) return false;
  for (std::size_t i = 0; i < v.size(); ++i)
    if (!std::isfinite(v[i])) return false;
  return true;
}
double trapz(const std::vector<double>& x, const std::vector<double>& y) {
  double s = 0.0;
  for (std::size_t i = 0; i + 1 < y.size(); ++i)
    s += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]);
  return s;
}
}

class MaxEntFrequencyGrid : public ::testing::TestWithParam<std::string> {};

TEST_P(MaxEntFrequencyGrid, FiniteNormalizedSpectrumAndStopCallbacks) {

  const double beta = 10.0;
  const int n_tau = 40;
  const int nom = 401;
  const double om_lo = -8.0, om_hi = 8.0;

  std::vector<double> tau(n_tau), om(nom), rho(nom);
  for (int i = 0; i < n_tau; ++i) tau[i] = beta * double(i) / double(n_tau - 1);
  for (int j = 0; j < nom; ++j) om[j] = om_lo + (om_hi - om_lo) * double(j) / double(nom - 1);
  const double dom = om[1] - om[0];
  for (int j = 0; j < nom; ++j)
    rho[j] = std::exp(-std::pow(om[j] - 1.5, 2) / 0.4)
           + std::exp(-std::pow(om[j] + 1.5, 2) / 0.4);
  double rho_int = 0.0;
  for (int j = 0; j + 1 < nom; ++j) rho_int += 0.5 * (rho[j] + rho[j + 1]) * dom;
  for (int j = 0; j < nom; ++j) rho[j] /= rho_int;

  std::vector<double> Gin(n_tau);
  for (int i = 0; i < n_tau; ++i) {
    double acc = 0.0;
    for (int j = 0; j < nom; ++j)
      acc += std::exp(-tau[i] * om[j]) / (1.0 + std::exp(-beta * om[j])) * rho[j];
    Gin[i] = -acc * dom;
  }
  std::vector<double> errors(n_tau, 1e-3);

  alps::params p;
  p["BETA"] = beta;
  p["NDAT"] = n_tau;
  p["NFREQ"] = GetParam() == "log" ? 401 : 400;
  p["NORM"] = 1.0;
  p["KERNEL"] = std::string("fermionic");
  p["DATASPACE"] = std::string("time");
  p["OMEGA_MIN"] = om_lo;
  p["OMEGA_MAX"] = om_hi;
  p["N_ALPHA"] = 12;
  p["ALPHA_MIN"] = 0.1;
  p["ALPHA_MAX"] = 100.0;
  p["DEFAULT_MODEL"] = std::string("flat");
  p["FREQUENCY_GRID"] = GetParam();
  p["MAX_IT"] = 40;
  p["PARTICLE_HOLE_SYMMETRY"] = 0;
  p["TEXT_OUTPUT"] = 0;
  p["SELF"] = true;
  p["DATA_IN_HDF5"] = 1;
  p["MAX_TIME"] = 600;
  p["VERBOSE"] = 0;

  alps::testing::TemporaryDirectory temporary;
  struct restore_directory {
    std::filesystem::path previous = std::filesystem::current_path();
    ~restore_directory() { std::filesystem::current_path(previous); }
  } restore;
  std::filesystem::current_path(temporary.path());
  const std::string base = (temporary.path() / "maxent").string();
  const std::string in_h5 = base + ".h5";
  const std::string out_h5 = base + ".custom.hdf5";
  p["DATA"] = in_h5;
  {
    alps::hdf5::archive ar(in_h5, "w");
    ar << alps::make_pvp("/Data", Gin);
    ar << alps::make_pvp("/Error", errors);
  }
  std::remove(out_h5.c_str());

  {
    MaxEntSimulation sim(p, out_h5);
    int callbacks = 0;
    ASSERT_TRUE((!sim.run([&]() { ++callbacks; return true; }))) << "initial stop request ignored";
    ASSERT_TRUE((callbacks == 1)) << "initial stop callback count changed";
    ASSERT_TRUE((!boost::filesystem::exists(out_h5))) << "stopped run created scientific output";

    bool propagated = false;
    try {
      sim.run([]() -> bool { throw std::logic_error("stop callback failure"); });
    } catch (std::logic_error const& error) {
      propagated = std::string(error.what()) == "stop callback failure";
    }
    ASSERT_TRUE((propagated)) << "stop callback exception not propagated";

    callbacks = 0;
    ASSERT_TRUE((sim.run([&]() { ++callbacks; return false; }))) << "complete run reported stopped";
    ASSERT_TRUE((callbacks == 2)) << "callback must run before and after the alpha sweep";

    callbacks = 0;
    ASSERT_TRUE((!sim.run([&]() { ++callbacks; return true; }))) << "completed run ignored stop callback";
    ASSERT_TRUE((callbacks == 1)) << "completed run must check the callback once";
  }

  // Construction and execution must not leak auxiliary files when text is off.
  EXPECT_EQ(std::distance(std::filesystem::directory_iterator(temporary.path()),
                          std::filesystem::directory_iterator{}), 2);
  std::vector<double> w, Aavg, Amax, Achi, widths, self;
  {
    alps::hdf5::archive ar(out_h5, "r");
    ar >> alps::make_pvp("/spectrum/omega", w);
    ar["/spectrum/bin_widths"] >> widths;
    ar["/spectrum/self/average"] >> self;
    ar >> alps::make_pvp("/spectrum/average", Aavg);
    ar >> alps::make_pvp("/spectrum/maximum", Amax);
    ar >> alps::make_pvp("/spectrum/chi", Achi);
  }

  ASSERT_EQ(widths.size(), w.size());
  ASSERT_EQ(self.size(), Aavg.size());
  for (std::size_t i=0; i<widths.size(); ++i) {
    EXPECT_GT(widths[i], 0.);
    EXPECT_NEAR(self[i], -std::acos(-1.)*Aavg[i], 1.e-14);
  }
  ASSERT_TRUE((!w.empty())) << "omega grid empty";
  ASSERT_TRUE((all_finite(w))) << "omega grid not finite";
  ASSERT_TRUE((all_finite(Aavg))) << "A_average not all-finite (NaN regression)";
  ASSERT_TRUE((all_finite(Amax))) << "A_maximum not all-finite (NaN regression)";
  ASSERT_TRUE((all_finite(Achi))) << "A_chi2 not all-finite (NaN regression)";

  if (all_finite(w) && all_finite(Aavg)) {
    const double sumrule = trapz(w, Aavg);
    std::fprintf(stderr, "sumrule = %.4f\n", sumrule);
    ASSERT_TRUE((std::abs(sumrule - 1.0) < 0.05)) << "sum rule violated (expected ~1.0)";
  }

}

INSTANTIATE_TEST_SUITE_P(AllGrids, MaxEntFrequencyGrid,
    ::testing::Values("linear", "Lorentzian", "half Lorentzian", "quadratic", "log"));

namespace {
struct GridParameters : ContiParameters {
  explicit GridParameters(const alps::params& p) : ContiParameters(p) {}
  const vector_type& knots() const { return t_array_; }
};
alps::params grid_parameters(const std::string& grid, int nfreq) {
  alps::params p;
  p["BETA"] = 2.;
  p["NDAT"] = 4;
  p["NFREQ"] = nfreq;
  p["NORM"] = 1.;
  p["OMEGA_MIN"] = -8.;
  p["OMEGA_MAX"] = 8.;
  p["DEFAULT_MODEL"] = std::string("flat");
  p["FREQUENCY_GRID"] = grid;
  for (int i=0; i<4; ++i) {
    p["X_" + std::to_string(i)] = -0.5;
    p["SIGMA_" + std::to_string(i)] = 0.01;
  }
  return p;
}
}

TEST(MaxEnt, EveryGridFillsAllIntegrationKnots)
{
  for (const std::string grid : {"linear", "Lorentzian", "half Lorentzian", "quadratic", "log"}) {
    for (const int nfreq : {1, 2, 3, 4, 31, 32, 2001}) {
      SCOPED_TRACE(grid + "/" + std::to_string(nfreq));
      GridParameters parameters(grid_parameters(grid, nfreq));
      const auto& knots = parameters.knots();
      ASSERT_EQ(knots.size(), nfreq + 1);
      EXPECT_DOUBLE_EQ(knots[0], 0.);
      EXPECT_DOUBLE_EQ(knots[nfreq], 1.);
      double width_sum = 0.;
      for (int i=0; i<nfreq; ++i) {
        ASSERT_TRUE(std::isfinite(knots[i+1]));
        ASSERT_GT(knots[i+1], knots[i]);
        EXPECT_LE(knots[i+1], 1.);
        width_sum += knots[i+1] - knots[i];
        if (grid != "half Lorentzian")
          EXPECT_NEAR(knots[i] + knots[nfreq-i], 1., 2e-13);
      }
      EXPECT_NEAR(width_sum, 1., 2e-13);
    }
  }
}

TEST(MaxEnt, GridParametersRejectUndefinedDomains)
{
  for (int count : {-1, 0})
    EXPECT_THROW(GridParameters(grid_parameters("linear", count)), std::invalid_argument);
  for (const std::string grid : {"Lorentzian", "half Lorentzian", "log", "quadratic"}) {
    auto p = grid_parameters(grid, 8);
    const std::string parameter = grid == "log" ? "LOG_MIN" : grid == "quadratic" ? "SPREAD" : "CUT";
    p[parameter] = 0.;
    EXPECT_THROW(GridParameters{p}, std::invalid_argument);
    p[parameter] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(GridParameters{p}, std::invalid_argument);
  }
}
