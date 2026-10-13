// SPDX-License-Identifier: MIT
// These resolvers require Parameters only, with no lattice, model or worker.
#include "../settings.h"
#include "../worms/settings.h"

#ifdef ALPS_SCHEDULER_WORKER_H
#error "Settings resolvers must not include scheduler workers"
#endif

#include <gtest/gtest.h>
#include <boost/lexical_cast.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

TEST(QmcSettings, TemperatureAliasPrecedence) {
  alps::Parameters p;
  p["Beta"] = 2;
  p["beta"] = 3;
  p["BETA"] = 4;
  p["T"] = 5;
  p["TEMPERATURE"] = 6;
  p["temperature"] = 8;
  const char* names[] = {"Beta", "beta", "BETA", "T", "TEMPERATURE", "temperature"};
  const double values[] = {2, 3, 4, 0.2, 1. / 6., 0.125};
  for (int i = 0; i < 6; ++i) {
    SCOPED_TRACE(names[i]);
    EXPECT_DOUBLE_EQ(qmc::resolve_inverse_temperature(p), values[i]);
    p.erase(names[i]);
  }
}

TEST(QmcSettings, TemperatureExpressionsUseParameters) {
  for (const char* name : {"Beta", "beta", "BETA", "T", "TEMPERATURE", "temperature"}) {
    SCOPED_TRACE(name);
    alps::Parameters p;
    p["U"] = 1.5;
    p[name] = "4*U";
    const bool inverse = std::string(name) == "Beta" || std::string(name) == "beta" ||
                         std::string(name) == "BETA";
    EXPECT_DOUBLE_EQ(qmc::resolve_inverse_temperature(p), inverse ? 6 : 1. / 6.);
  }
}

TEST(QmcSettings, MissingAndUnresolvedTemperatureKeepEvaluationErrors) {
  alps::Parameters p;
  EXPECT_THROW(qmc::resolve_inverse_temperature(p), std::exception);
  p["T"] = "4*unknown_temperature_parameter";
  EXPECT_THROW(qmc::resolve_inverse_temperature(p), std::exception);
  p["T"] = 2;
  p["Beta"] = "unknown_temperature_parameter";
  EXPECT_THROW(qmc::resolve_inverse_temperature(p), std::exception);
}

TEST(QmcSettings, NegativeAndZeroValuesRemainForConstructorValidation) {
  alps::Parameters p;
  p["T"] = -2;
  EXPECT_DOUBLE_EQ(qmc::resolve_inverse_temperature(p), -0.5);
  p["Beta"] = -3;
  EXPECT_DOUBLE_EQ(qmc::resolve_inverse_temperature(p), -3);
  p["Beta"] = 0;
  EXPECT_DOUBLE_EQ(qmc::resolve_inverse_temperature(p), 0);
  p.erase("Beta");
  p["T"] = 0;
  EXPECT_TRUE(std::isinf(qmc::resolve_inverse_temperature(p)));
}

TEST(WormSettings, ThermalizationUsesRequiredSignedInteger) {
  alps::Parameters p;
  EXPECT_THROW(worm::resolve_thermalization(p), std::runtime_error);
  p["THERMALIZATION"] = 0;
  EXPECT_EQ(worm::resolve_thermalization(p), 0);
  p["THERMALIZATION"] = 13;
  EXPECT_EQ(worm::resolve_thermalization(p), 13);
  // WRun checks sign after eta, so the parsing seam must retain this value.
  p["THERMALIZATION"] = -1;
  EXPECT_EQ(worm::resolve_thermalization(p), -1);
}

TEST(WormSettings, ThermalizationRejectsExpressionsAndFractions) {
  alps::Parameters p;
  p["U"] = 2;
  p["THERMALIZATION"] = "2*U";
  EXPECT_THROW(worm::resolve_thermalization(p), boost::bad_lexical_cast);
  p["THERMALIZATION"] = 2.5;
  EXPECT_THROW(worm::resolve_thermalization(p), boost::bad_lexical_cast);
}

TEST(WormSettings, ThermalizationKeepsInt32Range) {
  alps::Parameters p;
  p["THERMALIZATION"] = std::numeric_limits<std::int32_t>::max();
  EXPECT_EQ(worm::resolve_thermalization(p), std::numeric_limits<std::int32_t>::max());
  p["THERMALIZATION"] = "2147483648";
  EXPECT_THROW(worm::resolve_thermalization(p), boost::bad_lexical_cast);
  p["THERMALIZATION"] = "-2147483649";
  EXPECT_THROW(worm::resolve_thermalization(p), boost::bad_lexical_cast);
}
