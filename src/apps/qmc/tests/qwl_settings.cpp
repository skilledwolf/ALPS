// SPDX-License-Identifier: MIT
// The settings resolver is exercised without constructing a lattice or worker.
#include <gtest/gtest.h>
#include "../qwl/settings.h"

#ifdef ALPS_SCHEDULER_WORKER_H
#error "Settings resolvers must not include scheduler workers"
#endif
#include <limits>

TEST(QwlSettings, ZhouBhattDefaultsDependOnOrderWindow) {
  const auto schedule=qwl::resolve_wang_landau_schedule(alps::Parameters(),8,5,4);
  EXPECT_EQ(schedule.use_zhou_bhatt,1u);
  EXPECT_EQ(schedule.block_sweeps_total,5u);
  EXPECT_DOUBLE_EQ(schedule.logf,1.);
  EXPECT_DOUBLE_EQ(schedule.minimum_histogram,1.);
  EXPECT_DOUBLE_EQ(schedule.flatness_treshold,1E10);
}

TEST(QwlSettings, OriginalDefaultsDependOnCutoffAndSiteCount) {
  alps::Parameters p;
  p["USE_ZHOU_BHATT_METHOD"]=0;
  const auto schedule=qwl::resolve_wang_landau_schedule(p,500,501,4);
  EXPECT_EQ(schedule.use_zhou_bhatt,0u);
  EXPECT_EQ(schedule.block_sweeps_total,10000u);
  EXPECT_DOUBLE_EQ(schedule.logf,std::log(std::exp(500*std::log(4.)/10000)));
  EXPECT_DOUBLE_EQ(schedule.minimum_histogram,0.);
  EXPECT_DOUBLE_EQ(schedule.flatness_treshold,.2);
}

TEST(QwlSettings, ZhouBhattIgnoresOriginalScheduleParameters) {
  alps::Parameters p;
  p["INITIAL_MODIFICATION_FACTOR"]=4.;
  p["FLATNESS_TRESHOLD"]="0.125";
  p["BLOCK_SWEEPS"]="invalid";
  p["INITIAL_INCREASE_FACTOR"]="invalid";
  const auto schedule=qwl::resolve_wang_landau_schedule(p,8,3,4);
  EXPECT_EQ(schedule.block_sweeps_total,3u);
  EXPECT_DOUBLE_EQ(schedule.logf,std::log(4.));
  EXPECT_DOUBLE_EQ(schedule.minimum_histogram,1./std::log(4.));
  EXPECT_DOUBLE_EQ(schedule.flatness_treshold,.125);
}

TEST(QwlSettings, OriginalScheduleIgnoresModificationFactor) {
  alps::Parameters p;
  p["USE_ZHOU_BHATT_METHOD"]=0;
  p["BLOCK_SWEEPS"]=37;
  p["INITIAL_INCREASE_FACTOR"]=2.;
  p["INITIAL_MODIFICATION_FACTOR"]="invalid";
  p["FLATNESS_TRESHOLD"]=.75;
  const auto schedule=qwl::resolve_wang_landau_schedule(p,8,3,4);
  EXPECT_EQ(schedule.block_sweeps_total,37u);
  EXPECT_DOUBLE_EQ(schedule.logf,std::log(2.));
  EXPECT_DOUBLE_EQ(schedule.minimum_histogram,0.);
  EXPECT_DOUBLE_EQ(schedule.flatness_treshold,.75);
}

TEST(QwlSettings, MethodFlagUsesUnsignedNumericConversion) {
  alps::Parameters p;
  p["USE_ZHOU_BHATT_METHOD"]=2;
  EXPECT_EQ(qwl::resolve_wang_landau_schedule(p,8,3,4).use_zhou_bhatt,2u);
  p["USE_ZHOU_BHATT_METHOD"]=-1;
  EXPECT_EQ(qwl::resolve_wang_landau_schedule(p,8,3,4).use_zhou_bhatt,
            std::numeric_limits<boost::uint32_t>::max());
  p["USE_ZHOU_BHATT_METHOD"]="true";
  EXPECT_THROW(qwl::resolve_wang_landau_schedule(p,8,3,4),boost::bad_lexical_cast);
}

TEST(QwlSettings, BlockSweepsUsesUnsignedNumericConversion) {
  alps::Parameters p;
  p["USE_ZHOU_BHATT_METHOD"]=0;
  p["BLOCK_SWEEPS"]=-1;
  p["INITIAL_INCREASE_FACTOR"]=2.;
  EXPECT_EQ(qwl::resolve_wang_landau_schedule(p,8,3,4).block_sweeps_total,
            std::numeric_limits<unsigned int>::max());
  p["BLOCK_SWEEPS"]="2*U";
  p["U"]=4;
  EXPECT_THROW(qwl::resolve_wang_landau_schedule(p,8,3,4),boost::bad_lexical_cast);
}

TEST(QwlSettings, RetainsHistoricFlatnessParameterSpelling) {
  alps::Parameters p;
  p["FLATNESS_THRESHOLD"]="invalid";
  EXPECT_DOUBLE_EQ(qwl::resolve_wang_landau_schedule(p,8,3,4).flatness_treshold,1E10);
  p["FLATNESS_TRESHOLD"]="invalid";
  EXPECT_THROW(qwl::resolve_wang_landau_schedule(p,8,3,4),boost::bad_lexical_cast);
}

TEST(QwlSettings, DoesNotAddValidationForModificationFactor) {
  alps::Parameters p;
  p["INITIAL_MODIFICATION_FACTOR"]=1.;
  auto schedule=qwl::resolve_wang_landau_schedule(p,8,3,4);
  EXPECT_DOUBLE_EQ(schedule.logf,0.);
  EXPECT_TRUE(std::isinf(schedule.minimum_histogram));
  p["INITIAL_MODIFICATION_FACTOR"]=-1.;
  schedule=qwl::resolve_wang_landau_schedule(p,8,3,4);
  EXPECT_TRUE(std::isnan(schedule.logf));
  EXPECT_TRUE(std::isnan(schedule.minimum_histogram));
}
