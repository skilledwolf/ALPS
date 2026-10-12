// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include "default_model.hpp"
#include "control_validation.hpp"
#include <fstream>
namespace {
alps::params prior_controls() { alps::params p; p["OMEGA_MIN"] = 0.; p["OMEGA_MAX"] = 2.; return p; }
struct constant_density : Model { double operator()(double) override { return 1.; } };
}
TEST(MaxEntPrior, CumulativeMappingAndInverseMatchConstantDensity) {
    auto p = prior_controls(); p["BLOW_UP"] = 2.;
    GeneralDefaultModel prior(p, boost::shared_ptr<Model>(new constant_density));
    for (double t : {0., .0123, .25, .7, 1.}) {
        EXPECT_NEAR(prior.x(t), 2*t, 2.e-13);
        EXPECT_NEAR(prior.omega(prior.x(t)), 2*t, 2.e-13);
    }
}
TEST(MaxEntPrior, TabulationChecksShapeCoverageAndInterpolation) {
    const auto p = prior_controls();
    TabFunction prior(p, {-1., 0., 2., 3.}, {0., 1., 3., 0.});
    EXPECT_EQ(prior(0.),1.); EXPECT_EQ(prior(1.),2.); EXPECT_EQ(prior(2.),3.);
    EXPECT_THROW(prior(-2.), std::invalid_argument);
    EXPECT_THROW((TabFunction(p, {}, {})), std::invalid_argument);
    EXPECT_THROW((TabFunction(p, {0.,2.}, {0.,0.})), std::invalid_argument);
    EXPECT_THROW((TabFunction(p, {0.,0.,2.}, {1.,1.,1.})), std::invalid_argument);
    EXPECT_THROW((TabFunction(p, {0.,1.}, {1.,1.})), std::invalid_argument);
    EXPECT_THROW((TabFunction(p, {0.,2.}, {1.,-1.})), std::invalid_argument);
    alps::testing::TemporaryDirectory tmp; auto q = p;
    q["DEFAULT_MODEL"] = (tmp.path()/"prior.dat").string();
    for (const char* content : {"", "0 1\n2\n", "bad\n"}) {
        { std::ofstream out(q["DEFAULT_MODEL"].cast<std::string>()); out << content; }
        EXPECT_THROW((TabFunction(q,"DEFAULT_MODEL")), std::invalid_argument);
    }
}
TEST(MaxEntPrior, RejectsInvalidAnalyticPriorsAndRunControls) {
    auto p = prior_controls(); p["SIGMA"] = 0.;
    EXPECT_THROW(Gaussian{p}, std::invalid_argument);
    p["OMEGA_MIN"] = 2.; EXPECT_THROW(FlatDefaultModel{p}, std::invalid_argument);
    p = prior_controls(); p["N_ALPHA"] = 1; p["ALPHA_MIN"] = .1; p["ALPHA_MAX"] = 10.;
    EXPECT_THROW(maxent_detail::validate_controls(p), std::invalid_argument);
    p["N_ALPHA"] = 4; EXPECT_NO_THROW(maxent_detail::validate_controls(p));
    p["T"] = .5; p["BETA"] = 3.;
    EXPECT_THROW(maxent_detail::validate_controls(p), std::invalid_argument);
    p["BETA"] = 2.; p["ENFORCE_NORMALIZATION"] = true; p["SIGMA_NORMALIZATION"] = -1.;
    EXPECT_THROW(maxent_detail::validate_controls(p), std::invalid_argument);
}
