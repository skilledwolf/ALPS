// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include "maxent_parms.hpp"
#include <fstream>

namespace {
alps::params controls() {
    alps::params p; p["NDAT"] = 4; p["NFREQ"] = 4; p["T"] = .5;
    p["OMEGA_MAX"] = 4.; p["FREQUENCY_GRID"] = "linear"; return p;
}
alps::maxent::data valid() { return {{-.1, -.2, -.2, -.1}, {.01, .01, .01, .01}, {}, {0., .5, 1., 2.}}; }
struct covariance_probe : ContiParameters {
    covariance_probe(const alps::params& p, const alps::maxent::data& d) : ContiParameters(p,d) {
        vector_type frequencies(4); for (int i=0; i<4; ++i) frequencies[i] = i+.5;
        setup_kernel(p, 4, frequencies);
    }
};
}
TEST(MaxEntData, OwningValuesTemperatureOnlyAndDefaultNorm) {
    auto d = valid(); const auto p = controls(); ContiParameters c(p,d);
    d.values[0] = 99.; EXPECT_EQ(c.y(0), -.1); EXPECT_EQ(c.T(), .5);
    d = valid(); d.tau[1] = 3.; EXPECT_THROW(ContiParameters(p,d), std::invalid_argument);
}
TEST(MaxEntData, RejectsInvalidShapesAndNumbers) {
    auto d = valid(); d.values.resize(3); EXPECT_THROW(alps::maxent::validate_data(d), std::invalid_argument);
    d = valid(); d.errors.resize(3); EXPECT_THROW(alps::maxent::validate_data(d), std::invalid_argument);
    d = valid(); d.errors[0] = 0; EXPECT_THROW(alps::maxent::validate_data(d), std::invalid_argument);
    d = valid(); d.values[0] = std::numeric_limits<double>::infinity(); EXPECT_THROW(alps::maxent::validate_data(d), std::invalid_argument);
    d = valid(); d.covariance.resize(16); EXPECT_THROW(alps::maxent::validate_data(d), std::invalid_argument);
    d.errors.clear(); d.covariance[1] = .5; EXPECT_THROW(alps::maxent::validate_data(d), std::invalid_argument);
    d.covariance.assign(15, 0.); EXPECT_THROW(alps::maxent::validate_data(d), std::invalid_argument);
}
TEST(MaxEntData, CovarianceMustHaveNonnegativeUsableSpectrum) {
    auto d = valid(); d.errors.clear(); d.covariance.assign(16,0.);
    auto p = controls();
    EXPECT_THROW((covariance_probe(p,d)), std::invalid_argument);
    for (int i=0; i<4; ++i) d.covariance[5*i] = 1.;
    d.covariance[1] = d.covariance[4] = 2.;
    EXPECT_THROW((covariance_probe(p,d)), std::invalid_argument);
    d.covariance[1] = d.covariance[4] = .2;
    EXPECT_NO_THROW((covariance_probe(p,d)));
}
TEST(MaxEntData, ShortArchivesAndMalformedIndexedTextFail) {
    alps::testing::TemporaryDirectory tmp;
    auto p = controls(); const auto path = (tmp.path()/"input.h5").string();
    p["DATA"] = path; p["DATA_IN_HDF5"] = true;
    { alps::hdf5::archive ar(path,"w"); ar["/Data"] << std::vector<double>{1,2}; ar["/Error"] << valid().errors; }
    EXPECT_THROW(maxent_detail::legacy_data(p), std::invalid_argument);
    p["DATA_IN_HDF5"] = false; p["DATA"] = (tmp.path()/"input.txt").string();
    for (const char* text : {"-1 0 1\n", "0 1 1\n0 1 1\n", "0 1\n", "0 1 1\n", "4 1 1\n"}) {
        { std::ofstream out(p["DATA"].cast<std::string>()); out << text; }
        EXPECT_THROW(maxent_detail::legacy_data(p), std::invalid_argument);
    }
    { std::ofstream out(p["DATA"].cast<std::string>()); out << "3 4 .1\n1 2 .1\n0 1 .1\n2 3 .1\n"; }
    EXPECT_EQ(maxent_detail::legacy_data(p).values, (std::vector<double>{1,2,3,4}));
    p["NDAT"] = -1; EXPECT_THROW(maxent_detail::legacy_data(p), std::invalid_argument);
}
TEST(MaxEntData, MoreDataThanFrequenciesUsesRectangularSvd) {
    auto p = controls(); p["NDAT"] = 8;
    for(int i=0;i<8;++i) { p["X_"+std::to_string(i)] = -.1; p["SIGMA_"+std::to_string(i)] = .01; }
    EXPECT_NO_THROW(MaxEntParameters{p});
}
