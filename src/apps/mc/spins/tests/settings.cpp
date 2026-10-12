// SPDX-License-Identifier: MIT
// Unit tests for spinmc::resolve_settings. Settings resolve from parameters
// alone, without a lattice, a scheduler or Monte Carlo sweeps.
#include "settings.h"

#include <alps/expression.h>
#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace {

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

std::string resolution_error(const std::string& text) {
    try {
        spinmc::resolve_settings<3>(parameters(text));
    } catch (const std::exception& error) {
        return error.what();
    }
    return "";
}

TEST(SpinmcSettings, AppliesDefaults) {
    const auto s = spinmc::resolve_settings<3>(parameters("T=2;SWEEPS=95"));
    EXPECT_DOUBLE_EQ(s.beta, 0.5);
    EXPECT_EQ(s.sweeps, 95u);
    EXPECT_EQ(s.thermalization, 9u); // SWEEPS/10 in integer arithmetic
    EXPECT_FALSE(s.quantum_convention);
    EXPECT_DOUBLE_EQ(s.g, 1.);
    EXPECT_EQ(s.print_sweeps, 0);
    EXPECT_FALSE(s.use_error_limit);
    EXPECT_EQ(s.update, spinmc::update_request::automatic);
    EXPECT_EQ(s.update_text, "");
    EXPECT_EQ(s.default_coupling, "1.0");
    for (int i = 0; i < 3; ++i)
        EXPECT_DOUBLE_EQ(s.h[i], 0.);
}

TEST(SpinmcSettings, PrefersBetaAndRequiresTemperatureFirst) {
    EXPECT_DOUBLE_EQ(spinmc::resolve_settings<3>(parameters("beta=2;T=100;SWEEPS=1")).beta, 2.);
    EXPECT_EQ(resolution_error("SWEEPS=1"), "parameter T not defined");
    EXPECT_EQ(resolution_error("J=1"), "parameter T not defined");
    EXPECT_EQ(resolution_error("T=1"), "parameter SWEEPS not defined");
}

TEST(SpinmcSettings, EvaluatesFieldAndGyromagneticExpressions) {
    const auto p = parameters(R"(T=1;SWEEPS=10;J=0.5;h="0.2*J";g="2*J";THERMALIZATION=3)");
    const auto s = spinmc::resolve_settings<3>(p);
    EXPECT_DOUBLE_EQ(s.h[0], 0.);
    EXPECT_DOUBLE_EQ(s.h[1], 0.);
    EXPECT_DOUBLE_EQ(s.h[2], 0.1); // one value is the last component
    EXPECT_DOUBLE_EQ(s.g, 1.);
    EXPECT_EQ(s.thermalization, 3u);
    EXPECT_EQ(s.default_coupling, "0.5");

    const auto full = spinmc::resolve_settings<3>(parameters(R"(T=1;SWEEPS=10;h="1 2 3")"));
    EXPECT_DOUBLE_EQ(full.h[0], 1.);
    EXPECT_DOUBLE_EQ(full.h[1], 2.);
    EXPECT_DOUBLE_EQ(full.h[2], 3.);
    EXPECT_EQ(resolution_error(R"(T=1;SWEEPS=10;h="1 2")"), "too few input values");
}

TEST(SpinmcSettings, ClassifiesUpdateRequests) {
    const auto update = [](const std::string& value) {
        return spinmc::resolve_settings<1>(parameters("T=1;SWEEPS=1;UPDATE=" + value));
    };
    EXPECT_EQ(update(R"("local")").update, spinmc::update_request::local);
    EXPECT_EQ(update(R"("cluster")").update, spinmc::update_request::cluster);
    EXPECT_EQ(update(R"("")").update, spinmc::update_request::automatic);
    const auto unknown = update(R"("fast")");
    EXPECT_EQ(unknown.update, spinmc::update_request::unrecognized);
    EXPECT_EQ(unknown.update_text, "fast");
}

TEST(SpinmcSettings, RequiresBothErrorLimitParameters) {
    EXPECT_FALSE(spinmc::resolve_settings<1>(
        parameters(R"(T=1;SWEEPS=1;ERROR_VARIABLE="Energy")")).use_error_limit);
    EXPECT_TRUE(spinmc::resolve_settings<1>(
        parameters(R"(T=1;SWEEPS=1;ERROR_VARIABLE="Energy";ERROR_LIMIT=0.01)")).use_error_limit);
}

TEST(SpinmcSettings, LeavesParametersUnchanged) {
    const auto p = parameters("T=1;SWEEPS=10");
    spinmc::resolve_settings<3>(p);
    EXPECT_EQ(p.size(), 2u);
    EXPECT_FALSE(p.defined("UPDATE"));
}

// Current behavior, pinned so that a change is deliberate: expression
// evaluation during resolution uses the process-wide disorder generator.
// A random g consumes a draw from it.
TEST(SpinmcSettings, RandomExpressionsDrawFromTheDisorderGenerator) {
    alps::Disorder::seed(5);
    const double first = alps::Disorder::random();
    const double second = alps::Disorder::random();

    alps::Disorder::seed(5);
    const auto s = spinmc::resolve_settings<1>(parameters(R"p(T=1;SWEEPS=1;g="random()")p"));
    EXPECT_DOUBLE_EQ(s.g, first);
    EXPECT_DOUBLE_EQ(alps::Disorder::random(), second);
}

// Constructing the parameter evaluator also reseeds that generator when
// DISORDERSEED differs from its last seed, even for a non-random g.
TEST(SpinmcSettings, EvaluationReseedsTheDisorderGenerator) {
    alps::Disorder::seed(9);
    const double after_seed_9 = alps::Disorder::random();

    alps::Disorder::seed(5);
    spinmc::resolve_settings<1>(parameters("T=1;SWEEPS=1;DISORDERSEED=9;g=2"));
    EXPECT_DOUBLE_EQ(alps::Disorder::random(), after_seed_9);
}

} // namespace
