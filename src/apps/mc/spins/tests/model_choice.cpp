// SPDX-License-Identifier: MIT
// Unit tests for spinmc::choose_simulation, which reads parameters only.
#include "model_choice.h"

#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using spinmc::simulation_kind;

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

simulation_kind choose(const std::string& text) {
    return spinmc::choose_simulation(parameters(text));
}

std::string choice_error(const std::string& text) {
    try {
        choose(text);
    } catch (const std::exception& error) {
        return error.what();
    }
    return "";
}

TEST(SpinmcModelChoice, CountsValuesSeparatedBySpaces) {
    EXPECT_EQ(spinmc::count_elements(""), 0);
    EXPECT_EQ(spinmc::count_elements("1"), 1);
    EXPECT_EQ(spinmc::count_elements("  1   2 "), 2);
    EXPECT_EQ(spinmc::count_elements("1\t2"), 1);
    EXPECT_EQ(spinmc::count_elements("2 * 0.5"), 3);
}

TEST(SpinmcModelChoice, TakesTheLargestCountOfTheKnownCouplingNames) {
    EXPECT_EQ(spinmc::find_dominant_matrix_string(parameters("T=1")), 0);
    EXPECT_EQ(spinmc::find_dominant_matrix_string(parameters(R"p(J=1;D3="1 2 3")p")), 3);
    EXPECT_EQ(spinmc::find_dominant_matrix_string(parameters(R"p(J=1;J'="1 2")p")), 2);
    EXPECT_EQ(spinmc::find_dominant_matrix_string(parameters(R"p(J=1;J4="1 2 3")p")), 1);
    EXPECT_EQ(spinmc::find_dominant_matrix_string(parameters(R"p(J=1;D4="1 2 3")p")), 1);
}

TEST(SpinmcModelChoice, ChoosesMomentAndMatrixType) {
    EXPECT_EQ(choose(R"p(MODEL="Ising";J=1)p"), simulation_kind::ising);
    EXPECT_EQ(choose(R"p(MODEL="O(4)";J="1 2")p"), simulation_kind::o4);
    EXPECT_EQ(choose(R"p(MODEL="Heisenberg";J=1)p"), simulation_kind::heisenberg_identity);
    EXPECT_EQ(choose(R"p(MODEL="Heisenberg";J=1;D="1 2 3")p"),
              simulation_kind::heisenberg_diagonal);
    EXPECT_EQ(choose(R"p(MODEL="Heisenberg";J="1 0 1 0 0 1")p"),
              simulation_kind::heisenberg_square);
    EXPECT_EQ(choose(R"p(MODEL="Heisenberg";J="1 0 0 0 1 0 0 0 1")p"),
              simulation_kind::heisenberg_square);
    EXPECT_EQ(choose(R"p(MODEL="XY";J=1)p"), simulation_kind::xy_identity);
    EXPECT_EQ(choose(R"p(MODEL="XY";J="1 2")p"), simulation_kind::xy_diagonal);
    EXPECT_EQ(choose(R"p(MODEL="XY";J="1 0 1")p"), simulation_kind::xy_square);
    EXPECT_EQ(choose(R"p(MODEL="XY";J="1 0 0 1")p"), simulation_kind::xy_square);
    EXPECT_EQ(choose(R"p(MODEL="Potts";q=3)p"), simulation_kind::potts3);
    EXPECT_EQ(choose(R"p(MODEL="Potts";q=4)p"), simulation_kind::potts4);
    EXPECT_EQ(choose(R"p(MODEL="Potts";q=10)p"), simulation_kind::potts10);
}

TEST(SpinmcModelChoice, ReportsUnsupportedChoices) {
    EXPECT_EQ(choice_error(R"p(MODEL="Ising")p"),
              "invalid parameters: cannot produce model Ising without J, J', J0-J3, D, D' or D0-D3");
    EXPECT_EQ(choice_error(R"p(MODEL="Ising";J="1 2";D=0)p"),
              R"(invalid parameters: cannot produce model Ising from 2 values in J="1 2")");
    EXPECT_EQ(choice_error(R"p(MODEL="Heisenberg";D="1 2")p"),
              R"(invalid parameters: cannot produce model Heisenberg from 2 values in D="1 2")");
    EXPECT_EQ(choice_error(R"p(MODEL="XY";J="1 2 3 4 5")p"),
              R"(invalid parameters: cannot produce model XY from 5 values in J="1 2 3 4 5")");
    EXPECT_EQ(choice_error(R"p(MODEL="Potts")p"), "parameter q not defined");
    EXPECT_EQ(choice_error(R"p(MODEL="Potts";q=5)p"),
              "5-state Potts model not implemented in file factory.C");
    EXPECT_EQ(choice_error("J=1"), "parameter MODEL not defined");
    EXPECT_EQ(choice_error(R"p(MODEL="Clock")p"), "Clock model not implemented in file factory.C");
}

TEST(SpinmcModelChoice, LeavesParametersUnchanged) {
    const auto p = parameters(R"p(MODEL="Heisenberg";J=1)p");
    spinmc::choose_simulation(p);
    EXPECT_EQ(p.size(), 2u);
}

} // namespace
