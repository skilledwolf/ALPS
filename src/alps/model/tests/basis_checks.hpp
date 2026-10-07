// SPDX-License-Identifier: MIT
#pragma once

#include <gtest/gtest.h>
#include <alps/model.h>
#include <alps/parser/xslt_path.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace model_test {
inline alps::ModelLibrary library() {
  std::ifstream input(alps::search_xml_library_path("models.xml"));
  if (!input) throw std::runtime_error("Cannot open the model library");
  return alps::ModelLibrary(input);
}

inline alps::SiteBasisDescriptor<short> descriptor(
    const std::string& name, const alps::Parameters& parameters = {}) {
  auto result = library().get_site_basis(name);
  result.set_parameters(parameters);
  return result;
}

using States = std::vector<std::vector<double>>;

template<class State>
void expect_states(const alps::site_basis<short, State>& basis,
                   const std::vector<std::string>& names, const States& expected) {
  ASSERT_EQ(basis.basis().size(), names.size());
  ASSERT_EQ(basis.size(), expected.size());
  ASSERT_FALSE(basis.empty());
  EXPECT_TRUE(basis.check_sort());
  for (std::size_t column = 0; column < names.size(); ++column)
    EXPECT_EQ(basis.basis()[column].name(), names[column]);
  for (std::size_t row = 0; row < expected.size(); ++row) {
    SCOPED_TRACE(::testing::Message() << "state " << row);
    ASSERT_EQ(expected[row].size(), names.size());
    EXPECT_EQ(basis.index(basis[row]), row);
    for (std::size_t column = 0; column < names.size(); ++column)
      EXPECT_EQ(alps::get_quantumnumber(basis[row], column),
                alps::half_integer<short>(expected[row][column])) << names[column];
  }
}

// Enumerate the historical spinful-boson contract in integer quantum numbers:
// 0 <= N <= cutoff, 0 <= J <= N*spin, and -J <= Jz <= J, in that order.
inline States spinful_states(int cutoff, int spin) {
  States states;
  for (int n = 0; n <= cutoff; ++n)
    for (int j = 0; j <= n * spin; ++j)
      for (int m = -j; m <= j; ++m)
        states.push_back({double(n), double(j), double(m)});
  return states;
}
} // namespace model_test
