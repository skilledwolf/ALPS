// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// Legacy-file adaptation of reconciliation ef4b4ff36 tests/native_run.cpp.
// The four-site periodic Heisenberg spectrum is the independent oracle;
// compare eigenvalues and labels, not arbitrary degenerate eigenvectors.
#pragma once
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <fstream>
#include <map>

template<class Kernel> struct probe : Kernel {
  explicit probe(const boost::filesystem::path& path):Kernel(alps::ProcessList{},path) {}
  const auto& labels() const { return this->quantumnumbervalues_; }
  std::vector<double> energies() const {
    std::vector<double> result;
    for(const auto& sector:this->eigenvalues_) result.insert(result.end(),sector.begin(),sector.end());
    std::sort(result.begin(),result.end()); return result;
  }
};
namespace {
boost::filesystem::path input(const std::filesystem::path& directory,
                              std::map<std::string,std::string> extra={}) {
  std::map<std::string,std::string> p{{"MODEL","spin"},{"LATTICE","chain lattice"},
    {"L","4"},{"J","1"},{"CONSERVED_QUANTUMNUMBERS","Sz"},{"TRANSLATION_SYMMETRY","0"},
    {"MODEL_LIBRARY",std::string(ALPS_ED_TEST_RESOURCES)+"/models.xml"},
    {"LATTICE_LIBRARY",std::string(ALPS_ED_TEST_RESOURCES)+"/lattices.xml"}};
  for(const auto& item:extra) p[item.first]=item.second;
  const auto path=directory/"input.xml"; std::ofstream out(path);
  out<<"<SIMULATION><PARAMETERS>";
  for(const auto& item:p) out<<"<PARAMETER name=\""<<item.first<<"\">"<<item.second<<"</PARAMETER>";
  out<<"</PARAMETERS></SIMULATION>"; return path.string();
}
void near(const std::vector<double>& actual,const std::vector<double>& expected) {
  ASSERT_EQ(actual.size(),expected.size());
  for(std::size_t i=0;i<expected.size();++i) EXPECT_NEAR(actual[i],expected[i],2.e-10);
}
}
