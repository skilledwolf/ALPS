// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
// A unit-spin self bond is exactly constant and cannot affect its trajectory.
#include "ising.h"
#include "xy.h"
#include "heisenberg.h"
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <fstream>
namespace {
template<class Worker> void compare(const std::string& magnetization) {
    alps::testing::TemporaryDirectory dir;const auto path=(dir.path()/"loops.xml").string();
    {std::ofstream out(path);out<<"<LATTICES><GRAPH name=\"loop\" vertices=\"1\"><VERTEX id=\"1\"/>"
      "<EDGE source=\"1\" target=\"1\"/></GRAPH>"
      "<GRAPH name=\"plain\" vertices=\"1\"><VERTEX id=\"1\"/></GRAPH></LATTICES>";}
    for(double j:{-.6,.6,1.e20}) {
        alps::Parameters p;p["LATTICE_LIBRARY"]=path;p["LATTICE"]="loop";
        p["T"]=1.2;p["H"]=.4;p["J"]=j;p["SEED"]=173;p["WORKER_SEED"]=173;p["DISORDERSEED"]=173;
        p["BASE_NAME"]="unused";p["CLONE_ID"]=0;p["SWEEPS"]=1000;p["THERMALIZATION"]=0;
        Worker loop(p);alps::ObservableSet a;loop.init_observables(p,a);
        p["LATTICE"]="plain";Worker plain(p);alps::ObservableSet b;plain.init_observables(p,b);
        for(unsigned step=0;step<256;++step) {
            loop.run(a);plain.run(b);
            EXPECT_DOUBLE_EQ(alps::RealObsevaluator(a[magnetization]).mean(),alps::RealObsevaluator(b[magnetization]).mean());
            if(std::abs(j)<1)EXPECT_NEAR(alps::RealObsevaluator(a["Energy"]).mean()-alps::RealObsevaluator(b["Energy"]).mean(),-j,1.e-12);
        }
    }
}
}
TEST(SimpleLegacySelfLoops, Ising){compare<ising_worker>("Magnetization Density");}
TEST(SimpleLegacySelfLoops, XY){compare<xy_worker>("Magnetization Density X");}
TEST(SimpleLegacySelfLoops, Heisenberg){compare<heisenberg_worker>("Magnetization Density Z");}
