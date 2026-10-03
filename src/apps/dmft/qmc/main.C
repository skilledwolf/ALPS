/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2026 by Emanuel Gull <gull@phys.columbia.edu>
 *               2005 - 2009 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/


#include "hirschfyesim.h"
#include "selfconsistency.h"
#include "externalsolver.h"
#include "hilberttransformer.h"
#include "interaction_expansion_choice.h"
#include "interaction_expansion/interaction_expansion.hpp"
#include "run_config.h"
#include <alps/utility/copyright.hpp>
#include <memory>
#include <iostream>

namespace {
void print_citation() {
  std::cout << "ALPS DMFT framework for the single site impurity problem.       "<<std::endl;
  std::cout << "  For further information see the ALPS DMFT paper:              "<<std::endl;
  std::cout << "  Computer Physics Communications 182, 1078 (2011)              "<<std::endl;
  std::cout << "                                                                "<<std::endl;
  std::cout << "  copyright (c) 2005-2026 by the ALPS collaboration.            "<<std::endl;
  std::cout << std::endl;
  alps::print_copyright(std::cout);

  std::cout << "****************************************************************"<<std::endl;
  std::cout << "* Recommended citation in scientific publications:             *"<<std::endl;
  std::cout << "* This code used the ALPS [1] DMFT framework [2]               *"<<std::endl;
  std::cout << "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011)            *"<<std::endl;
  std::cout << "****************************************************************"<<std::endl;
}
// Band moments set by the Hilbert transformers are recorded as derived values.
void record_derived(alps::run_configuration& run) {
  for (const auto& [name, value] : run.parameters)
    if (!run.origins.count("parameters." + name)) run.origins["parameters." + name] = "derived";
}
}

int main(int argc, char** argv) {
  try {
    bool validate = false, show_schema = false;
    std::string filename;
    for (int i = 1; i < argc; ++i) {
      const std::string arg(argv[i]);
      if (arg == "--help" || arg == "-h") {
        std::cout << "Usage: dmft [--validate|--schema] run.toml\n";
        return 0;
      }
      if (arg == "--schema") show_schema = true;
      else if (arg == "--validate") validate = true;
      else if (arg.empty() || arg.front() == '-') throw std::invalid_argument("Unknown option: " + arg);
      else if (filename.empty()) filename = arg;
      else throw std::invalid_argument("Expected one TOML run file");
    }
    if (show_schema) { std::cout << alps::dmft::schema_for_run(filename); return 0; }
    if (filename.empty()) throw std::invalid_argument("No TOML run file specified");
    auto run = alps::dmft::load_run(filename);
    auto& p = run.parameters;
    auto initial_parameters = p;
    if (p.exists("H_INIT")) initial_parameters["H"] = p["H_INIT"];
    const auto ready = [&] {
      record_derived(run);
      if (validate) std::cout << "Valid DMFT configuration: " << filename << '\n';
      else print_citation();
      return !validate;
    };
    if (run.execution["loop"].as<std::string>() == "tau") {
      SemicircleHilbertTransformer transform(p);
      auto initial = transform.initial_G0(initial_parameters, run.input);
      require_finite(initial, "DMFT initial Green function");
      if (!ready()) return 0;
      ExternalSolver solver(run);
      F_selfconsistency_loop(run, solver, transform, std::move(initial));
    } else {
      std::unique_ptr<FrequencySpaceHilbertTransformer> transform;
      if (p["SEMICIRCLE_HILBERT"].as<bool>()) transform.reset(new SemicircleFSHilbertTransformer(p));
      else transform.reset(new GeneralFSHilbertTransformer(p, run.input));
      auto initial = transform->initial_G0(initial_parameters, run.input);
      require_finite(initial, "DMFT initial Green function");
      if (validate) {
        boost::shared_ptr<FourierTransformer> fourier;
        FourierTransformer::generate_transformer(initial_parameters, fourier);
      }
      if (!ready()) return 0;
      std::unique_ptr<MatsubaraImpuritySolver> solver;
      alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> hf;
      alps::scheduler::BasicFactory<InteractionExpansionSim,HubbardInteractionExpansionRun> ss;
      alps::scheduler::BasicFactory<InteractionExpansionSim,MultiBandDensityHubbardInteractionExpansionRun> mb;
      switch (alps::dmft::selected_solver(run)) {
        case alps::dmft::solver_kind::hirsch_fye:
          solver.reset(new alps::ImpuritySolver(hf, run, argc, argv));
          break;
        case alps::dmft::solver_kind::interaction_expansion:
          // load_run rejects the unsupported dimensions.
          if (select_interaction_expansion(p["FLAVORS"].as<int>(), p["SITES"].as<int>()) ==
              interaction_expansion_choice::single_site_hubbard)
            solver.reset(new alps::ImpuritySolver(ss, run, argc, argv));
          else solver.reset(new alps::ImpuritySolver(mb, run, argc, argv));
          break;
        default:
          solver.reset(new ExternalSolver(run));
      }
      selfconsistency_loop_omega(run, *solver, *transform, std::move(initial));
    }
    alps::hdf5::archive archive(run.output["results"].as<std::string>(), "a");
    archive["/run_config"] << run;
    archive["/parameters"] << p;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "dmft: " << error.what() << '\n'; return 1;
  }
}
