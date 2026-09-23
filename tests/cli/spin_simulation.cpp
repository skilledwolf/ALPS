// SPDX-License-Identifier: MIT
#include "abstractspinsim.h"
#include "spinsim.h"
#include "factory.h"
#include <memory>
#include <sstream>
#include <alps/osiris/xdrdump.h>
#include <cmath>
#include <iomanip>
#include <stdexcept>

template<class Moment, class Matrix>
void check(std::string const& model, std::string const& name, std::string const& update, std::string const& coupling,
           std::string const& onsite, std::string const& field, std::string const& output, int q = 0) {
  alps::Parameters params;
  params["LATTICE_LIBRARY"] = ALPS_TEST_LATTICES;
  params["LATTICE"] = "chain lattice";
  params["L"] = 8;
  params["T"] = 1.75;
  params["SWEEPS"] = 128;
  params["THERMALIZATION"] = 0;
  params["SEED"] = 1729;
  params["UPDATE"] = update;
  params["MODEL"] = model;
  if (q) params["q"] = q;
  params["J"] = coupling;
  if (!onsite.empty()) params["D"] = onsite;
  if (!field.empty()) params["h"] = field;
  SpinSim<Moment, Matrix> sim(alps::ProcessList(), params, 0);
  std::unique_ptr<alps::scheduler::Worker> worker(SpinFactory().make_worker({}, params, 0));
  auto* selected = dynamic_cast<SpinSim<Moment, Matrix>*>(worker.get());
  if (!selected) throw std::runtime_error("factory chose the wrong spin or matrix type: " + name);
  for (int step = 0; step < 128; ++step) { sim.dostep(); selected->dostep(); }
  std::ostringstream expected, actual;
  expected << std::setprecision(17) << sim.get_measurements();
  actual << std::setprecision(17) << selected->get_measurements();
  if (actual.str() != expected.str()) throw std::runtime_error("factory changed spin measurements: " + name);
  alps::RealObsevaluator energy(sim.get_measurements()["Energy"]);
  alps::RealObsevaluator magnetization(sim.get_measurements()["Magnetization^2"]);
  if (energy.count() != 128 || !std::isfinite(energy.mean()) ||
      !std::isfinite(magnetization.mean()) || magnetization.mean() < 0)
    throw std::runtime_error("invalid spin measurements: " + name);
  // Optional output supports comparing complete statistics and serialized spin
  // configurations against another build with exactly the same random seed.
  if (!output.empty()) {
    std::cout << name << '\n' << std::setprecision(17) << sim.get_measurements();
    alps::OXDRFileDump dump(output + "/" + name + ".dump");
    sim.save(dump);
  }
}

int main(int argc, char** argv) {
  std::string output = argc > 1 ? argv[1] : "";
  for (std::string const update : {"local", "cluster"}) {
    check<IsingMoment, MIdMatrix<double, 1>>("Ising", "ising-" + update, update, "1", ".3", "", output);
    check<XYMoment, MIdMatrix<double, 2>>("XY", "xy-" + update, update, "1", ".2", "", output);
    check<HeisenbergMoment, MIdMatrix<double, 3>>("Heisenberg", "heisenberg-" + update, update, "1", "", "", output);
    check<ONMoment<4>, MIdMatrix<double, 4>>("O(4)", "o4-" + update, update, "1", "", "", output);
    check<PottsMoment<3>, MIdMatrix<double, 2>>("Potts", "potts3-" + update, update, "1", "", "", output, 3);
    check<PottsMoment<4>, MIdMatrix<double, 2>>("Potts", "potts4-" + update, update, "1", "", "", output, 4);
    check<PottsMoment<10>, MIdMatrix<double, 2>>("Potts", "potts10-" + update, update, "1", "", "", output, 10);
  }
  check<IsingMoment, MIdMatrix<double, 1>>("Ising", "ising-field", "local", "1", ".3", ".2", output);
  check<XYMoment, MIdMatrix<double, 2>>("XY", "xy-field", "local", "1", ".2", ".1 .2", output);
  check<XYMoment, DiagMatrix<double, 2>>("XY", "xy-diagonal", "local", "1 .8", ".2 .1", ".1 .2", output);
  check<XYMoment, SquareMatrix<double, 2>>("XY", "xy-square", "local", "1 .1 .1 .8", ".2 .1 .1 .3", ".1 .2", output);
  check<HeisenbergMoment, DiagMatrix<double, 3>>("Heisenberg", "heisenberg-diagonal", "local", "1 .8 .6", ".2 .1 .3", "", output);
  check<HeisenbergMoment, SquareMatrix<double, 3>>("Heisenberg", "heisenberg-square", "local", "1 .1 .2 .1 .8 .1 .2 .1 .6", "", "", output);
  check<XYMoment, SquareMatrix<double, 2>>("XY", "xy-triangular", "local", "1 .1 .8", "", "", output);
  check<HeisenbergMoment, SquareMatrix<double, 3>>("Heisenberg", "heisenberg-triangular", "local", "1 .1 .2 .8 .1 .6", "", "", output);
  check<XYMoment, DiagMatrix<double, 2>>("XY", "xy-onsite", "local", "1", ".2 .1", "", output);
  check<HeisenbergMoment, DiagMatrix<double, 3>>("Heisenberg", "heisenberg-onsite", "local", "1", ".2 .1 .3", "", output);
  for (std::string model : {"Ising", "XY", "Heisenberg"}) {
    alps::Parameters invalid;
    invalid["MODEL"] = model;
    invalid["J"] = "1 2 3 4 5";
    invalid["D"] = "0";
    bool rejected = false;
    try { std::unique_ptr<alps::scheduler::Worker> worker(SpinFactory().make_worker({}, invalid, 0)); }
    catch (std::runtime_error const&) { rejected = true; }
    if (!rejected) throw std::runtime_error("factory accepted an invalid matrix size: " + model);
  }
}
