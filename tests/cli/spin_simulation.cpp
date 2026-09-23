// SPDX-License-Identifier: MIT
#include "abstractspinsim.h"
#include "spinsim.h"
#include <alps/osiris/xdrdump.h>
#include <cmath>
#include <iomanip>
#include <stdexcept>

template<class Moment, class Matrix>
void check(std::string const& name, std::string const& update, std::string const& coupling,
           std::string const& onsite, std::string const& field, std::string const& output) {
  alps::Parameters params;
  params["LATTICE_LIBRARY"] = ALPS_TEST_LATTICES;
  params["LATTICE"] = "chain lattice";
  params["L"] = 8;
  params["T"] = 1.75;
  params["SWEEPS"] = 128;
  params["THERMALIZATION"] = 0;
  params["SEED"] = 1729;
  params["UPDATE"] = update;
  params["J"] = coupling;
  if (!onsite.empty()) params["D"] = onsite;
  if (!field.empty()) params["h"] = field;
  SpinSim<Moment, Matrix> sim(alps::ProcessList(), params, 0);
  for (int step = 0; step < 128; ++step) sim.dostep();
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
    check<IsingMoment, MIdMatrix<double, 1>>("ising-" + update, update, "1", ".3", "", output);
    check<XYMoment, MIdMatrix<double, 2>>("xy-" + update, update, "1", ".2", "", output);
    check<HeisenbergMoment, MIdMatrix<double, 3>>("heisenberg-" + update, update, "1", "", "", output);
    check<ONMoment<4>, MIdMatrix<double, 4>>("o4-" + update, update, "1", "", "", output);
    check<PottsMoment<3>, MIdMatrix<double, 2>>("potts3-" + update, update, "1", "", "", output);
    check<PottsMoment<4>, MIdMatrix<double, 2>>("potts4-" + update, update, "1", "", "", output);
    check<PottsMoment<10>, MIdMatrix<double, 2>>("potts10-" + update, update, "1", "", "", output);
  }
  check<IsingMoment, MIdMatrix<double, 1>>("ising-field", "local", "1", ".3", ".2", output);
  check<XYMoment, MIdMatrix<double, 2>>("xy-field", "local", "1", ".2", ".1 .2", output);
  check<XYMoment, DiagMatrix<double, 2>>("xy-diagonal", "local", "1 .8", ".2 .1", ".1 .2", output);
  check<XYMoment, SquareMatrix<double, 2>>("xy-square", "local", "1 .1 .1 .8", ".2 .1 .1 .3", ".1 .2", output);
  check<HeisenbergMoment, DiagMatrix<double, 3>>("heisenberg-diagonal", "local", "1 .8 .6", ".2 .1 .3", "", output);
  check<HeisenbergMoment, SquareMatrix<double, 3>>("heisenberg-square", "local", "1 .1 .2 .1 .8 .1 .2 .1 .6", "", "", output);
}
