// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include "run_config.hpp"
#include <alps/ctint.hpp>

namespace nb = nanobind;

NB_MODULE(ctint, module) {
  module.def("schema", [](const nb::dict &parameters) {
    return alps::ctint::schema(pyalps::params_from_dict(parameters));
  }, nb::arg("parameters") = nb::dict());
  module.def("prepare", [](const nb::dict &parameters, const nb::dict &input,
              const nb::dict &output, const nb::dict &execution) {
    auto run = pyalps::resolve_run(alps::ctint::schema(pyalps::params_from_dict(parameters)),
                                  parameters, input, output, execution);
    alps::ctint::prepare_run(run);
    return run;
  }, nb::arg("parameters"), nb::arg("input") = nb::dict(),
     nb::arg("output") = nb::dict(), nb::arg("execution") = nb::dict());
  module.def("solve", [](const alps::run_configuration &run) {
    pyalps::scoped_signal_handlers signal_handlers;
    alps::solvers::ctint(run);
  }, nb::arg("run"));
}
