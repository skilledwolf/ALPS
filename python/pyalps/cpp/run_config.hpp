// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "dict_to_params.hpp"
#include <alps/run_config.hpp>
#include "scoped_signal_handlers.hpp"
namespace pyalps {
inline alps::run_configuration resolve_run(
    std::string_view schema, const nb::dict &parameters, const nb::dict &input,
    const nb::dict &output, const nb::dict &execution, const std::string &base = {}) {
    alps::run_configuration supplied;
    supplied.parameters = params_from_dict(parameters);
    supplied.input = params_from_dict(input);
    supplied.output = params_from_dict(output);
    supplied.execution = params_from_dict(execution);
    return alps::resolve_run_configuration(supplied, schema, base);
}
template<class Prepare, class Solve>
void bind_configured_solver(nb::module_ &module, std::string_view schema,
                            Prepare prepare, Solve solve) {
    module.def("schema", [schema] { return std::string(schema); });
    module.def("prepare", [schema, prepare](const nb::dict &parameters, const nb::dict &input,
               const nb::dict &output, const nb::dict &execution) {
        auto run = resolve_run(schema, parameters, input, output, execution);
        prepare(run);
        return run;
    }, nb::arg("parameters"), nb::arg("input") = nb::dict(),
       nb::arg("output") = nb::dict(), nb::arg("execution") = nb::dict());
    module.def("solve", [solve](const alps::run_configuration &run) {
        scoped_signal_handlers signal_handlers;
        return solve(run);
    }, nb::arg("run"));
    module.def("solve", [schema, solve](const std::string &filename) {
        scoped_signal_handlers signal_handlers;
        return solve(alps::load_run_configuration(filename, schema));
    }, nb::arg("filename"));
    module.def("solve", [schema, solve](const nb::dict &parameters, const nb::dict &input,
               const nb::dict &output, const nb::dict &execution) {
        scoped_signal_handlers signal_handlers;
        return solve(resolve_run(schema, parameters, input, output, execution));
    }, nb::arg("parameters"), nb::arg("input") = nb::dict(),
       nb::arg("output") = nb::dict(), nb::arg("execution") = nb::dict());
}
}
