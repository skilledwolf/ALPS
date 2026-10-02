// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "dict_to_params.hpp"
#include "scoped_signal_handlers.hpp"
#include <alps/maxent.hpp>
NB_MODULE(maxent_c,module) {
    module.def("AnalyticContinuation",[](const nanobind::dict& parameters,const nanobind::dict& input,
                                         const std::string& output,int time_limit,bool text_output) {
        pyalps::scoped_signal_handlers signal_handlers;
        const auto data=alps::maxent::read_data(pyalps::params_from_dict(input));
        alps::solvers::maxent(pyalps::params_from_dict(parameters),data,output,time_limit,text_output);
    },nanobind::arg("parameters"),nanobind::arg("input"),nanobind::arg("output_file"),
      nanobind::arg("time_limit")=60,nanobind::arg("text_output")=false);
}
