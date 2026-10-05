/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "loop_worker.h"
#include <looper/evaluator_impl.h>
#include <looper/version.h>
#include <alps/parapack/parapack.h>
#include <alps/parapack/exchange.h>

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_COPYRIGHT(LOOPER_COPYRIGHT)
PARAPACK_SET_VERSION(LOOPER_VERSION_STRING)

namespace {
using continuous_time_worker = looper::loop_worker<looper::path_integral>;
using sse_worker = looper::loop_worker<looper::sse>;
using loop_evaluator = looper::evaluator<loop_config::measurement_set>;

PARAPACK_REGISTER_ALGORITHM(continuous_time_worker, "loop");
PARAPACK_REGISTER_ALGORITHM(continuous_time_worker, "loop; path integral");
PARAPACK_REGISTER_ALGORITHM(sse_worker, "loop; sse");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<continuous_time_worker>, "loop; exchange");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<continuous_time_worker>, "loop; path integral; exchange");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<sse_worker>, "loop; sse; exchange");
PARAPACK_REGISTER_EVALUATOR(loop_evaluator, "loop");
} // namespace
