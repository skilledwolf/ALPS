/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_CLONE_PROXY_H
#define PARAPACK_CLONE_PROXY_H

#include "clone.h"

#include <alps/parapack/detail/clone_proxy.hpp>

namespace alps {
using clone_proxy = detail::clone_proxy<clone, Parameters, parapack::option>;
#ifdef ALPS_HAVE_MPI
using clone_proxy_mpi = detail::clone_proxy_mpi<clone_mpi, Parameters,
  clone_create_msg_t, parapack::option>;
#endif
} // namespace alps

#endif // PARAPACK_CLONE_PROXY_H
