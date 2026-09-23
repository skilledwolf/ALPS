/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>,
*                            Ryo Igarashi <rigarash@issp.u-tokyo.ac.jp>,
*                            Haruhiko Matsuo <halm@rist.or.jp>,
*                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Yuichi Motoyama <yomichi@looper.t.u-tokyo.ac.jp>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef NGS_PARAPACK_CLONE_INFO_H
#define NGS_PARAPACK_CLONE_INFO_H

#include <alps/parapack/clone_info.h>
#include <alps/ngs/params.hpp>

namespace alps { namespace ngs_parapack {
using alps::clone_phase;
using alps::clone_info;
#ifdef ALPS_HAVE_MPI
using alps::clone_info_mpi;
#endif
} }

#endif
