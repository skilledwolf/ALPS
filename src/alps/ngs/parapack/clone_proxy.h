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

#ifndef NGS_PARAPACK_CLONE_PROXY_H
#define NGS_PARAPACK_CLONE_PROXY_H

#include <alps/ngs/parapack/clone.h>

#include <alps/parapack/detail/clone_proxy.hpp>

namespace alps { namespace ngs_parapack {
using clone_proxy = alps::detail::clone_proxy<clone, alps::params,
  dump_policy_t, clone_timer::duration_t>;
#ifdef ALPS_HAVE_MPI
using clone_proxy_mpi = alps::detail::clone_proxy_mpi<clone_mpi, alps::params,
  clone_create_msg_t, dump_policy_t, clone_timer::duration_t>;
#endif
}} // namespace alps::ngs_parapack

#endif // NGS_PARAPACK_CLONE_PROXY_H
