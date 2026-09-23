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

#ifndef NGS_PARAPACK_JOB_P_H
#define NGS_PARAPACK_JOB_P_H

#include "job.h"
#include <alps/parapack/detail/job_xml.hpp>
#include <alps/ngs/parapack/params_p.h>

namespace alps { namespace ngs_parapack {
using job_xml_writer = alps::parapack::detail::job_xml_writer<task>;
using job_task_xml_handler = alps::parapack::detail::job_task_xml_handler<task>;
using job_tasks_xml_handler = alps::parapack::detail::job_tasks_xml_handler<task>;
using alps::parapack::detail::filename_xml_handler;
using alps::parapack::detail::version_xml_handler;
}
}

#endif
