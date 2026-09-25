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

#ifndef NGS_PARAPACK_JOB_H
#define NGS_PARAPACK_JOB_H

#include <alps/ngs/parapack/clone_info.h>
#include <alps/parapack/detail/task_state.hpp>
#include <alps/parapack/integer_range.h>
#include <alps/parapack/logger.h>
#include <alps/parapack/types.h>
#include <alps/parser/xmlstream.h>
// #include <alps/scheduler.h>
#include <alps/config.h>
#include <boost/optional.hpp>
#include <deque>
#include <set>
#include <vector>

namespace alps {
namespace ngs_parapack {

class ALPS_DECL task : public alps::parapack::detail::task_state<task, alps::params> {
  typedef alps::parapack::detail::task_state<task, alps::params> base_type;
public:
  task() = default;
  explicit task(boost::filesystem::path const& file) : base_type(file) {}
  void load();
  void save(bool write_xml) const;
  void halt();
  void check_parameter(bool write_xml);
  void write_xml_archive(oxstream& os) const;

private:
  friend base_type;
  void clear_observable_cache() {  }
};


} // end namespace ngs_parapack
} // end namespace alps

#endif // NGS_PARAPACK_JOB_H
