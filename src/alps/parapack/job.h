/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2013 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_JOB_H
#define PARAPACK_JOB_H

#include "clone_info.h"
#include <alps/parapack/detail/task_state.hpp>
#include "integer_range.h"
#include "logger.h"
#include "option.h"
#include "types.h"
#include <alps/parser/xmlstream.h>
#include <alps/alea/observableset.h>
#include <boost/optional.hpp>
#include <deque>
#include <set>
#include <vector>

namespace alps {

class ALPS_DECL task : public alps::parapack::detail::task_state<task, Parameters> {
  typedef alps::parapack::detail::task_state<task, Parameters> base_type;
public:
  task() = default;
  explicit task(boost::filesystem::path const& file) : base_type(file) {}
  void load();
  void save(alps::parapack::option const& opt) const;
  void halt();
  void check_parameter(alps::parapack::option const& opt);
  void write_xml_archive(oxstream& os) const;
  void save_observable(alps::parapack::option const& opt) const;
  void evaluate(alps::parapack::option const& opt);

private:
  friend base_type;
  void clear_observable_cache() { obs_.clear(); }
  std::vector<ObservableSet> obs_;
};


} // end namespace alps

#endif // PARAPACK_JOB_H
