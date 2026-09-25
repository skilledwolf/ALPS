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

#include <alps/ngs/parapack/job.h>
#include <alps/ngs/parapack/clone.h>
#include <alps/parapack/filelock.h>
// #include "measurement.h"
#include <alps/ngs/parapack/simulation_p.h>

#include <alps/hdf5.hpp>

// some file (probably a python header) defines a tolower macro ...
#undef tolower
#undef toupper

#include <boost/config.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace alps {
namespace ngs_parapack {

//
// task
//

void task::load() {
  if (on_memory()) boost::throw_exception(std::logic_error("task::load() task already loaded"));
  params_ = alps::params();
  clone_info_.clear();
  boost::filesystem::path file_in = absolute(boost::filesystem::path(file_in_str_), basedir_);
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str_), basedir_);
  simulation_xml_handler handler(params_, clone_info_);
  XMLParser parser(handler);
  if (!exists(file_out)) {
    parser.parse(file_in);
  } else {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    parser.parse(file_out);
  }

  num_clones_ = (params_["NUM_CLONES"] | 1);
  restore_state();
}

void task::save(bool write_xml) const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str_), basedir_);
  {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    if (exists(file_out)) {
      alps::params params_tmp;
      std::deque<clone_info> clone_info_tmp;
      simulation_xml_handler handler(params_tmp, clone_info_tmp);
      XMLParser parser(handler);
      parser.parse(file_out);
      simulation_xml_writer(file_out, write_xml, true, params_, clone_info_);
    } else {
      simulation_xml_writer(file_out, write_xml, true, params_, clone_info_);
    }
  }
}

void task::halt() {
  halt_state();
  params_ = alps::params();
}

void task::check_parameter(bool write_xml) {
  // change in NUM_CLONES : keep old calculations and add new ones
  // change in SEED : ignored
  // change in other parameters : throw away all the old clones

  boost::filesystem::path file_in = absolute(boost::filesystem::path(file_in_str_), basedir_);
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str_), basedir_);
  if (exists(file_out)) {
    // loading parameters from *.in.xml
    alps::params params_in;
    simulation_parameters_xml_handler handler(params_in);
    XMLParser parser(handler);
    parser.parse(file_in);

    // loading from *.out.xml
    load();
    if (params_.defined("SEED")) params_in["SEED"] = params_["SEED"];

    bool changed = false;
    for (alps::params::const_iterator p = params_.begin(); p != params_.end(); ++p) {
      if (p->first != "NUM_CLONES")
        if (!params_in.defined(p->first) ||
            boost::lexical_cast<std::string>(params_in[p->first]) != boost::lexical_cast<std::string>(p->second))
          changed = true;
    }
    for (alps::params::const_iterator p = params_in.begin(); p != params_in.end(); ++p) {
      if (p->first != "NUM_CLONES")
        if (!params_.defined(p->first)) changed = true;
    }
    if (changed) {
      std::cout << "Info: parameters in " << logger::task(task_id_) << " have been changed. "
                << "All the clones are being thrown away." << std::endl;
      params_ = params_in;
      clone_info_.clear();
      clone_status_.clear();
      clone_master_.clear();
      running_.clear();
      suspended_.clear();
      finished_.clear();
    } else if (params_in.defined("NUM_CLONES")) {
      if ((params_.defined("NUM_CLONES") &&
           (int(params_in["NUM_CLONES"]) != int(params_["NUM_CLONES"]))) ||
          !params_.defined("NUM_CLONES")) {
        std::cout << "Info: number of clones in " << logger::task(task_id_) << " has been changed."
                  << std::endl;
        params_["NUM_CLONES"] = params_in["NUM_CLONES"];
        changed = true;
      }
    }

    if (changed) {
      num_clones_ = (params_["NUM_CLONES"] | 1);
      progress_ = calc_progress();
      boost::tie(weight_, dump_weight_) = calc_weight();
      status_ = calc_status();
      save(write_xml);
    }
    halt();
  }
}

void task::write_xml_archive(oxstream& os) const {
  os << alps::start_tag("SIMULATION")
     << attribute("id", task_id_+1)
     << attribute("status", task_status::to_string(status_))
     << attribute("progress", precision(progress() * 100, 3) + '%')
     << params_;
  for (unsigned int i = 0; i < clone_info_.size(); ++i) if (clone_info_[i].clone_id() == i) os << clone_info_[i];
  os << alps::end_tag("SIMULATION");
}

} // end namespace ngs_parapack
} // end namespace alps
