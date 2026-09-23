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

#include "job.h"
#include "clone.h"
#include "filelock.h"
#include "measurement.h"
#include "simulation_p.h"

#include <alps/hdf5.hpp>

// some file (probably a python header) defines a tolower macro ...
#undef tolower
#undef toupper

#include <boost/config.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace alps {

//
// task
//

void task::load() {
  if (on_memory()) boost::throw_exception(std::logic_error("task::load() task already loaded"));
  params_.clear();
  obs_.clear();
  clone_info_.clear();
  boost::filesystem::path file_in = absolute(boost::filesystem::path(file_in_str_), basedir_);
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str_), basedir_);
  simulation_xml_handler handler(params_, obs_, clone_info_);
  XMLParser parser(handler);
  if (!exists(file_out)) {
    parser.parse(file_in);
  } else {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    parser.parse(file_out);
  }

  num_clones_ = range_type(params_.value_or_default("NUM_CLONES", "1"), params_);
  restore_state();
}

void task::save(alps::parapack::option const& opt) const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str_), basedir_);
  {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    if (exists(file_out)) {
      Parameters params_tmp;
      std::vector<ObservableSet> obs_tmp;
      std::deque<clone_info> clone_info_tmp;
      simulation_xml_handler handler(params_tmp, obs_tmp, clone_info_tmp);
      XMLParser parser(handler);
      parser.parse(file_out);
      simulation_xml_writer(file_out, opt.write_xml, true, params_, obs_tmp, clone_info_);
    } else {
      simulation_xml_writer(file_out, opt.write_xml, true, params_, obs_, clone_info_);
    }
  }
}

void task::save_observable(alps::parapack::option const& opt) const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str_), basedir_);
  {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    if (exists(file_out)) {
      Parameters params_tmp;
      std::vector<ObservableSet> obs_tmp;
      std::deque<clone_info> clone_info_tmp;
      simulation_xml_handler handler(params_tmp, obs_tmp, clone_info_tmp);
      XMLParser parser(handler);
      parser.parse(file_out);
      if (obs_.size() > 1 && !params_tmp.defined("NUM_REPLICAS"))
        params_tmp["NUM_REPLICAS"] = obs_.size();
      simulation_xml_writer(file_out, opt.write_xml, true, params_tmp, obs_, clone_info_tmp);
      if (obs_.size() == 1) {
        if (opt.dump_format == dump_format::hdf5) {
          #pragma omp critical (hdf5io)
          {
            boost::filesystem::path file = absolute(boost::filesystem::path(base_ + ".out.h5"),
              basedir_);
            hdf5::archive h5(file.string(), "a");
            h5["/parameters"] << params_tmp;
            h5["/simulation/results"] << obs_[0];
            // for (std::size_t n = 0; n < oss.size(); ++n)
            //   h5["/simulation/realizations/0/clones/" +
            //      boost::lexical_cast<std::string>(n) + "/results"] << oss[n][0];
          }
        } else {
          boost::filesystem::path file = absolute(boost::filesystem::path(base_ + ".out.xdr"),
            basedir_);
          OXDRFileDump dp(file);
          dp << params_tmp << obs_[0];
        }
      } else {
        for (std::size_t i = 0; i < obs_.size(); ++i) {
          Parameters p = params_tmp;
          if (!p.defined("REPLICA")) p["REPLICA"] = i+1;
          if (opt.dump_format == dump_format::hdf5) {
            #pragma omp critical (hdf5io)
            {
              boost::filesystem::path file = absolute(boost::filesystem::path(
                base_ + ".replica" + boost::lexical_cast<std::string>(i+1) + ".h5"), basedir_);
              hdf5::archive h5(file.string(), "a");
              h5["/parameters"] << p;
              h5["/simulation/results"] << obs_[i];
              // for (std::size_t n = 0; n < oss.size(); ++n)
              //   h5["/simulation/realizations/0/clones/" +
              //      boost::lexical_cast<std::string>(n) + "/results"] << oss[n][i];
            }
          } else {
            boost::filesystem::path file = absolute(boost::filesystem::path(
              base_ + ".replica" + boost::lexical_cast<std::string>(i+1) + ".xdr"), basedir_);
            OXDRFileDump dp(file);
            dp << p << obs_[i];
          }
        }
      }
    } else {
      boost::throw_exception(std::logic_error("task::save_observable()"));
    }
  }
}

void task::halt() {
  halt_state();
  params_.clear();
  obs_.clear();
}

void task::check_parameter(alps::parapack::option const& opt) {
  // change in NUM_CLONES : keep old calculations and add new ones
  // change in SEED : ignored
  // change in other parameters : throw away all the old clones

  boost::filesystem::path file_in = absolute(boost::filesystem::path(file_in_str_), basedir_);
  boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str_), basedir_);
  if (exists(file_out)) {
    // loading parameters from *.in.xml
    Parameters params_in;
    simulation_parameters_xml_handler handler(params_in);
    XMLParser parser(handler);
    parser.parse(file_in);

    // loading from *.out.xml
    load();
    if (params_.defined("SEED")) params_in["SEED"] = params_["SEED"];

    bool changed = false;
    BOOST_FOREACH(Parameter const& p, params_) {
      if (p.key() != "NUM_CLONES")
        if (!params_in.defined(p.key()) || params_in[p.key()] != p.value()) changed = true;
    }
    BOOST_FOREACH(Parameter const& p, params_in) {
      if (p.key() != "NUM_CLONES")
        if (!params_.defined(p.key())) changed = true;
    }
    if (changed) {
      std::cout << "Info: parameters in " << logger::task(task_id_) << " have been changed. "
                << "All the clones are being thrown away." << std::endl;
      params_ = params_in;
      obs_.clear();
      clone_info_.clear();
      clone_status_.clear();
      clone_master_.clear();
      running_.clear();
      suspended_.clear();
      finished_.clear();
    } else if (params_in.defined("NUM_CLONES")) {
      if ((params_.defined("NUM_CLONES") && params_in["NUM_CLONES"] != params_["NUM_CLONES"]) ||
          !params_.defined("NUM_CLONES")) {
        std::cout << "Info: number of clones in " << logger::task(task_id_) << " has been changed."
                  << std::endl;
        params_["NUM_CLONES"] = params_in["NUM_CLONES"];
        changed = true;
      }
    }

    if (changed) {
      num_clones_ = range_type(params_.value_or_default("NUM_CLONES", "1"), params_);
      progress_ = calc_progress();
      boost::tie(weight_, dump_weight_) = calc_weight();
      status_ = calc_status();
      save(opt);
    }
    halt();
  }
}

void task::evaluate(alps::parapack::option const& opt) {
  std::cout << "evaluating " << file_out_str() << std::endl;

  if (!on_memory()) load();

  // bool same_weight = params_.defined("EVALUATE_CLONES_WITH_SAME_WEIGHT") &&
  //   static_cast<bool>(alps::evaluate("EVALUATE_CLONES_WITH_SAME_WEIGHT", params_));
  bool only_finished = params_.defined("EVALUATE_ONLY_FINISHED_CLONES") &&
    static_cast<bool>(alps::evaluate("EVALUATE_ONLY_FINISHED_CLONES", params_));

  obs_.clear();
  std::vector<cid_t> clones;
  BOOST_FOREACH(cid_t cid, finished_) clones.push_back(cid);
  if (!only_finished) {
    BOOST_FOREACH(cid_t cid, running_) clones.push_back(cid);
    BOOST_FOREACH(cid_t cid, suspended_) clones.push_back(cid);
  }
  std::sort(clones.begin(), clones.end());

  alps::Parameters p = params_;
  p["DIR_NAME"] = basedir_.string();
  p["BASE_NAME"] = base_;

  boost::shared_ptr<parapack::abstract_evaluator> evaluator
    = parapack::evaluator_factory::instance()->make_evaluator(p);

  std::cout << "  loading clones: ";
  // std::vector<std::vector<ObservableSet> > oss;
  BOOST_FOREACH(cid_t cid, clones) {
    std::cout << (cid+1) << ' ' << std::flush;
    // std::vector<ObservableSet> os;
    if (clone_info_[cid].checkpoints().size() == 1) {
      boost::filesystem::path dump_h5 =
        absolute(boost::filesystem::path(clone_info_[cid].checkpoints()[0] + ".h5"), basedir_);
      boost::filesystem::path dump_xdr =
        absolute(boost::filesystem::path(clone_info_[cid].checkpoints()[0] + ".xdr"), basedir_);
      std::vector<ObservableSet> o;
      bool success;
      if (exists(dump_h5)) {
        #pragma omp critical (hdf5io)
        {
          hdf5::archive h5(dump_h5);
          success = load_observable(h5, cid, o);
        }
      } else {
        IXDRFileDump dp(dump_xdr);
        success = load_observable(dp, o);
      }
      if (success) {
        // evaluator->load(o, os);
        evaluator->load(o, obs_);
      }
    } else {
      for (unsigned int w = 0; w < clone_info_[cid].checkpoints().size(); ++w) {
        boost::filesystem::path dump_h5 =
          absolute(boost::filesystem::path(clone_info_[cid].checkpoints()[w] + ".h5"), basedir_);
        boost::filesystem::path dump_xdr =
          absolute(boost::filesystem::path(clone_info_[cid].checkpoints()[w] + ".xdr"), basedir_);
        std::vector<ObservableSet> o;
        bool success;
        if (exists(dump_h5)) {
          #pragma omp critical (hdf5io)
          {
            hdf5::archive h5(dump_h5);
            success = load_observable(h5, cid, w, o);
          }
        } else {
          IXDRFileDump dp(dump_xdr);
          success = load_observable(dp, o);
        }
        if (success) {
          // evaluator->load(o, os);
          evaluator->load(o, obs_);
        }
      }
    }
    // oss.push_back(os);
  }
  std::cout << std::endl;
  if (clones.size() > 0) {
    evaluator->evaluate(obs_);
    save_observable(opt);
  }
  halt();
}

void task::write_xml_archive(oxstream& os) const {
  os << alps::start_tag("SIMULATION")
     << attribute("id", task_id_+1)
     << attribute("status", task_status::to_string(status_))
     << attribute("progress", precision(progress() * 100, 3) + '%')
     << params_;
  if (obs_.size() == 1) {
    obs_[0].write_xml(os);
  } else {
    for (unsigned int i = 0; i < obs_.size(); ++i) obs_[i].write_xml_with_id(os, i+1);
  }
  for (unsigned int i = 0; i < clone_info_.size(); ++i) if (clone_info_[i].clone_id() == i) os << clone_info_[i];
  os << alps::end_tag("SIMULATION");
}

} // end namespace alps
