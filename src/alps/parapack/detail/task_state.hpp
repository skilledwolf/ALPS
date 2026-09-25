// Copyright (C) 1997-2013 by Synge Todo <wistaria@comp-phys.org>
// SPDX-License-Identifier: MIT
// Shared task bookkeeping for the classic and NGS schedulers.
#ifndef ALPS_PARAPACK_DETAIL_TASK_STATE_HPP
#define ALPS_PARAPACK_DETAIL_TASK_STATE_HPP

#include <alps/parapack/clone_info.h>
#include <alps/parapack/integer_range.h>
#include <alps/parapack/logger.h>
#include <alps/parapack/util.h>
#include <alps/parser/xmlstream.h>
#include <boost/optional.hpp>
#include <boost/regex.hpp>
#include <deque>
#include <set>

namespace alps { namespace parapack { namespace detail {
template<class Task> class job_task_xml_handler;
template<class Task> class job_tasks_xml_handler;

template<class Derived, class Parameters>
class task_state {
public:
  typedef integer_range<uint32_t> range_type;

  task_state() : status_(task_status::Undefined) {}
  task_state(boost::filesystem::path const& file) : status_(task_status::Undefined) {
    basedir_ = file.parent_path();
    file_in_str_ = file.filename().string();
    file_out_str_ = regex_replace(file.filename().string(), boost::regex("\\.in\\.xml$"), ".out.xml");
    if (file_in_str_ == file_out_str_) {
      file_in_str_ = regex_replace(file.filename().string(), boost::regex("\\.out\\.xml$"), ".in.xml");
      file_out_str_ = file.filename().string();
    }
    base_ = regex_replace(file_out_str_, boost::regex("\\.out\\.xml$"), "");
  }

  task_status_t status() const { return status_; }
  tid_t task_id() const { return task_id_; }
  bool on_memory() const {
    return status_ == task_status::Ready || status_ == task_status::Running ||
      status_ == task_status::Continuing || status_ == task_status::Idling;
  }

  range_type const& num_clones() const {
    if (!on_memory()) boost::throw_exception(std::logic_error(
      "task::num_clones() task not loaded"));
    return num_clones_;
  }
  uint32_t num_running() const {
    if (!on_memory()) boost::throw_exception(std::logic_error(
      "task::num_running() task not loaded"));
    return running_.size();
  }
  uint32_t num_suspended() const {
    if (!on_memory()) boost::throw_exception(std::logic_error(
      "task::num_suspended() task not loaded"));
    return suspended_.size();
  }
  uint32_t num_finished() const {
    if (!on_memory()) boost::throw_exception(std::logic_error(
      "task::num_finished() task not loaded"));
    return finished_.size();
  }
  uint32_t num_started() const {
    return num_running() + num_suspended() + num_finished();
  }
  double progress() const { return progress_; }
  double weight() const { return weight_; }

  bool can_dispatch() const {
    return num_suspended() > 0 || num_started() < num_clones_.max BOOST_PREVENT_MACRO_SUBSTITUTION ();
  }

  template<typename PROXY, typename GROUP>
  boost::optional<cid_t>
  dispatch_clone(PROXY& proxy, GROUP const& group) {
    if (!on_memory()) derived().load();
    if (!can_dispatch()) {
      if (num_running() == 0) derived().halt();
      return boost::optional<cid_t>();
    }
    bool is_new;
    cid_t cid;
    Process master = group.master();
    if (num_suspended()) {
      is_new = false;
      cid = *suspended_.begin();
      suspended_.erase(suspended_.begin());
      clone_status_[cid] = clone_status::Running;
      clone_master_[cid] = master;
    } else {
      is_new = true;
      cid = clone_status_.size();
      clone_status_.push_back(clone_status::Running);
      clone_master_.push_back(master);
      clone_info_.push_back(clone_info(cid));
    }
    std::cout << logger::header() << (is_new ? "dispatching a new " : "resuming a suspended ")
              << logger::clone(task_id_, cid) << " on " << logger::group(group) << std::endl;
    running_.insert(cid);
    proxy.start(task_id_, cid, group, params_, base_, is_new);
    boost::tie(weight_, dump_weight_) = calc_weight();
    status_ = calc_status();
    if (!is_new) report(proxy, cid);
    return boost::optional<cid_t>(cid);
  }

  bool is_running(cid_t cid) const {
    return clone_status_[cid] == clone_status::Running;
  }

  template<typename PROXY>
  void checkpoint(PROXY& proxy, cid_t cid) {
    if (clone_status_[cid] == clone_status::Running) {
      proxy.checkpoint(clone_master_[cid]);
      if (proxy.is_local(clone_master_[cid])) {
        clone_info const& info = proxy.info(clone_master_[cid]);
        std::cout << logger::header() << "regular checkpoint: "
                  << logger::clone(task_id_, cid) << " is " << info.phase()
                  << " (" << precision(info.progress() * 100, 3) << "% done)\n";
        info_updated(cid, info);
        // save(write_xml);
      }
    }
  }

  template<typename PROXY>
  void report(PROXY& proxy, cid_t cid) {
    if (clone_status_[cid] == clone_status::Running) {
      proxy.update_info(clone_master_[cid]);
      if (proxy.is_local(clone_master_[cid])) {
        clone_info const& info = proxy.info(clone_master_[cid]);
        std::cout << logger::header() << "progress report: "
                  << logger::clone(task_id_, cid) << " is " << info.phase()
                  << " (" << precision(info.progress() * 100, 3) << "% done)\n";
      }
    }
  }

  template<typename PROXY, typename Options>
  void suspend_remote_clones(PROXY& proxy, Options const& opt) {
    if (on_memory()) {
      BOOST_FOREACH(cid_t cid, running_) {
        if (clone_status_[cid] == clone_status::Running && !proxy.is_local(clone_master_[cid])) {
          clone_status_[cid] = clone_status::Stopping;
          proxy.suspend(clone_master_[cid]);
        }
      }
      if (num_running() == 0) {
        derived().save(opt);
        derived().halt();
      }
    }
  }

  template<typename PROXY, typename GROUP, typename Options>
  void suspend_clone(PROXY& proxy, Options const& opt, cid_t cid, GROUP const& g) {
    if (clone_status_[cid] == clone_status::Running && proxy.is_local(clone_master_[cid])) {
      clone_status_[cid] = clone_status::Stopping;
      proxy.suspend(clone_master_[cid]);
      clone_info const& info = proxy.info(clone_master_[cid]);
      clone_suspended(cid, g, info);
      proxy.destroy(clone_master_[cid]);
      if (num_running() == 0) {
        derived().save(opt);
        derived().halt();
      }
    }
  }

  template<typename PROXY, typename GROUP, typename Options>
  void halt_clone(PROXY& proxy, Options const& opt, cid_t cid, GROUP const& g) {
    if (clone_status_[cid] == clone_status::Idling) {
      std::cout << logger::header() << logger::clone(task_id_, cid) << " finished"
                << " on " << logger::group(g) << std::endl;
      clone_status_[cid] = clone_status::Stopping;
      proxy.halt(clone_master_[cid]);
      if (proxy.is_local(clone_master_[cid])) {
        clone_info const& info = proxy.info(clone_master_[cid]);
        clone_halted(cid, info);
      }
      proxy.destroy(clone_master_[cid]);
      if (num_running() == 0) {
        derived().save(opt);
        derived().halt();
      }
    }
  }

  void info_updated(cid_t cid, clone_info const& info) {
    if (clone_status_[cid] == clone_status::Running) {
      clone_info_[cid] = info;
      if (info.progress() >= 1) clone_status_[cid] = clone_status::Idling;
    }
  }
  template<typename GROUP>
  void clone_suspended(cid_t cid, GROUP const& g, clone_info const& info) {
    if (clone_status_[cid] != clone_status::Stopping)
      boost::throw_exception(std::logic_error("clone is not stopping"));
    std::cout << logger::header() << logger::clone(task_id_, cid)
              << " suspended (" << precision(info.progress() * 100, 3)
              << "% done)" << " on " << logger::group(g) << std::endl;
    clone_info_[cid] = info;
    clone_status_[cid] = clone_status::Suspended;
    running_.erase(cid);
    suspended_.insert(cid);
    boost::tie(weight_, dump_weight_) = calc_weight();
  }
  void clone_halted(cid_t cid, clone_info const& info) {
    info_updated(cid, info);
    clone_halted(cid);
  }
  void clone_halted(cid_t cid) {
    if (clone_status_[cid] != clone_status::Stopping)
      boost::throw_exception(std::logic_error("clone is not stopping"));
    clone_status_[cid] = clone_status::Finished;
    running_.erase(cid);
    finished_.insert(cid);
    progress_ = calc_progress();
    status_ = calc_status();
    boost::tie(weight_, dump_weight_) = calc_weight();

    derived().clear_observable_cache();
  }

  void write_xml_summary(oxstream& os) const {
    os << start_tag("TASK")
       << attribute("id", task_id_+1)
       << attribute("status", task_status::to_string(status_))
       << attribute("progress", precision(progress() * 100, 3) + '%')
       << attribute("weight", precision(dump_weight_, 3))
       << start_tag("INPUT")
       << attribute("file", file_in_str_)
       << end_tag("INPUT")
       << start_tag("OUTPUT")
       << attribute("file", file_out_str_)
       << end_tag("OUTPUT")
       << end_tag("TASK");
  }

  std::string const& file_in_str() const { return file_in_str_; }
  std::string const& file_out_str() const { return file_out_str_; }

protected:
  void restore_state() {
    clone_status_.clear();
    clone_master_.clear();
    running_.clear();
    suspended_.clear();
    finished_.clear();
    BOOST_FOREACH(clone_info const& info, clone_info_) {
      if (info.progress() < 1){
        clone_status_.push_back(clone_status::Suspended);
        suspended_.insert(info.clone_id());
      } else {
        clone_status_.push_back(clone_status::Finished);
        finished_.insert(info.clone_id());
      }
      clone_master_.push_back(Process());
    }

    status_ = task_status::Ready;
    progress_ = calc_progress();
    boost::tie(weight_, dump_weight_) = calc_weight();
    status_ = calc_status();
  }

  void halt_state() {
    if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
    if (running_.size()) boost::throw_exception(std::logic_error("running clone exists"));

    switch (status_) {
    case task_status::Ready:
      status_ = task_status::NotStarted;
      break;
    case task_status::Running:
      status_ = task_status::Suspended;
      break;
    case task_status::Continuing:
      status_ = task_status::Finished;
      break;
    case task_status::Idling:
      status_ = task_status::Completed;
      break;
    default:
      boost::throw_exception(std::logic_error("unknown task_status"));
    }

    clone_status_.clear();
    clone_master_.clear();
    clone_info_.clear();
    running_.clear();
    suspended_.clear();
    finished_.clear();
  }

  double calc_progress() const {
    if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
    return (double)(num_finished()) / num_clones().min BOOST_PREVENT_MACRO_SUBSTITUTION ();
  }
  std::pair<double, double> calc_weight() const {
    if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
    double w;
    if (num_suspended() > 0)
      w = 4.0;
    else if (num_started() == 0)
      // NotStarted
      w = 3.0;
    else if (num_started() < num_clones().min BOOST_PREVENT_MACRO_SUBSTITUTION ())
      // Running
      w = 2.0 - (double)(num_started()) / num_clones().min BOOST_PREVENT_MACRO_SUBSTITUTION ();
    else
      // Continuing
      w = 1.0 - (double)(num_started()) / num_clones().max BOOST_PREVENT_MACRO_SUBSTITUTION ();
    double d = (num_running() + num_suspended() > 0) ? 4.0 : w;
    return std::make_pair(w, d);
  }
  task_status_t calc_status() const {
    if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
    if (num_started() == 0)
      return task_status::Ready;
    else if (num_finished() < num_clones_.min BOOST_PREVENT_MACRO_SUBSTITUTION ())
      return task_status::Running;
    else if (num_finished() < num_clones_.max BOOST_PREVENT_MACRO_SUBSTITUTION ())
      return task_status::Continuing;
    else
      return task_status::Idling;
  }

  template<class> friend class job_task_xml_handler;
  template<class> friend class job_tasks_xml_handler;
  Derived& derived() { return static_cast<Derived&>(*this); }

  // these parameters will be read from job xml file
  task_status_t status_;
  tid_t task_id_;
  double progress_;
  double weight_;
  double dump_weight_;
  boost::filesystem::path basedir_;
  std::string file_in_str_;
  std::string file_out_str_;
  std::string base_;

  // these parameters will be read from simulation xml file
  Parameters params_;

  range_type num_clones_;
  std::deque<clone_info> clone_info_;
  std::deque<clone_status_t> clone_status_;
  std::deque<Process> clone_master_;
  std::set<cid_t> running_;   // id of running clones
  std::set<cid_t> suspended_; // id of suspended clones
  std::set<cid_t> finished_;  // id of finished clones
};

} } }
#endif
