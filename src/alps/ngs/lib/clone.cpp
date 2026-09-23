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

#include <alps/ngs/parapack/clone.h>
// #include <alps/parapack/clone.h>
#include <alps/parapack/logger.h>
#include <boost/filesystem/operations.hpp>

#include <alps/parapack/detail/clone_workflow.hpp>

namespace alps {
namespace ngs_parapack {

//
// clone
//

clone::clone(boost::filesystem::path const& basedir, dump_policy_t dump_policy, 
  clone_timer::duration_t const& check_interval, tid_t tid, cid_t cid, alps::params const& p,
  std::string const& base, bool is_new)
  : task_id_(tid), clone_id_(cid), params_(p), basedir_(basedir), dump_policy_(dump_policy),
    timer_(check_interval) {
  params_["DIR_NAME"] = basedir_.string();
  params_["BASE_NAME"] = base;

  params_["TASK_ID"] = task_id_ + 1;
  params_["CLONE_ID"] = clone_id_ + 1;

  info_ = clone_info(clone_id_, params_, base);
  params_["WORKER_SEED"] = info_.worker_seed();
  params_["DISORDER_SEED"] = info_.disorder_seed();

  worker_ = ngs_parapack::worker_factory::make_worker(params_);
  if (!is_new) {
    // NGS stores the worker state in the HDF5 checkpoint itself.
    bool exists =
      boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile_h5()), basedir_));
    if (exists) {
      this->load();
    } else {
      std::cerr << logger::header() << "warning: dump file not found. Restarting "
                << logger::clone(task_id_, clone_id_) << std::endl;
      is_new = true;
    }
  }

  detail::clone_workflow::initialize<false>(*this, is_new,
    [this] { return worker_->fraction_completed(); }, "running");
}

clone::~clone() {}

void clone::run(boost::function<bool ()> const& stop_callback,
  boost::function<void (double)> const& progress_callback) {
  detail::clone_workflow::step<false>(*this,
    [this] { return worker_->fraction_completed(); },
    [&] { worker_->run(stop_callback, progress_callback); return false; });
}

bool clone::halted() const { return !worker_; }

clone_info const& clone::info() const { return info_; }

void clone::load() {
  detail::clone_workflow::load_ngs<false>(*this);
}

void clone::save() const{
  detail::clone_workflow::save_ngs<false>(*this);
}

void clone::load(hdf5::archive & ar) {
  ar >> make_pvp("log/alps", info_);
}

void clone::save(hdf5::archive & ar) const {
  ar << make_pvp("log/alps", info_);
}

void clone::checkpoint() {
  detail::clone_workflow::checkpoint<false>(*this);
}

void clone::suspend() {
  detail::clone_workflow::suspend<false>(*this);
}

void clone::do_halt() {
  detail::clone_workflow::halt<false>(*this);
}

void clone::output() const{
  std::cout << params_;
  // BOOST_FOREACH(ObservableSet const& m, measurements_) std::cout << m;
}

#ifdef ALPS_HAVE_MPI

//
// clone_mpi
//

clone_mpi::clone_mpi(boost::mpi::communicator const& ctrl, boost::mpi::communicator const& work,
  boost::filesystem::path const& basedir, dump_policy_t dump_policy,
  clone_timer::duration_t const& check_interval, clone_create_msg_t const& msg)
  : ctrl_(ctrl), work_(work), task_id_(msg.task_id), clone_id_(msg.clone_id),
    group_id_(msg.group_id), params_(msg.p), basedir_(basedir),
    dump_policy_(dump_policy), timer_(check_interval) {
  bool is_new = msg.is_new;
  params_["DIR_NAME"] = basedir.string();
  params_["BASE_NAME"] = msg.base;

  params_["TASK_ID"] = task_id_ + 1;
  params_["CLONE_ID"] = clone_id_ + 1;

  info_ = clone_info_mpi(work_, clone_id_, params_, msg.base);
  params_["WORKER_SEED"] = info_.worker_seed();
  params_["DISORDER_SEED"] = info_.disorder_seed();

  if (work_.size() > 1)
    worker_ = alps::ngs_parapack::parallel_worker_factory::make_worker(work_, params_);
  else
    worker_ = alps::ngs_parapack::worker_factory::make_worker(params_);
  if (!is_new) {
    bool exists = 
      boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile_h5()), basedir_));
    exists = boost::mpi::all_reduce(work_, exists, boost::mpi::bitwise_and<bool>());
    if (exists) {
      this->load();
    } else {
      if (work_.rank() == 0) {
        std::cerr << logger::header() << "warning: dump file not found. Restarting "
                  << logger::clone(task_id_, clone_id_) << std::endl;
      }
      is_new = true;
    }
  }

  detail::clone_workflow::initialize<true>(*this, is_new,
    [this] { return worker_->fraction_completed(); }, "running");
}

clone_mpi::~clone_mpi() {}

void clone_mpi::run(boost::function<bool ()> const& stop_callback,
  boost::function<void (double)> const& progress_callback) {
  detail::clone_workflow::run_mpi(*this,
    [this] { return worker_->fraction_completed(); },
    [&] { worker_->run(stop_callback, progress_callback); return false; });
}

void clone_mpi::checkpoint() {
  detail::clone_workflow::request(*this, mcmp_tag::clone_checkpoint);
}

void clone_mpi::suspend() {
  detail::clone_workflow::request(*this, mcmp_tag::clone_suspend);
}

void clone_mpi::do_checkpoint() {
  detail::clone_workflow::checkpoint<true>(*this);
}

void clone_mpi::do_suspend() {
  detail::clone_workflow::suspend<true>(*this);
}

void clone_mpi::do_halt() {
  detail::clone_workflow::halt<true>(*this);
}

bool clone_mpi::halted() const { return !worker_; }

clone_info const& clone_mpi::info() const { return info_; }

void clone_mpi::load() {
  detail::clone_workflow::load_ngs<true>(*this);
}

void clone_mpi::save() const{
  detail::clone_workflow::save_ngs<true>(*this);
}

void clone_mpi::load(hdf5::archive & ar) {
  ar >> make_pvp("log/alps", info_);
}

void clone_mpi::save(hdf5::archive & ar) const {
  ar << make_pvp("log/alps", info_);
}

void clone_mpi::output() const{
  std::cout << params_;
  // BOOST_FOREACH(ObservableSet const& m, measurements_) std::cout << m;
}

void clone_mpi::send_info(mcmp_tag_t tag) {
  if (work_.rank() == 0) {
    if (ctrl_.rank() == 0) std::cerr << "Error: sending to myself in clone_mpi::send_info()\n";
    ctrl_.send(0, tag, clone_info_msg_t(task_id_, clone_id_, group_id_, info_));
  }
}

void clone_mpi::send_halted() {
  if (work_.rank() == 0) {
    if (ctrl_.rank() == 0) std::cerr << "Error: sending to myself in clone_mpi::send_halted()\n";
    ctrl_.send(0, mcmp_tag::clone_halt, clone_halt_msg_t(task_id_, clone_id_, group_id_));
  }
}

#endif // ALPS_HAVE_MPI

} // end namespace ngs_parapack
} // end namespace alps
