/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2014 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "clone.h"
#include "logger.h"
#include <boost/filesystem/operations.hpp>

#include <alps/parapack/detail/clone_workflow.hpp>

namespace alps {

void save_observable(alps::hdf5::archive & ar, std::string const& prefix,
                     std::vector<ObservableSet> const& obs) {
  if (obs.size() == 1)
    ar[prefix] << obs[0];
  else
    for (std::size_t m = 0; m < obs.size(); ++m)
      ar[prefix + "/sections/" + boost::lexical_cast<std::string>(m)] << obs[m];
}

void save_observable(alps::hdf5::archive & ar, std::vector<ObservableSet> const& obs) {
  save_observable(ar, "simulation/results", obs);
}

void save_observable(alps::hdf5::archive & ar, cid_t cid, std::vector<ObservableSet> const& obs) {
  save_observable(ar, "simulation/realizations/" + boost::lexical_cast<std::string>(0) +
                  "/clones/" + boost::lexical_cast<std::string>(cid) + "/results", obs);
}

void save_observable(alps::hdf5::archive & ar, cid_t cid, int rank,
                     std::vector<ObservableSet> const& obs) {
  save_observable(ar, "simulation/realizations/" + boost::lexical_cast<std::string>(0) +
                  "/clones/" + boost::lexical_cast<std::string>(cid) +
                  "/workers/" + boost::lexical_cast<std::string>(rank) + "/results", obs);
}

bool load_observable(alps::hdf5::archive & ar, std::string const& prefix,
                     std::vector<ObservableSet>& obs) {
  obs.clear();
  if (ar.is_group(prefix)) {
    if (!ar.is_group(prefix + "/sections/0")) {
      obs.resize(1);
      ar[prefix] >> obs[0];
    } else {
      for (int m = 0; ; ++m) {
        std::string p = prefix + "/sections/" + boost::lexical_cast<std::string>(m);
        if (ar.is_group(p)) {
          obs.push_back(ObservableSet());
          ar[p] >> obs[m];
        } else {
          break;
        }
      }
    }
    return true;
  } else {
    return false;
  }
}

bool load_observable(alps::hdf5::archive & ar, std::vector<ObservableSet>& obs) {
  return load_observable(ar, "simulation/results", obs);
}

bool load_observable(alps::hdf5::archive & ar, std::string const& prefix, cid_t cid,
                     std::string const& suffix,  std::vector<ObservableSet>& obs) {
  if (!ar.is_group(prefix)) return false;
  std::vector<std::string> children = ar.list_children(prefix);
  if (children.size() > 1) {
    std::cerr << logger::header() << "warning: more than one clone is found in in prefix = "
              << prefix << ". Only the first one will be loaded." << std::endl;
  }
  if (children.size() > 0) {
    if (children[0] != boost::lexical_cast<std::string>(cid)) {
      std::cerr << logger::header() << "warning: try to find clone id = " << cid << ", but found "
                << children[0] << " in prefix = " << prefix << std::endl;
    }
    return load_observable(ar, prefix + "/" + children[0] + "/" + suffix, obs);
  }
  return false;
}

bool load_observable(alps::hdf5::archive & ar, cid_t cid, std::vector<ObservableSet>& obs) {
  std::string prefix = "simulation/realizations/" + boost::lexical_cast<std::string>(0) +
    "/clones";
  if (ar.is_group(prefix))
    return load_observable(ar, prefix, cid, "results", obs);
  else
    return false;
}

bool load_observable(alps::hdf5::archive & ar, cid_t cid, int rank,
                       std::vector<ObservableSet>& obs) {
  return load_observable(ar, "simulation/realizations/" + boost::lexical_cast<std::string>(0) +
                         "/clones", cid, "workers/" + boost::lexical_cast<std::string>(rank) +
                         "/results", obs);
}

bool load_observable(alps::IDump& dp, std::vector<ObservableSet>& obs) {
  Parameters ptmp;
  clone_info itmp;
  dp >> ptmp >> itmp >> obs;
  return true;
}
  

//
// clone
//

clone::clone(boost::filesystem::path const& basedir, alps::parapack::option opt, tid_t tid,
  cid_t cid, Parameters const& params, std::string const& base, bool is_new)
  : task_id_(tid), clone_id_(cid), params_(params), basedir_(basedir),
    dump_format_(opt.dump_format), dump_policy_(opt.dump_policy), timer_(opt.check_interval) {
  params_["DIR_NAME"] = basedir_.string();
  params_["BASE_NAME"] = base;

  params_["TASK_ID"] = task_id_ + 1;
  params_["CLONE_ID"] = clone_id_ + 1;

  info_ = clone_info(clone_id_, params_, base);
  params_["WORKER_SEED"] = info_.worker_seed();
  params_["DISORDER_SEED"] = info_.disorder_seed();

  worker_ = parapack::worker_factory::make_worker(params_);
  if (!is_new) {
    bool exists = 
      boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile()), basedir_)) &&
      (boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile_h5()), basedir_)) ||
       boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile_xdr()), basedir_)));
    if (exists) {
      this->load();
    } else {
      std::cerr << logger::header() << "warning: dump file not found. Restarting "
                << logger::clone(task_id_, clone_id_) << std::endl;
      is_new = true;
    }
  }
  if (is_new) worker_->init_observables(params_, measurements_);

  if (is_new && worker_->is_thermalized()) { // no thermalization steps
    BOOST_FOREACH(alps::ObservableSet& m, measurements_) { m.reset(true); }
  }

  detail::clone_workflow::initialize<false>(*this, is_new,
    [this] { return worker_->progress(); }, worker_->is_thermalized() ? "running" : "equilibrating");
}

clone::~clone() {}

void clone::run() {
  detail::clone_workflow::step<false>(*this,
    [this] { return worker_->progress(); },
    [this] { return detail::clone_workflow::classic_step(*worker_, measurements_); });
}

bool clone::halted() const { return !worker_; }

clone_info const& clone::info() const { return info_; }

void clone::load() {
  detail::clone_workflow::load_classic<false>(*this);
}

void clone::save() const{
  detail::clone_workflow::save_classic<false>(*this);
}

void clone::load(hdf5::archive & ar) {
  ar["parameters"] >> params_;
  ar["log/alps"] >> info_;
  load_observable(ar, clone_id_, measurements_);
}

void clone::save(hdf5::archive & ar) const {
  ar["parameters"] << params_;
  ar["log/alps"] << info_;
  save_observable(ar, clone_id_, measurements_);
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
  BOOST_FOREACH(ObservableSet const& m, measurements_) std::cout << m;
}

#ifdef ALPS_HAVE_MPI

//
// clone_mpi
//

clone_mpi::clone_mpi(boost::mpi::communicator const& ctrl, boost::mpi::communicator const& work,
  boost::filesystem::path const& basedir, alps::parapack::option opt, clone_create_msg_t const& msg)
  : ctrl_(ctrl), work_(work), task_id_(msg.task_id), clone_id_(msg.clone_id),
    group_id_(msg.group_id), params_(msg.params), basedir_(basedir), dump_format_(opt.dump_format),
    dump_policy_(opt.dump_policy), timer_(opt.check_interval) {
  bool is_new = msg.is_new;
  params_["DIR_NAME"] = basedir;
  params_["BASE_NAME"] = msg.base;

  params_["TASK_ID"] = task_id_ + 1;
  params_["CLONE_ID"] = clone_id_ + 1;

  info_ = clone_info_mpi(work_, clone_id_, params_, msg.base);
  params_["WORKER_SEED"] = info_.worker_seed();
  params_["DISORDER_SEED"] = info_.disorder_seed();

  if (work_.size() > 1)
    worker_ = alps::parapack::parallel_worker_factory::make_worker(work_, params_);
  else
    worker_ = alps::parapack::worker_factory::make_worker(params_);
  if (!is_new) {
    bool exists = 
      boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile()), basedir_)) &&
      (boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile_h5()), basedir_)) ||
       boost::filesystem::exists(absolute(boost::filesystem::path(info_.dumpfile_xdr()), basedir_)));
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
  if (is_new) worker_->init_observables(params_, measurements_);

  if (is_new && worker_->is_thermalized()) { // no thermalization steps
    BOOST_FOREACH(alps::ObservableSet& m, measurements_) { m.reset(true); }
  }

  detail::clone_workflow::initialize<true>(*this, is_new,
    [this] { return worker_->progress(); }, worker_->is_thermalized() ? "running" : "equilibrating");
}

clone_mpi::~clone_mpi() {}

void clone_mpi::run() {
  detail::clone_workflow::run_mpi(*this,
    [this] { return worker_->progress(); },
    [this] { return detail::clone_workflow::classic_step(*worker_, measurements_); });
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
  detail::clone_workflow::load_classic<true>(*this);
}

void clone_mpi::save() const{
  detail::clone_workflow::save_classic<true>(*this);
}

void clone_mpi::load(hdf5::archive & ar) {
  ar["parameters"] >> params_;
  ar["log/alps"] >> info_;
  if (work_.size() == 1)
    load_observable(ar, clone_id_, measurements_);
  else
    load_observable(ar, clone_id_, work_.rank(), measurements_);
}

void clone_mpi::save(hdf5::archive & ar) const {
  ar["parameters"] << params_;
  ar["log/alps"] << info_;
  if (work_.size() == 1)
    save_observable(ar, clone_id_, measurements_);
  else
    save_observable(ar, clone_id_, work_.rank(), measurements_);
}

void clone_mpi::output() const{
  std::cout << params_;
  BOOST_FOREACH(ObservableSet const& m, measurements_) std::cout << m;
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

} // end namespace alps
