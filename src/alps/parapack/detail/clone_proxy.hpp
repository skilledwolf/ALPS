// Copyright (C) 1997-2012 by Synge Todo and the ALPS contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include <alps/parapack/clone_info.h>
#include <alps/parapack/process.h>
#include <boost/filesystem/path.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <stdexcept>
#include <tuple>

namespace alps { namespace detail {

// Both worker runtimes use the same local/remote proxy protocol. Only the
// concrete clone, parameter type, and constructor options differ.
template<class Clone, class Parameters, class... Options>
class clone_proxy {
public:
  clone_proxy(Clone*& clone_ptr, boost::filesystem::path const& basedir, Options... options)
    : clone_ptr_(clone_ptr), basedir_(basedir), options_(options...) {}

  bool is_local(Process const&) const { return true; }
  void start(tid_t tid, cid_t cid, thread_group const&, Parameters const& params,
    std::string const& base, bool is_new) {
    clone_ptr_ = std::apply([&](auto const&... options) {
      return new Clone(basedir_, options..., tid, cid, params, base, is_new);
    }, options_);
  }
  clone_info const& info(Process const&) const {
    if (!clone_ptr_)
      boost::throw_exception(std::logic_error("clone_proxy::info()"));
    return clone_ptr_->info();
  }
  void checkpoint(Process const&) { if (clone_ptr_) clone_ptr_->checkpoint(); }
  void update_info(Process const&) const {}
  void suspend(Process const&) { if (clone_ptr_) clone_ptr_->suspend(); }
  void halt(Process const&) {}
  void destroy(Process const&) {
    delete clone_ptr_;
    clone_ptr_ = nullptr;
  }

private:
  Clone*& clone_ptr_;
  boost::filesystem::path basedir_;
  std::tuple<Options...> options_;
};

#ifdef ALPS_HAVE_MPI
template<class Clone, class Parameters, class CreateMessage, class... Options>
class clone_proxy_mpi {
public:
  clone_proxy_mpi(Clone*& clone_ptr, boost::mpi::communicator const& comm_ctrl,
    boost::mpi::communicator const& comm_work, boost::filesystem::path const& basedir,
    Options... options)
    : clone_ptr_(clone_ptr), comm_ctrl_(comm_ctrl), comm_work_(comm_work), basedir_(basedir),
      options_(options...) {}

  bool is_local(Process const& proc) const { return proc == 0; }
  void start(tid_t tid, cid_t cid, process_group const& procs, Parameters const& params,
    std::string const& base, bool is_new) const {
    CreateMessage msg(tid, cid, procs.group_id, params, base, is_new);
    bool worker_on_master = false;
    for (Process p : procs.process_list) {
      if (p == 0) worker_on_master = true;
      else comm_ctrl_.send(p, mcmp_tag::clone_create, msg);
    }
    if (worker_on_master)
      clone_ptr_ = std::apply([&](auto const&... options) {
        return new Clone(comm_ctrl_, comm_work_, basedir_, options..., msg);
      }, options_);
  }

  clone_info const& info(Process const& proc) const {
    if (proc != 0 || !clone_ptr_) {
      std::cerr << "clone_proxy_mpi::info()\n";
      boost::throw_exception(std::logic_error("clone_proxy_mpi::info()"));
    }
    return clone_ptr_->info();
  }
  void checkpoint(Process const& proc) {
    if (proc == 0) {
      if (clone_ptr_) clone_ptr_->checkpoint();
    } else comm_ctrl_.send(proc, mcmp_tag::clone_checkpoint);
  }
  void update_info(Process const& proc) const {
    if (proc != 0) comm_ctrl_.send(proc, mcmp_tag::clone_info);
  }
  void suspend(Process const& proc) {
    if (proc == 0) {
      if (clone_ptr_) clone_ptr_->suspend();
    } else comm_ctrl_.send(proc, mcmp_tag::clone_suspend);
  }
  void halt(Process const& proc) {
    if (proc != 0) comm_ctrl_.send(proc, mcmp_tag::clone_halt);
  }
  void destroy(Process const& proc) {
    if (proc == 0) {
      delete clone_ptr_;
      clone_ptr_ = nullptr;
    }
  }

private:
  Clone*& clone_ptr_;
  boost::mpi::communicator comm_ctrl_, comm_work_;
  boost::filesystem::path basedir_;
  std::tuple<Options...> options_;
};
#endif

}} // namespace alps::detail
