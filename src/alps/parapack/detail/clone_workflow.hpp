// Copyright (C) 1997-2014 by Synge Todo <wistaria@comp-phys.org>
// SPDX-License-Identifier: MIT
#pragma once

#include <alps/parapack/clone_timer.h>
#include <alps/parapack/types.h>
#include <alps/hdf5.hpp>
#include <alps/osiris.h>
#include <boost/filesystem/operations.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <stdexcept>

namespace alps { namespace detail {

// Runtime adapters provide worker progress and one update; the workflow owns
// lifecycle transitions and the ordering of MPI control messages/collectives.
struct clone_workflow {
  template<bool MPI, class Clone> static bool worker_dump(Clone const& clone) {
    bool enabled = clone.dump_policy_ == dump_policy::All ||
      (clone.dump_policy_ == dump_policy::RunningOnly && clone.info_.progress() < 1);
    if constexpr (MPI) broadcast(clone.work_, enabled, 0);
    return enabled;
  }

  template<bool MPI, class Clone> static void load_classic(Clone& clone) {
    auto worker_path = absolute(boost::filesystem::path(clone.info_.dumpfile()), clone.basedir_);
    auto h5_path = absolute(boost::filesystem::path(clone.info_.dumpfile_h5()), clone.basedir_);
    if (exists(h5_path)) {
      #pragma omp critical (hdf5io)
      {
        hdf5::archive ar(h5_path.string());
        ar["/"] >> clone;
      }
    } else {
      IXDRFileDump dump(absolute(boost::filesystem::path(clone.info_.dumpfile_xdr()), clone.basedir_));
      dump >> clone.params_ >> clone.info_ >> clone.measurements_;
    }
    if (worker_dump<MPI>(clone)) {
      if constexpr (MPI)
        worker_path = absolute(boost::filesystem::path(clone.info_.dumpfile()), clone.basedir_);
      IXDRFileDump dump(worker_path);
      clone.worker_->load_worker(dump);
    }
  }

  template<bool MPI, class Clone> static void save_classic(Clone const& clone) {
    if (clone.dump_format_ == dump_format::hdf5) {
      #pragma omp critical (hdf5io)
      {
        auto path = absolute(boost::filesystem::path(clone.info_.dumpfile_h5()), clone.basedir_);
        hdf5::archive ar(path.string(), "a");
        ar["/"] << clone;
      }
    } else if (clone.dump_format_ == dump_format::xdr) {
      OXDRFileDump dump(absolute(boost::filesystem::path(clone.info_.dumpfile_xdr()), clone.basedir_));
      dump << clone.params_ << clone.info_ << clone.measurements_;
    }
    auto path = absolute(boost::filesystem::path(clone.info_.dumpfile()), clone.basedir_);
    if (worker_dump<MPI>(clone)) {
      OXDRFileDump dump(path);
      clone.worker_->save_worker(dump);
    } else if (exists(path)) remove(path);
  }

  template<bool MPI, class Clone> static void load_ngs(Clone& clone) {
    auto path = absolute(boost::filesystem::path(clone.info_.dumpfile_h5()), clone.basedir_);
    bool enabled = worker_dump<MPI>(clone);
    #pragma omp critical (hdf5io)
    {
      hdf5::archive ar(path.string());
      ar["/"] >> clone;
      if (enabled) clone.worker_->load_worker(ar);
    }
  }

  template<bool MPI, class Clone> static void save_ngs(Clone const& clone) {
    auto path = absolute(boost::filesystem::path(clone.info_.dumpfile_h5()), clone.basedir_);
    bool enabled = worker_dump<MPI>(clone);
    #pragma omp critical (hdf5io)
    {
      hdf5::archive ar(path.string(), "a");
      ar["/"] << clone;
      if (enabled) clone.worker_->save_worker(ar);
    }
  }

  template<bool MPI, class Clone> static bool leader(Clone const& clone) {
    if constexpr (MPI) return clone.work_.rank() == 0;
    else return true;
  }

  template<bool MPI, class Clone> static void complete(Clone& clone) {
    if constexpr (MPI) {
      if (clone.group_id_ != 0) {
        clone.send_info(mcmp_tag::clone_info);
        return;
      }
    }
    clone.do_halt();
  }

  template<bool MPI, class Clone, class Progress>
  static void initialize(Clone& clone, bool is_new, Progress progress, char const* phase) {
    double previous = [&] {
      if constexpr (MPI) return clone.info_.progress();
      else return progress();
    }();
    if (is_new || previous < 1) clone.info_.start(phase);
    bool finished = leader<MPI>(clone) && is_new && progress() >= 1;
    // Halting a local MPI group saves and synchronizes every worker, including
    // when a newly constructed worker is already complete.
    if constexpr (MPI) broadcast(clone.work_, finished, 0);
    if (finished) {
      if (leader<MPI>(clone)) {
        clone.info_.set_progress(progress());
        clone.info_.stop();
      }
      complete<MPI>(clone);
    }
    if (leader<MPI>(clone) && !is_new) clone.timer_.reset(progress());
    clone.loops_ = 1;
  }

  template<class Worker, class Measurements>
  static bool classic_step(Worker& worker, Measurements& measurements) {
    bool thermalized = worker.is_thermalized();
    worker.run(measurements);
    if (!thermalized && worker.is_thermalized()) {
      for (auto& measurement : measurements) measurement.reset(true);
      return true;
    }
    return false;
  }

  template<bool MPI, class Clone, class Progress, class Step>
  static void step(Clone& clone, Progress progress, Step advance) {
    for (clone_timer::loops_t i = 0; i < clone.loops_; ++i) {
      double previous = progress();
      if (advance() && leader<MPI>(clone)) {
        clone.info_.stop();
        clone.info_.start("running");
      }
      if (previous < 1 && progress() >= 1) {
        if (leader<MPI>(clone)) {
          clone.info_.set_progress(progress());
          clone.info_.stop();
        }
        complete<MPI>(clone);
        return;
      }
    }
    if (leader<MPI>(clone)) {
      clone.info_.set_progress(progress());
      clone.loops_ = clone.timer_.next_loops(clone.loops_);
    }
    if constexpr (MPI) broadcast(clone.work_, clone.loops_, 0);
  }

  template<bool MPI, class Clone> static void checkpoint(Clone& clone) {
    if (leader<MPI>(clone) && clone.info_.progress() < 1) clone.info_.stop();
    clone.save();
  }

  template<bool MPI, class Clone> static void suspend(Clone& clone) {
    if (leader<MPI>(clone)) clone.info_.stop();
    clone.save();
    if constexpr (MPI) clone.work_.barrier();
    clone.worker_.reset();
  }

  template<bool MPI, class Clone> static void halt(Clone& clone) {
    if (leader<MPI>(clone) && clone.info_.progress() < 1) {
      if constexpr (MPI) std::cerr << "clone is not finished\n";
      boost::throw_exception(std::logic_error("clone is not finished"));
    }
    clone.save();
    if constexpr (MPI) clone.work_.barrier();
    clone.worker_.reset();
  }

#ifdef ALPS_HAVE_MPI
  template<class Clone, class Progress, class Step>
  static void run_mpi(Clone& clone, Progress progress, Step advance) {
    int tag = clone.info_.progress() < 1 ? mcmp_tag::do_step : mcmp_tag::do_nothing;
    if (leader<true>(clone) && clone.ctrl_.iprobe(0, boost::mpi::any_tag)) {
      boost::mpi::status status = clone.ctrl_.recv(0, boost::mpi::any_tag);
      tag = status.tag();
      if (tag == mcmp_tag::clone_suspend && clone.info_.progress() >= 1)
        tag = mcmp_tag::clone_halt;
    }
    broadcast(clone.work_, tag, 0);
    switch (tag) {
    case mcmp_tag::do_step: step<true>(clone, progress, advance); break;
    case mcmp_tag::clone_info: clone.send_info(mcmp_tag::clone_info); break;
    case mcmp_tag::clone_checkpoint:
      checkpoint<true>(clone);
      clone.send_info(mcmp_tag::clone_checkpoint);
      break;
    case mcmp_tag::clone_suspend:
      suspend<true>(clone);
      clone.send_info(mcmp_tag::clone_suspend);
      break;
    case mcmp_tag::clone_halt:
      halt<true>(clone);
      clone.send_halted();
      break;
    case mcmp_tag::do_nothing: break;
    default:
      if (leader<true>(clone))
        std::cerr << "Warning: ignoring a message with an unknown tag " << tag << std::endl;
    }
  }

  template<class Clone> static void request(Clone& clone, mcmp_tag_t tag) {
    if (leader<true>(clone)) {
      int message = tag;
      broadcast(clone.work_, message, 0);
      if (tag == mcmp_tag::clone_checkpoint) checkpoint<true>(clone);
      else suspend<true>(clone);
    }
  }
#endif
};

}} // namespace alps::detail
