#include <chrono>
#include <thread>
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

#include "parapack.h"
#include "clone.h"
#include "clone_proxy.h"
#include "filelock.h"
#include "job_p.h"
#include "logger.h"
#include "queue.h"
#include "staging.h"
#include "version.h"

#include <alps/config.h>
#include <alps/utility/copyright.hpp>
#include <alps/utility/os.hpp>
#include <alps/utility/vmusage.hpp>
#include <alps/osiris/comm.h>

// some file (probably a python header) defines a tolower macro ...
#undef tolower
#undef toupper

#include <boost/config.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>

#include <ctime>
#include <iostream>
#include <string>
#include <vector>


#ifdef _OPENMP
# include <omp.h>
#endif

#include <alps/parapack/scheduler_workflow.hpp>
#include <alps/parapack/detail/job_workflow.hpp>

namespace alps {

namespace parapack {

int evaluate(int argc, char **argv) {
  #ifndef BOOST_NO_EXCEPTIONS
  try {
  #endif

    option opt(argc, argv, /* for_evaluate = */ true);
    if (!opt.valid) {
      std::cerr << "Error: unknown command line option(s)\n";
      opt.print(std::cerr);
      return 127;
    }
    if (opt.show_help) {
      opt.print(std::cout);
      return 0;
    }
    if (opt.show_license) {
      print_copyright(std::cout);
      print_license(std::cout);
      return 0;
    }

    BOOST_FOREACH(std::string const& file_str, opt.jobfiles) {
      boost::filesystem::path file = absolute(boost::filesystem::path(file_str)).lexically_normal();
      if (!exists(file)) {
        std::cerr << "Error: file not found: " << file << std::endl;
        return 127;
      }
      boost::filesystem::path basedir = file.parent_path();
      std::string file_in_str;
      std::string file_out_str;
      std::vector<task> tasks;

      std::cout << logger::header() << "starting evaluation on "
                << alps::hostname() << std::endl;
      int t = load_filename(file, file_in_str, file_out_str);
      if (t == 1) {
        // process all tasks
        boost::filesystem::path file_in = absolute(boost::filesystem::path(file_in_str), basedir);
        boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str), basedir);
        std::string simname;
        load_tasks(file_in, file_out, basedir, simname, tasks, false, opt);
        std::cout << "  master input file  = " << file_in.string() << std::endl
                  << "  master output file = " << file_out.string() << std::endl;
        print_taskinfo(std::cout, tasks, opt.task_range);
        BOOST_FOREACH(task& t, tasks) {
          if (!opt.task_range.valid() || opt.task_range.is_included(t.task_id()+1))
            t.evaluate(opt);
        }
      } else {
        // process one task
        task t(file);
        if (!opt.task_range.valid() || opt.task_range.is_included(t.task_id()+1))
          t.evaluate(opt);
      }
      std::cout << logger::header() << "all tasks evaluated\n";
    }
    return 0;

  #ifndef BOOST_NO_EXCEPTIONS
  }
  catch (const std::exception& excp) {
    std::cerr << excp.what() << std::endl;
  }
  catch (...) {
    std::cerr << "known exception occurred!" << std::endl;
  }
  return 127;
  #endif
}

void print_copyright(std::ostream& os) {
  worker_factory::print_copyright(os);
  os << std::endl << "using " << parapack_copyright() << std::endl;
  alps::print_copyright(os);
}

void print_license(std::ostream& os) {
  alps::print_license(os);
}

std::string alps_version() {
  return version_string() + "; configured on " + config_host() +
    " by " + config_user() + "; compiled on " + compile_date();
}

void print_taskinfo(std::ostream& os, std::vector<alps::task> const& tasks,
  task_range_t const& task_range) {
  alps::parapack::detail::print_taskinfo(os, tasks, task_range);
}

int load_filename(boost::filesystem::path const& file, std::string& file_in_str,
  std::string& file_out_str) {
  bool is_master;
  alps::filename_xml_handler handler(file_in_str, file_out_str, is_master);
  alps::XMLParser parser(handler);
  parser.parse(file);
  if (is_master) {
    if (file_out_str.empty())
      file_out_str = file.filename().string();
    if (file_in_str.empty())
      file_in_str = regex_replace(file_out_str, boost::regex("\\.out\\.xml$"), ".in.xml");
  }
  return is_master ? 1 : 2;
}

void load_version(boost::filesystem::path const& file,
  std::vector<std::pair<std::string, std::string> >& versions) {
  alps::version_xml_handler handler(versions);
  alps::XMLParser parser(handler);
  parser.parse(file);
}

void load_tasks(boost::filesystem::path const& file_in, boost::filesystem::path const& file_out,
  boost::filesystem::path const& basedir, std::string& simname, std::vector<alps::task>& tasks,
  bool check_parameter, alps::parapack::option const& opt) {
  alps::parapack::detail::load_tasks(file_in, file_out, basedir, simname, tasks,
    check_parameter, opt);
}

void load_tasks(boost::filesystem::path const& file_in, boost::filesystem::path const& file_out,
  boost::filesystem::path const& basedir, std::string& simname, std::vector<alps::task>& tasks) {
  alps::parapack::detail::load_tasks(file_in, file_out, basedir, simname, tasks);
}

void save_tasks(boost::filesystem::path const& file, std::string const& simname,
  std::string const& file_in_str, std::string const& file_out_str,
  std::vector<alps::task>& tasks) {
  alps::job_xml_writer(file, simname, file_in_str, file_out_str, alps_version(),
    worker_factory::version(), tasks, true);
}

int run_sequential(int argc, char **argv) {

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  alps::ParameterList parameterlist;
  std::cin >> parameterlist;

#ifdef _OPENMP
  // set default number of threads to 1
  char* p = getenv("OMP_NUM_THREADS");
  if (p == 0 && omp_get_max_threads() != 1) omp_set_num_threads(1);
#endif

  for (std::size_t i = 0; i < parameterlist.size(); ++i) {
    alps::Parameters p = parameterlist[i];
    const std::clock_t tm = std::clock();
    if (!p.defined("DIR_NAME")) p["DIR_NAME"] = ".";
    if (!p.defined("BASE_NAME")) p["BASE_NAME"] = "task" + boost::lexical_cast<std::string>(i+1);
    if (!p.defined("CLONE_ID")) p["CLONE_ID"] = 1;
    if (!p.defined("SEED")) p["SEED"] = static_cast<unsigned int>(std::time(0));
    p["WORKER_SEED"] = p["SEED"];
    p["DISORDER_SEED"] = p["SEED"];
    std::cout << "[input parameters]\n" << p << std::flush;
    std::vector<alps::ObservableSet> obs;
    boost::shared_ptr<alps::parapack::abstract_worker>
      worker = worker_factory::make_worker(p);
    worker->init_observables(p, obs);
    bool thermalized = worker->is_thermalized();
    if (thermalized) for (std::size_t r = 0; r < obs.size(); ++r) obs[r].reset(true);
    while (worker->progress() < 1.0) {
      worker->run(obs);
      if (!thermalized && worker->is_thermalized()) {
        for (std::size_t r = 0; r < obs.size(); ++r) obs[r].reset(true);
        thermalized = true;
      }
    }
    std::vector<alps::ObservableSet> obs_out;
    boost::shared_ptr<alps::parapack::abstract_evaluator>
      evaluator = evaluator_factory::make_evaluator(p);
    evaluator->load(obs, obs_out);
    evaluator->evaluate(obs_out);
    std::cerr << "[speed]\nelapsed time = "
              << static_cast<double>(std::clock() - tm) / CLOCKS_PER_SEC << " sec\n";
    std::cout << "[results]\n";
    if (obs_out.size() == 1) {
      std::cout << obs_out[0];
    } else {
      for (std::size_t r = 0; r < obs_out.size(); ++r)
        std::cout << "[[replica " << r << "]]\n" << obs_out[r];
    }
    std::cout << std::flush;
  }

#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (const std::exception& excp) {
    std::cerr << excp.what() << std::endl;
    alps::comm_exit(true);
    return 127; }
  catch (...) {
    std::cerr << "Unknown exception occurred!" << std::endl;
    alps::comm_exit(true);
    return 127; }
#endif
  return 0;
}

namespace {
struct scheduler_runtime {
  static constexpr auto start_sgl = &alps::parapack::start_sgl;
  static constexpr auto run_sequential = &alps::parapack::run_sequential;
#ifdef ALPS_HAVE_MPI
  static constexpr auto start_mpi = &alps::parapack::start_mpi;
  static constexpr auto run_sequential_mpi = &alps::parapack::run_sequential_mpi;
#endif
  using task_type = alps::task;
  using clone_type = alps::clone;
  static constexpr int error_code = 127;
  static constexpr bool supports_evaluation = true;
  static constexpr bool reports_vmusage = true;
  static constexpr bool validate_thread_layout = true;
  static constexpr auto print_copyright = &alps::parapack::print_copyright;
  static constexpr auto print_license = &alps::parapack::print_license;
  static constexpr auto load_filename = &alps::parapack::load_filename;
  static constexpr auto save_tasks = &alps::parapack::save_tasks;
  static constexpr auto print_taskinfo = &alps::parapack::print_taskinfo;
  template<class... Args> static void load_tasks(Args&&... args) {
    alps::parapack::load_tasks(std::forward<Args>(args)...);
  }
  static alps::parapack::option const& task_options(alps::parapack::option const& opt) {
    return opt;
  }
  template<class Clone> static void run(Clone& worker) {
    worker.run();
  }
  struct clone_proxy_type : alps::clone_proxy {
    clone_proxy_type(clone_type*& worker, boost::filesystem::path const& basedir, alps::parapack::option const& opt)
      : alps::clone_proxy(worker, basedir, opt) {}
  };
#ifdef ALPS_HAVE_MPI
  using clone_mpi_type = alps::clone_mpi;
  using clone_info_msg_t_type = alps::clone_info_msg_t;
  using clone_halt_msg_t_type = alps::clone_halt_msg_t;
  using clone_create_msg_t_type = alps::clone_create_msg_t;
  static int receive_source(boost::mpi::status const& status) {
    return status.source();
  }
  struct clone_proxy_mpi_type : alps::clone_proxy_mpi {
    clone_proxy_mpi_type(clone_mpi_type*& worker, boost::mpi::communicator const& ctrl,
                        boost::mpi::communicator const& work, boost::filesystem::path const& basedir,
                        alps::parapack::option const& opt)
      : alps::clone_proxy_mpi(worker, ctrl, work, basedir, opt) {}
  };
  static clone_mpi_type* make_clone(boost::mpi::communicator const& ctrl,
                                  boost::mpi::communicator const& work, boost::filesystem::path const& basedir,
                                  alps::parapack::option const& opt, clone_create_msg_t_type const& msg) {
    return new clone_mpi_type(ctrl, work, basedir, opt, msg);
  }
#endif
};
} // namespace

int start(int argc, char** argv) {
  return alps::parapack::detail::start<scheduler_runtime>(argc, argv);
}

int start_sgl(int argc, char** argv) {
  return alps::parapack::detail::start_sgl<scheduler_runtime>(argc, argv);
}

#ifdef ALPS_HAVE_MPI

int run_sequential_mpi(int argc, char** argv) {

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

#ifdef _OPENMP
  // set default number of threads to 1
  char* p = getenv("OMP_NUM_THREADS");
  if (p == 0 && omp_get_max_threads() != 1) omp_set_num_threads(1);
#endif

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  alps::ParameterList parameterlist;
  if (world.rank() == 0) std::cin >> parameterlist;
  broadcast(world, parameterlist, 0);

  for (int i = 0; i < parameterlist.size(); ++i) {
    alps::Parameters p = parameterlist[i];
    world.barrier();
    const std::clock_t tm = std::clock();
    if (!p.defined("DIR_NAME")) p["DIR_NAME"] = ".";
    if (!p.defined("BASE_NAME")) p["BASE_NAME"] = "task" + boost::lexical_cast<std::string>(i+1);
    if (!p.defined("SEED")) p["SEED"] = static_cast<unsigned int>(std::time(0));
    p["WORKER_SEED"] = static_cast<unsigned int>(p["SEED"]) ^ (world.rank() << 11);
    p["DISORDER_SEED"] = p["SEED"];
    if (world.rank() == 0) std::cout << "[input parameters]\n" << p << std::flush;
    std::vector<alps::ObservableSet> obs;
    boost::shared_ptr<alps::parapack::abstract_worker>
      worker = parallel_worker_factory::make_worker(world, p);
    worker->init_observables(p, obs);
    bool thermalized = worker->is_thermalized();
    if (thermalized) {
      BOOST_FOREACH(alps::ObservableSet& o, obs) { o.reset(true); }
    }
    while (worker->progress() < 1.0) {
      worker->run(obs);
      if (!thermalized && worker->is_thermalized()) {
        BOOST_FOREACH(alps::ObservableSet& o, obs) { o.reset(true); }
        thermalized = true;
      }
    }
    world.barrier();
    if (world.rank() == 0) {
      std::cerr << "[speed]\nelapsed time = "
                << static_cast<double>(std::clock() - tm) / CLOCKS_PER_SEC << " sec" << std::endl;
    }
    std::vector<alps::ObservableSet> obs_out;
    boost::shared_ptr<alps::parapack::abstract_evaluator>
      evaluator = evaluator_factory::make_evaluator(p);
    evaluator->load(obs, obs_out);
    evaluator->evaluate(obs_out);
    for (int r = 0; r < world.size(); ++r) {
      if (world.rank() == r) {
        std::cout << "[results " << r << "]\n";
        if (obs_out.size() == 1) {
          std::cout << obs_out[0];
        } else {
          for (int i = 0; i < obs_out.size(); ++i)
            std::cout << "[[replica " << i << "]]\n" << obs_out[i];
        }
        std::cout << std::flush;
      }
      world.barrier();
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (const std::exception& excp) {
    std::cerr << excp.what() << std::endl;
    alps::comm_exit(true);
    return 127; }
  catch (...) {
    std::cerr << "Unknown exception occurred!" << std::endl;
    alps::comm_exit(true);
    return 127; }
#endif
  return 0;
}

#else // ALPS_HAVE_MPI

int run_sequential_mpi(int argc, char** argv) {
  std::cerr << "This program has not been compiled for use with MPI\n";
  return 127;
}

#endif // ALPS_HAVE_MPI

#ifdef ALPS_HAVE_MPI

int start_mpi(int argc, char** argv) {
  return alps::parapack::detail::start_mpi<scheduler_runtime>(argc, argv);
}

#else // ALPS_HAVE_MPI

int start_mpi(int, char**) {
  boost::throw_exception(std::runtime_error(
    "This program has not been compiled for use with MPI"));
  return 127;
}

#endif // ALPS_HAVE_MPI

} // end namespace parapack
} // end namespace alps
