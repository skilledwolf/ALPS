#include <chrono>
#include <thread>
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

#include <alps/ngs/parapack/parapack.h>
#include <alps/ngs/parapack/clone.h>
#include <alps/ngs/parapack/clone_proxy.h>
#include <alps/ngs/parapack/job_p.h>
#include <alps/parapack/filelock.h>
#include <alps/parapack/logger.h>
#include <alps/parapack/queue.h>
#include <alps/parapack/staging.h>
#include <alps/parapack/version.h>

#include <alps/config.h>
#include <alps/stop_callback.hpp>
#include <alps/utility/copyright.hpp>
#include <alps/utility/os.hpp>
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

namespace alps {

namespace ngs_parapack {

int start_impl(int argc, char **argv) {
  #ifndef BOOST_NO_EXCEPTIONS
  try {
  #endif

    alps::parapack::option opt(argc, argv);
    if (!opt.valid) {
      std::cerr << "Error: unknown command line option(s)\n";
      opt.print(std::cerr);
      return -1;
    }
    int ret;
    if (opt.jobfiles.size() == 0) {
      if (!opt.use_mpi) {
        if (opt.show_help) {
          opt.print(std::cout);
          return 0;
        }
        if (opt.show_license) {
          print_copyright(std::cout);
          print_license(std::cout);
          return 0;
        }
        ret = run_sequential(argc, argv);
      } else {
#ifdef ALPS_HAVE_MPI
        if (opt.show_help || opt.show_license) {
          boost::mpi::environment env(argc, argv);
          boost::mpi::communicator world;
          if (world.rank() == 0) {
            if (opt.show_help) {
              opt.print(std::cout);
            } else {
              print_copyright(std::cout);
              print_license(std::cout);
            }
          }
          return 0;
        }
        ret = run_sequential_mpi(argc, argv);
#else
        std::cerr << "ERROR: MPI is not supported\n";
        return -1;
#endif
      }
    } else {
      if (!opt.use_mpi) {
        ret = start_sgl(argc, argv);
      } else {
#ifdef ALPS_HAVE_MPI
        ret = start_mpi(argc, argv);
#else
        std::cerr << "ERROR: MPI is not supported\n";
        return -1;
#endif
      }
    }
    return ret;

  #ifndef BOOST_NO_EXCEPTIONS
  }
  catch (const std::exception& excp) {
    std::cerr << excp.what() << std::endl;
  }
  catch (...) {
    std::cerr << "Unknown exception occurred!" << std::endl;
  }
  return -1;
  #endif
}

// int evaluate(int argc, char **argv) {
//   #ifndef BOOST_NO_EXCEPTIONS
//   try {
//   #endif

//     evaluate_option opt(argc, argv);
//     if (!opt.valid) {
//       std::cerr << "Error: unknown command line option(s)\n";
//       opt.print(std::cerr);
//       return -1;
//     }
//     if (opt.show_help) {
//       opt.print(std::cout);
//       return 0;
//     }
//     if (opt.show_license) {
//       print_copyright(std::cout);
//       print_license(std::cout);
//       return 0;
//     }

//     BOOST_FOREACH(std::string const& file_str, opt.jobfiles) {
//       boost::filesystem::path file = absolute(boost::filesystem::path(file_str)).normalize();
//       if (!exists(file)) {
//         std::cerr << "Error: file not found: " << file << std::endl;
//         return -1;
//       }
//       boost::filesystem::path basedir = file.parent_path();
//       std::string file_in_str;
//       std::string file_out_str;
//       std::vector<task> tasks;

//       std::cout << logger::header() << "starting evaluation on "
//                 << alps::hostname() << std::endl;
//       int t = load_filename(file, file_in_str, file_out_str);
//       if (t == 1) {
//         // process all tasks
//         boost::filesystem::path file_in = absolute(boost::filesystem::path(file_in_str), basedir);
//         boost::filesystem::path file_out = absolute(boost::filesystem::path(file_out_str), basedir);
//         std::string simname;
//         load_tasks(file_in, file_out, basedir, simname, tasks, false, opt.write_xml);
//         std::cout << "  master input file  = " << file_in.string() << std::endl
//                   << "  master output file = " << file_out.string() << std::endl;
//         print_taskinfo(std::cout, tasks, opt.task_range);
//         BOOST_FOREACH(task& t, tasks) {
//           if (!opt.task_range.valid() || opt.task_range.is_included(t.task_id()+1))
//             t.evaluate(opt.write_xml);
//         }
//       } else {
//         // process one task
//         task t(file);
//         if (!opt.task_range.valid() || opt.task_range.is_included(t.task_id()+1))
//           t.evaluate(opt.write_xml);
//       }
//       std::cout << logger::header() << "all tasks evaluated\n";
//     }
//     return 0;

//   #ifndef BOOST_NO_EXCEPTIONS
//   }
//   catch (const std::exception& excp) {
//     std::cerr << excp.what() << std::endl;
//   }
//   catch (...) {
//     std::cerr << "known exception occurred!" << std::endl;
//   }
//   return -1;
//   #endif
// }

void print_copyright(std::ostream& os) {
  // worker_factory::print_copyright(os);
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

void print_taskinfo(std::ostream& os, std::vector<alps::ngs_parapack::task> const& tasks,
  task_range_t const& task_range) {
  uint32_t num_new = 0;
  uint32_t num_running = 0;
  uint32_t num_continuing = 0;
  uint32_t num_suspended = 0;
  uint32_t num_finished = 0;
  uint32_t num_completed = 0;
  uint32_t num_skipped = 0;
  BOOST_FOREACH(alps::ngs_parapack::task const& t, tasks) {
    if (!task_range.valid() || task_range.is_included(t.task_id()+1)) {
      switch (t.status()) {
      case alps::task_status::NotStarted :
        ++num_new;
        break;
      case alps::task_status::Running :
        ++num_running;
        break;
      case alps::task_status::Continuing :
        ++num_continuing;
        break;
      case alps::task_status::Suspended :
        ++num_suspended;
        break;
      case alps::task_status::Finished :
        ++num_finished;
        break;
      case alps::task_status::Completed :
        ++num_completed;
        break;
      default :
        break;
      }
    } else {
      ++num_skipped;
    }
  }
  os << logger::header() << "task status: "
     << "total number of tasks = " << tasks.size() << std::endl
     << "  new = " << num_new
     << ", running = " << num_running
     << ", continuing = " << num_continuing
     << ", suspended = " << num_suspended
     << ", finished = " << num_finished
     << ", completed = " << num_completed
     << ", skipped = " << num_skipped << std::endl;
}

int load_filename(boost::filesystem::path const& file, std::string& file_in_str,
  std::string& file_out_str) {
  bool is_master;
  alps::ngs_parapack::filename_xml_handler handler(file_in_str, file_out_str, is_master);
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
  alps::ngs_parapack::version_xml_handler handler(versions);
  alps::XMLParser parser(handler);
  parser.parse(file);
}

void load_tasks(boost::filesystem::path const& file_in, boost::filesystem::path const& file_out,
  boost::filesystem::path const& basedir, std::string& simname, std::vector<alps::ngs_parapack::task>& tasks,
  bool check_parameter, bool write_xml) {
  tasks.clear();
  if (!exists(file_out)) {
    alps::ngs_parapack::job_tasks_xml_handler handler(simname, tasks, basedir);
    alps::XMLParser parser(handler);
    parser.parse(file_in);
  } else {
    alps::ngs_parapack::job_tasks_xml_handler handler_out(simname, tasks, basedir);
    alps::XMLParser parser_out(handler_out);
    parser_out.parse(file_out);

    if (check_parameter) {
      std::vector<alps::ngs_parapack::task> tasks_in;
      alps::ngs_parapack::job_tasks_xml_handler handler_in(simname, tasks_in, basedir);
      alps::XMLParser parser_in(handler_in);
      parser_in.parse(file_in);

      int nc = std::min BOOST_PREVENT_MACRO_SUBSTITUTION (tasks_in.size(), tasks.size());
      for (int i = 0; i < nc; ++i) {
        if (tasks_in[i].file_in_str() != tasks[i].file_in_str() ||
            tasks_in[i].file_out_str() != tasks[i].file_out_str()) {
          std::cout << "Info: input/output XML filename of " << logger::task(i)
                    << " has been modified" << std::endl;
          tasks[i] = tasks_in[i];
        }
        tasks[i].check_parameter(write_xml);
      }
      if (tasks_in.size() > tasks.size()) {
        std::cout << "Info: number of parameter sets has been increased from " << tasks.size()
                  << " to " << tasks_in.size() << std::endl;
        for (int i = tasks.size(); i < tasks_in.size(); ++i) tasks.push_back(tasks_in[i]);
      } else if (tasks_in.size() < tasks.size()) {
        std::cout << "Info: number of parameter sets has been decreased from " << tasks.size()
                  << " to " << tasks_in.size() << std::endl;
        tasks.resize(tasks_in.size());
      }
    }
  }
}


void save_tasks(boost::filesystem::path const& file, std::string const& simname,
  std::string const& file_in_str, std::string const& file_out_str,
  std::vector<alps::ngs_parapack::task>& tasks) {
  alps::ngs_parapack::job_xml_writer(file, simname, file_in_str, file_out_str, alps_version(),
    "unknown", tasks, true);
}

int run_sequential(int argc, char **argv) {

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  alps::ParameterList parameterlist;
  std::cin >> parameterlist;

#ifdef _OPENMP
  // set default number of threads to 1
  char* nt = getenv("OMP_NUM_THREADS");
  if (nt == 0 && omp_get_max_threads() != 1) omp_set_num_threads(1);
#endif

  for (std::size_t i = 0; i < parameterlist.size(); ++i) {
    alps::params p;
    for (alps::Parameters::const_iterator itr = parameterlist[i].begin(); itr != parameterlist[i].end(); ++itr) {
      p[itr->key()] = itr->value();
    }
    const std::clock_t tm = std::clock();
    if (!p.defined("DIR_NAME")) p["DIR_NAME"] = ".";
    if (!p.defined("BASE_NAME")) p["BASE_NAME"] = "task" + boost::lexical_cast<std::string>(i+1);
    if (!p.defined("SEED")) p["SEED"] = static_cast<unsigned int>(std::time(0));
    p["WORKER_SEED"] = p["SEED"];
    p["DISORDER_SEED"] = p["SEED"];
    std::cout << "[input parameters]\n" << p << std::flush;
    boost::shared_ptr<alps::ngs_parapack::abstract_worker>
      worker = alps::ngs_parapack::worker_factory::make_worker(p);
    while (worker->fraction_completed() < 1.0) {
      worker->run(stop_callback(0), boost::function<void (double)>());
    }
    std::cerr << "[speed]\nelapsed time = "
              << static_cast<double>(std::clock() - tm) / CLOCKS_PER_SEC << " sec\n";
    // std::cout << "[results]\n" << collect_results(*worker);
    std::cout << std::flush;
  }

#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (const std::exception& excp) {
    std::cerr << excp.what() << std::endl;
    alps::comm_exit(true);
    return -1; }
  catch (...) {
    std::cerr << "Unknown exception occurred!" << std::endl;
    alps::comm_exit(true);
    return -1; }
#endif
  return 0;
}

namespace {
struct scheduler_runtime {
  using task_type = alps::ngs_parapack::task;
  using clone_type = alps::ngs_parapack::clone;
  static constexpr int error_code = -1;
  static constexpr bool supports_evaluation = false;
  static constexpr bool reports_vmusage = false;
  static constexpr bool validate_thread_layout = false;
  static constexpr auto print_copyright = &alps::ngs_parapack::print_copyright;
  static constexpr auto print_license = &alps::ngs_parapack::print_license;
  static constexpr auto load_filename = &alps::ngs_parapack::load_filename;
  static constexpr auto save_tasks = &alps::ngs_parapack::save_tasks;
  static constexpr auto print_taskinfo = &alps::ngs_parapack::print_taskinfo;
  template<class... Args> static void load_tasks(Args&&... args) {
    alps::ngs_parapack::load_tasks(std::forward<Args>(args)...);
  }
  static bool task_options(alps::parapack::option const& opt) {
    return opt.write_xml;
  }
  template<class Clone> static void run(Clone& worker) {
    worker.run(stop_callback(0), boost::function<void (double)>());
  }
  struct clone_proxy_type : alps::ngs_parapack::clone_proxy {
    clone_proxy_type(clone_type*& worker, boost::filesystem::path const& basedir, alps::parapack::option const& opt)
      : alps::ngs_parapack::clone_proxy(worker, basedir, opt.dump_policy, opt.check_interval) {}
  };
#ifdef ALPS_HAVE_MPI
  using clone_mpi_type = alps::ngs_parapack::clone_mpi;
  using clone_info_msg_t_type = alps::ngs_parapack::clone_info_msg_t;
  using clone_halt_msg_t_type = alps::ngs_parapack::clone_halt_msg_t;
  using clone_create_msg_t_type = alps::ngs_parapack::clone_create_msg_t;
  static int receive_source(boost::mpi::status const& status) {
    return boost::mpi::any_source;
  }
  struct clone_proxy_mpi_type : alps::ngs_parapack::clone_proxy_mpi {
    clone_proxy_mpi_type(clone_mpi_type*& worker, boost::mpi::communicator const& ctrl,
                        boost::mpi::communicator const& work, boost::filesystem::path const& basedir,
                        alps::parapack::option const& opt)
      : alps::ngs_parapack::clone_proxy_mpi(worker, ctrl, work, basedir, opt.dump_policy, opt.check_interval) {}
  };
  static clone_mpi_type* make_clone(boost::mpi::communicator const& ctrl,
                                  boost::mpi::communicator const& work, boost::filesystem::path const& basedir,
                                  alps::parapack::option const& opt, clone_create_msg_t_type const& msg) {
    return new clone_mpi_type(ctrl, work, basedir, opt.dump_policy, opt.check_interval, msg);
  }
#endif
};
} // namespace

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
    alps::params p;
    for (alps::Parameters::const_iterator itr = parameterlist[i].begin(); itr != parameterlist[i].end(); ++itr) {
      p[itr->key()] = itr->value();
    }
    world.barrier();
    const std::clock_t tm = std::clock();
    if (!p.defined("DIR_NAME")) p["DIR_NAME"] = ".";
    if (!p.defined("BASE_NAME")) p["BASE_NAME"] = "task" + boost::lexical_cast<std::string>(i+1);
    if (!p.defined("SEED")) p["SEED"] = static_cast<unsigned int>(std::time(0));
    p["WORKER_SEED"] = static_cast<unsigned int>(p["SEED"]) ^ (world.rank() << 11);
    p["DISORDER_SEED"] = p["SEED"];
    if (world.rank() == 0) std::cout << "[input parameters]\n" << p << std::flush;
    boost::shared_ptr<alps::ngs_parapack::abstract_worker>
      worker = parallel_worker_factory::make_worker(world, p);
    while (worker->fraction_completed() < 1.0) {
      worker->run(stop_callback(0), boost::function<void (double)>());
    }
    world.barrier();
    if (world.rank() == 0) {
      std::cerr << "[speed]\nelapsed time = "
                << static_cast<double>(std::clock() - tm) / CLOCKS_PER_SEC << " sec" << std::endl;
    }
    for (int r = 0; r < world.size(); ++r) {
      if (world.rank() == r) {
        // std::cout << "[results " << r << "]\n" << collect_results(*worker);
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
    return -1; }
  catch (...) {
    std::cerr << "Unknown exception occurred!" << std::endl;
    alps::comm_exit(true);
    return -1; }
#endif
  return 0;
}

#else // ALPS_HAVE_MPI

int run_sequential_mpi(int argc, char** argv) {
  std::cerr << "This program has not been compiled for use with MPI\n";
  return -1;
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
  return -1;
}

#endif // ALPS_HAVE_MPI

} // end namespace parapack
} // end namespace alps
