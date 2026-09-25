// Copyright (C) 1997-2012 by Synge Todo and the ALPS contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include <alps/parapack/detail/job_xml.hpp>
#include <alps/parapack/logger.h>
#include <alps/parapack/option.h>
#include <algorithm>

namespace alps { namespace parapack { namespace detail {

template<class Task>
void print_taskinfo(std::ostream& os, std::vector<Task> const& tasks,
  task_range_t const& task_range) {
  uint32_t num_new = 0;
  uint32_t num_running = 0;
  uint32_t num_continuing = 0;
  uint32_t num_suspended = 0;
  uint32_t num_finished = 0;
  uint32_t num_completed = 0;
  uint32_t num_skipped = 0;
  for (Task const& t : tasks) {
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

// The XML task representation is shared; each runtime supplies the task type
// and its parameter-check options.
template<class Task>
void read_tasks(boost::filesystem::path const& file, boost::filesystem::path const& basedir,
  std::string& simname, std::vector<Task>& tasks) {
  tasks.clear();
  job_tasks_xml_handler<Task> handler(simname, tasks, basedir);
  XMLParser parser(handler);
  parser.parse(file);
}

template<class Task>
void load_tasks(boost::filesystem::path const& file_in, boost::filesystem::path const& file_out,
  boost::filesystem::path const& basedir, std::string& simname, std::vector<Task>& tasks) {
  read_tasks(exists(file_out) ? file_out : file_in, basedir, simname, tasks);
}

template<class Task, class Options>
void load_tasks(boost::filesystem::path const& file_in, boost::filesystem::path const& file_out,
  boost::filesystem::path const& basedir, std::string& simname, std::vector<Task>& tasks,
  bool check_parameter, Options const& options) {
  bool restarting = exists(file_out);
  read_tasks(restarting ? file_out : file_in, basedir, simname, tasks);
  if (!restarting || !check_parameter) return;

  std::vector<Task> tasks_in;
  read_tasks(file_in, basedir, simname, tasks_in);
  std::size_t common = (std::min)(tasks_in.size(), tasks.size());
  for (std::size_t i = 0; i < common; ++i) {
    if (tasks_in[i].file_in_str() != tasks[i].file_in_str() ||
        tasks_in[i].file_out_str() != tasks[i].file_out_str()) {
      std::cout << "Info: input/output XML filename of " << logger::task(i)
                << " has been modified" << std::endl;
      tasks[i] = tasks_in[i];
    }
    tasks[i].check_parameter(options);
  }
  if (tasks_in.size() > tasks.size()) {
    std::cout << "Info: number of parameter sets has been increased from " << tasks.size()
              << " to " << tasks_in.size() << std::endl;
    tasks.insert(tasks.end(), tasks_in.begin() + common, tasks_in.end());
  } else if (tasks_in.size() < tasks.size()) {
    std::cout << "Info: number of parameter sets has been decreased from " << tasks.size()
              << " to " << tasks_in.size() << std::endl;
    tasks.resize(tasks_in.size());
  }
}

}}} // namespace alps::parapack::detail
