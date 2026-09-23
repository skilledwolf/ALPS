// SPDX-License-Identifier: MIT
#include <alps/parapack/parapack.h>
#include <alps/ngs/parapack/parapack.h>
#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>
#include <stdexcept>

void require(bool condition) {
  if (!condition) throw std::runtime_error("job reconciliation contract failed");
}

void write_job(boost::filesystem::path const& path,
  std::initializer_list<std::pair<std::string, std::string>> entries) {
  std::ofstream out(path.string());
  out << "<JOB name=\"test\">";
  for (auto const& [name, status] : entries)
    out << "<TASK status=\"" << status << "\"><INPUT file=\"" << name
        << ".in.xml\"/><OUTPUT file=\"" << name << ".out.xml\"/></TASK>";
  out << "</JOB>";
}

template<class Task, class Load, class Report>
void check(Load load, Report report) {
  auto dir = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path();
  boost::filesystem::create_directory(dir);
  struct cleanup {
    boost::filesystem::path path;
    ~cleanup() { boost::filesystem::remove_all(path); }
  } cleanup{dir};
  auto input = dir / "job.in.xml", output = dir / "job.out.xml";
  std::string name;
  std::vector<Task> tasks;
  auto read = [&](bool reconcile) { load(input, output, dir, name, tasks, reconcile); };

  write_job(input, {{"a", "new"}, {"b", "new"}});
  read(true);
  require(name == "test" && tasks.size() == 2 && tasks[1].task_id() == 1);
  write_job(output, {{"a", "completed"}, {"b", "running"}});
  write_job(input, {{"a", "new"}, {"renamed", "new"}, {"c", "new"}});
  read(false);
  require(tasks.size() == 2 && tasks[0].status() == alps::task_status::Completed);
  require(tasks[1].file_in_str() == "b.in.xml");
  read(true);
  require(tasks.size() == 3 && tasks[0].status() == alps::task_status::Completed);
  require(tasks[1].file_in_str() == "renamed.in.xml" && tasks[1].status() == alps::task_status::NotStarted);
  require(tasks[2].task_id() == 2 && tasks[2].file_in_str() == "c.in.xml");
  write_job(input, {{"replacement", "new"}});
  read(true);
  require(tasks.size() == 1 && tasks[0].file_out_str() == "replacement.out.xml");
  write_job(input, {});
  read(true);
  require(tasks.empty());

  write_job(output, {{"a", "new"}, {"b", "running"}, {"c", "finished"}, {"d", "completed"}});
  read(false);
  std::ostringstream all, selected;
  report(all, tasks, alps::task_range_t());
  report(selected, tasks, alps::task_range_t(2, 3));
  require(all.str().find("new = 1, running = 0, continuing = 0, suspended = 1, finished = 1, completed = 1, skipped = 0") != std::string::npos);
  require(selected.str().find("new = 0, running = 0, continuing = 0, suspended = 1, finished = 1, completed = 0, skipped = 2") != std::string::npos);
}

int main() {
  char program[] = "job_workflow";
  char* argv[] = {program};
  alps::parapack::option options(1, argv);
  check<alps::task>([&](auto const& in, auto const& out, auto const& dir, auto& name, auto& tasks, bool reconcile) {
    alps::parapack::load_tasks(in, out, dir, name, tasks, reconcile, options);
  }, alps::parapack::print_taskinfo);
  check<alps::ngs_parapack::task>([](auto const& in, auto const& out, auto const& dir, auto& name, auto& tasks, bool reconcile) {
    alps::ngs_parapack::load_tasks(in, out, dir, name, tasks, reconcile, true);
  }, alps::ngs_parapack::print_taskinfo);
}
