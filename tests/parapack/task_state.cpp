// SPDX-License-Identifier: MIT
#include <alps/parapack/detail/task_state.hpp>
#include <alps/ngs/params.hpp>
#include <stdexcept>

void require(bool value) { if (!value) throw std::runtime_error("task state contract failed"); }

template<class Params>
struct test_task : alps::parapack::detail::task_state<test_task<Params>, Params> {
  using base = alps::parapack::detail::task_state<test_task<Params>, Params>;
  int saved = 0, invalidated = 0;
  void load() {
    this->task_id_ = 0;
    this->num_clones_ = typename base::range_type(1, 2);
    this->restore_state();
  }
  void save(bool) { ++saved; }
  void halt() { this->halt_state(); }
  void clear_observable_cache() { ++invalidated; }
};

struct proxy {
  alps::clone_info data{0};
  bool local = true, fresh = false;
  int destroyed = 0;
  proxy() { data.start("measurement"); }
  template<class Params>
  void start(alps::tid_t, alps::cid_t cid, alps::thread_group const&, Params const&,
             std::string const&, bool is_new) {
    fresh = is_new;
    if (fresh) { data = alps::clone_info(cid); data.start("measurement"); }
  }
  bool is_local(alps::Process) { return local; }
  alps::clone_info const& info(alps::Process) { return data; }
  void checkpoint(alps::Process) {}
  void update_info(alps::Process) {}
  void suspend(alps::Process) {}
  void halt(alps::Process) {}
  void destroy(alps::Process) { ++destroyed; }
};

template<class Params> void check() {
  using namespace alps;
  test_task<Params> task;
  proxy worker;
  thread_group group;
  auto cid = task.dispatch_clone(worker, group);
  require(cid && *cid == 0 && worker.fresh && task.num_running()==1);
  bool threw = false;
  try { task.halt(); } catch (std::logic_error const&) { threw = true; }
  require(threw);
  worker.data.set_progress(.4);
  task.checkpoint(worker, *cid);
  task.suspend_clone(worker, true, *cid, group);
  require(task.status()==task_status::Suspended && !task.on_memory() && task.saved==1);
  // Restore metadata as a resumed task would read it from its checkpoint.
  task.load();
  task.halt(); // Ready -> NotStarted
  require(task.status()==task_status::NotStarted);

  test_task<Params> resumed;
  resumed.load();
  auto first = resumed.dispatch_clone(worker, group);
  worker.data.set_progress(.4);
  worker.local = false;
  resumed.suspend_remote_clones(worker, true);
  resumed.clone_suspended(*first, group, worker.data);
  require(resumed.num_suspended()==1 && resumed.num_running()==0);
  worker.local = true;
  auto again = resumed.dispatch_clone(worker, group);
  require(again==first && !worker.fresh && resumed.num_suspended()==0);
  worker.data.set_progress(1);
  resumed.info_updated(*again,worker.data);
  resumed.halt_clone(worker,true,*again,group);
  require(resumed.status()==task_status::Finished && resumed.invalidated==1);

  test_task<Params> complete;
  for (int n=0;n<2;++n) {
    // Keep both clones running so finishing the first does not unload metadata.
    require(bool(complete.dispatch_clone(worker, group)));
  }
  for (int n=0;n<2;++n) {
    clone_info info(n); info.start("measurement"); info.set_progress(1);
    complete.info_updated(n,info);
    worker.data=info;
    complete.halt_clone(worker,true,n,group);
  }
  require(complete.status()==task_status::Completed && complete.invalidated==2);
}

int main() { check<alps::Parameters>(); check<alps::params>(); }
