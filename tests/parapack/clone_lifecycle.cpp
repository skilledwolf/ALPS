// SPDX-License-Identifier: MIT
#include <alps/parapack/clone.h>
#include <alps/ngs/parapack/clone.h>
#include <boost/filesystem/operations.hpp>
#include <memory>
#include <stdexcept>

namespace {
void require(bool ok, char const* message) {
  if (!ok) throw std::runtime_error(message);
}

class classic_worker : public alps::parapack::abstract_worker {
public:
  explicit classic_worker(alps::Parameters const& p)
    : count_(p.value_or_default("INITIAL_COUNT", 0)) {}
#ifdef ALPS_HAVE_MPI
  classic_worker(boost::mpi::communicator const&, alps::Parameters const& p) : classic_worker(p) {}
#endif
  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) override {
    obs << alps::RealObservable("Count");
  }
  void run(alps::ObservableSet& obs) override { obs["Count"] << double(++count_); }
  void load(alps::IDump& dump) override { dump >> count_; }
  void save(alps::ODump& dump) const override { dump << count_; }
  bool is_thermalized() const override { return count_ >= 2; }
  double progress() const override { return count_ / 6.; }
private:
  int count_;
};

class ngs_worker {
public:
  explicit ngs_worker(alps::params const& p) : count_(p["INITIAL_COUNT"] | 0) {}
#ifdef ALPS_HAVE_MPI
  ngs_worker(boost::mpi::communicator const&, alps::params const& p) : ngs_worker(p) {}
#endif
  bool run(boost::function<bool()> const& stop, boost::function<void(double)> const& progress) {
    if (stop()) return false;
    ++count_;
    progress(fraction_completed());
    return fraction_completed() >= 1;
  }
  void load(alps::hdf5::archive& ar) { ar["/test/count"] >> count_; }
  void save(alps::hdf5::archive& ar) const { ar["/test/count"] << count_; }
  double fraction_completed() const { return count_ / 6.; }
private:
  int count_;
};

template<class Factory, class Run, class Verify>
void exercise(Factory make, Run run, Verify verify, bool leader, int steps) {
  auto first = make(true);
  for (int i=0; i<steps; ++i) run(*first);
  if (leader) first->checkpoint(); else run(*first);
  require(!first->halted(), "checkpoint halted the worker");
  if (leader) first->suspend(); else run(*first);
  require(first->halted(), "suspend kept the worker alive");
  auto restored = make(false);
  if (leader) require(restored->info().progress() == steps/6., "resume lost checkpoint progress");
  int calls = 0;
  while (!restored->halted() && calls++ < 10) run(*restored);
  require(restored->halted(), "resumed worker did not finish");
  if (leader) require(restored->info().progress() == 1., "completion lost final progress");
  verify(restored->info());
}
}

int main(int argc, char** argv) {
  bool mpi = argc > 1;
  bool leader = true;
  int rank = 0;
#ifdef ALPS_HAVE_MPI
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  leader = !mpi || world.rank() == 0;
  rank = world.rank();
#else
  if (mpi) return 2;
#endif
  auto directory = boost::filesystem::temp_directory_path() /
    boost::filesystem::unique_path("alps-clone-%%%%-%%%%");
#ifdef ALPS_HAVE_MPI
  if (mpi) {
    std::string path = directory.string();
    broadcast(world, path, 0);
    directory = path;
  }
#endif
  if (leader) boost::filesystem::create_directory(directory);
#ifdef ALPS_HAVE_MPI
  if (mpi) world.barrier();
#endif
  alps::parapack::worker_factory::instance()->register_worker<classic_worker>("lifecycle");
  alps::ngs_parapack::worker_factory::register_worker<ngs_worker>();
#ifdef ALPS_HAVE_MPI
  alps::parapack::parallel_worker_factory::instance()->register_worker<classic_worker>("lifecycle");
  alps::ngs_parapack::parallel_worker_factory::register_worker<ngs_worker>();
#endif
  char name[] = "clone_lifecycle";
  char* args[] = {name};
  alps::parapack::option opt(1, args);
  opt.dump_policy = alps::dump_policy::All;
  opt.check_interval = boost::posix_time::microseconds(0); // One update per call.
  alps::Parameters p;
  p["ALGORITHM"] = "lifecycle"; p["SEED"] = 42;
  alps::params ngs_params;
  ngs_params["SEED"] = 42;

  for (int steps : {1,3}) { // Checkpoint before and after thermalization.
    for (auto format : {alps::dump_format::hdf5, alps::dump_format::xdr}) {
      opt.dump_format = format;
      std::string base = "classic-" + std::to_string(steps) + "-" + std::to_string(format);
      auto make = [&](bool fresh) -> std::unique_ptr<alps::abstract_clone> {
#ifdef ALPS_HAVE_MPI
        if (mpi) return std::make_unique<alps::clone_mpi>(world, world, directory, opt,
          alps::clone_create_msg_t(0,0,0,p,base,fresh));
#endif
        return std::make_unique<alps::clone>(directory,opt,0,0,p,base,fresh);
      };
      exercise(make, [](auto& clone) { clone.run(); }, [&](auto const& info) {
        std::vector<alps::ObservableSet> obs;
        if (format == alps::dump_format::hdf5) {
          alps::hdf5::archive ar((directory/info.dumpfile_h5()).string());
          bool loaded = mpi ? alps::load_observable(ar,0,rank,obs) : alps::load_observable(ar,0,obs);
          require(loaded, "missing checkpoint measurements");
        } else {
          alps::IXDRFileDump dump(directory/info.dumpfile_xdr());
          alps::load_observable(dump,obs);
        }
        alps::RealObsevaluator count(obs.at(0)["Count"]);
        require(count.count() == 4 && count.mean() == 4.5, "resume changed thermalized samples");
      }, leader, steps);
    }
    int stops=0, reports=0;
    boost::function<bool()> stop = [&] { ++stops; return false; };
    boost::function<void(double)> progress = [&](double) { ++reports; };
    std::string base = "ngs-" + std::to_string(steps);
    auto make = [&](bool fresh) -> std::unique_ptr<alps::ngs_parapack::abstract_clone> {
#ifdef ALPS_HAVE_MPI
      if (mpi) return std::make_unique<alps::ngs_parapack::clone_mpi>(world,world,directory,
        opt.dump_policy,opt.check_interval,alps::ngs_parapack::clone_create_msg_t(0,0,0,ngs_params,base,fresh));
#endif
      return std::make_unique<alps::ngs_parapack::clone>(directory,opt.dump_policy,opt.check_interval,
        0,0,ngs_params,base,fresh);
    };
    exercise(make, [&](auto& clone) { clone.run(stop,progress); }, [&](auto const& info) {
      alps::hdf5::archive ar((directory/info.dumpfile_h5()).string());
      int count=0;
      ar["/test/count"] >> count;
      require(count == 6 && stops == 6 && reports == 6, "resume changed worker state or callbacks");
    }, leader, steps);
  }
  p["INITIAL_COUNT"] = 6;
  ngs_params["INITIAL_COUNT"] = 6;
#ifdef ALPS_HAVE_MPI
  if (mpi) {
    alps::clone_mpi finished(world,world,directory,opt,
      alps::clone_create_msg_t(0,0,0,p,"complete-classic",true));
    alps::ngs_parapack::clone_mpi finished_ngs(world,world,directory,opt.dump_policy,opt.check_interval,
      alps::ngs_parapack::clone_create_msg_t(0,0,0,ngs_params,"complete-ngs",true));
    require(finished.halted() && finished_ngs.halted(), "already complete MPI workers did not halt");
  } else
#endif
  {
    alps::clone finished(directory,opt,0,0,p,"complete-classic",true);
    alps::ngs_parapack::clone finished_ngs(directory,opt.dump_policy,opt.check_interval,
      0,0,ngs_params,"complete-ngs",true);
    require(finished.halted() && finished_ngs.halted(), "already complete workers did not halt");
  }
#ifdef ALPS_HAVE_MPI
  if (mpi) world.barrier();
#endif
  if (leader) boost::filesystem::remove_all(directory);
}
