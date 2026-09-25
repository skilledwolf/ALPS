// SPDX-License-Identifier: MIT
#include <alps/parapack/exchange.h>
#ifdef ALPS_HAVE_MPI
#include <alps/parapack/exchange_multi.h>
#endif
#include <alps/osiris/xdrdump.h>
#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <iterator>
#include <stdexcept>

struct walker {
  using weight_parameter_type = double;
  explicit walker(alps::Parameters const& p) : energy(static_cast<unsigned>(p["WORKER_SEED"]) % 7) {}
#ifdef ALPS_HAVE_MPI
  walker(boost::mpi::communicator const&, alps::Parameters const& p) : walker(p) {}
#endif
  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Energy");
  }
  void set_beta(double b) { beta = b; }
  void run(alps::ObservableSet& obs) {
    energy = std::sin(energy + beta);
    obs["Energy"] << energy;
  }
  double weight_parameter() const { return energy; }
  static double log_weight(double energy, double beta) { return energy * beta; }
  void save(alps::ODump& dp) const { dp << energy << beta; }
  void load(alps::IDump& dp) { dp >> energy >> beta; }
  double energy, beta = 0;
};

template<class Worker>
std::string snapshot(Worker const& worker, boost::filesystem::path const& path) {
  { alps::OXDRFileDump dump(path); worker.save_worker(dump); }
  std::ifstream in(path.string(), std::ios::binary);
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

template<class Worker>
void check(Worker& uninterrupted, Worker& resumed, alps::Parameters const& p, int ranks_per_replica = 1) {
  std::vector<alps::ObservableSet> obs;
  uninterrupted.init_observables(p, obs);
  for (int i = 0; i < 32; ++i) uninterrupted.run(obs);
  auto path = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("alps-exchange-%%%%-%%%%");
  snapshot(uninterrupted, path);
  { alps::IXDRFileDump dump(path); resumed.load_worker(dump); }
  auto restored_obs = obs;
  for (int i = 0; i < 32; ++i) {
    uninterrupted.run(obs);
    resumed.run(restored_obs);
  }
  auto expected = snapshot(uninterrupted, path);
  auto actual = snapshot(resumed, path);
  boost::filesystem::remove(path);
  if (actual != expected) throw std::runtime_error("Exchange checkpoint changed the random stream or replica state");
  unsigned long long samples = 0;
  for (std::size_t p = 0; p < obs.size(); ++p) {
    alps::RealObsevaluator a(obs[p]["Energy"]), b(restored_obs[p]["Energy"]);
    samples += a.count();
    if (a.count() != b.count() || (a.count() && a.mean() != b.mean()))
      throw std::runtime_error("Exchange restart changed measurements");
  }
#ifdef ALPS_HAVE_MPI
  samples = boost::mpi::all_reduce(boost::mpi::communicator(), samples, std::plus<unsigned long long>());
#endif
  if (samples != 64 * static_cast<int>(p["NUM_REPLICAS"]) * ranks_per_replica)
    throw std::runtime_error("Exchange lost replicas while distributing workers");
  // Deterministic fingerprint also allows comparison with the pre-refactor implementation.
  uint64_t fingerprint = 14695981039346656037ull;
  for (unsigned char c : actual) fingerprint = (fingerprint ^ c) * 1099511628211ull;
  std::cout << "checkpoint " << fingerprint << '\n';
}

int main(int argc, char** argv) {
#ifdef ALPS_HAVE_MPI
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;
#endif
  for (int replicas : {4, 5}) for (int random : {0, 1})
  for (std::string mode : {"disabled", "none", "rate", "population"}) {
    alps::Parameters p;
    p["WORKER_SEED"] = 42;
    p["DISORDER_SEED"] = 7;
    p["BETA_MIN"] = 0.5;
    p["BETA_MAX"] = 2;
    p["NUM_REPLICAS"] = replicas;
    p["SWEEPS"] = 128;
    p["THERMALIZATION"] = 8;
    p["RANDOM_EXCHANGE"] = random;
    p["NO_EXCHANGE"] = mode == "disabled";
    p["OPTIMIZE_TEMPERATURE"] = mode == "rate" || mode == "population";
    p["OPTIMIZATION_TYPE"] = mode;
    p["INITIAL_BLOCK_SWEEPS"] = 8;
    p["OPTIMIZATION_ITERATIONS"] = 1;
    p["PROCESS_PER_WORKER"] = argc > 1 ? std::stoi(argv[1]) : 1;
#ifdef ALPS_HAVE_MPI
    alps::parapack::parallel_exchange_worker<walker> parallel(comm, p), parallel_copy(comm, p);
    check(parallel, parallel_copy, p);
    alps::parapack::multiple_parallel_exchange_worker<walker> multiple(comm, p), multiple_copy(comm, p);
    check(multiple, multiple_copy, p, static_cast<int>(p["PROCESS_PER_WORKER"]));
#else
    alps::parapack::single_exchange_worker<walker> serial(p), serial_copy(p);
    check(serial, serial_copy, p);
#endif
  }
}
