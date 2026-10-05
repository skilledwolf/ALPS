// SPDX-License-Identifier: MIT
// Exercise the actual spatial worker, independently of its historical goldens.
#include "../../../../tutorials/00-examples/parapack/multiple/ising.h"
#include <boost/mpi/environment.hpp>
#include <array>
#include <cstdio>
#include <iostream>

// Make application-level blocking sends synchronous: halo correctness must not
// depend on MPI's eager message buffers. The old send-then-receive ring hangs.
extern "C" int MPI_Send(const void* buffer, int count, MPI_Datatype type, int peer,
                         int tag, MPI_Comm communicator) {
  return PMPI_Ssend(buffer, count, type, peer, tag, communicator);
}

namespace {
class memory_output : public alps::OXDRDump {
public:
  memory_output() : file_(std::tmpfile()) {
    if (!file_) throw std::runtime_error("Cannot create temporary checkpoint stream");
    xdrstdio_create(&xdr_, file_, XDR_ENCODE);
  }
  ~memory_output() { xdr_destroy(&xdr_); std::fclose(file_); }
  std::vector<char> bytes() const {
    std::fflush(file_); const auto size = std::ftell(file_); std::rewind(file_);
    std::vector<char> result(size);
    if (std::fread(result.data(), 1, result.size(), file_) != result.size())
      throw std::runtime_error("Cannot read temporary checkpoint stream");
    return result;
  }
private:
  FILE* file_;
};
class memory_input : public alps::IXDRDump {
public:
  explicit memory_input(std::vector<char> const& data) : file_(std::tmpfile()) {
    if (!file_) throw std::runtime_error("Cannot create temporary checkpoint stream");
    if (std::fwrite(data.data(), 1, data.size(), file_) != data.size())
      throw std::runtime_error("Cannot write temporary checkpoint stream");
    std::rewind(file_); xdrstdio_create(&xdr_, file_, XDR_DECODE);
  }
  ~memory_input() { xdr_destroy(&xdr_); std::fclose(file_); }
private:
  FILE* file_;
};
void require(mpi::communicator const& comm, bool valid, char const* reason) {
  if (!mpi::all_reduce(comm, valid, std::logical_and<bool>())) throw std::runtime_error(reason);
}
alps::Parameters parameters(mpi::communicator const& comm, int length, double coupling = 1) {
  alps::Parameters p;
  p["L"] = length; p["J"] = coupling; p["T"] = 2.;
  p["SWEEPS"] = 100000; p["THERMALIZATION"] = 2000;
  // Independent local streams, as supplied by the production clone scheduler.
  p["WORKER_SEED"] = 2873 + 1009 * comm.rank(); p["DISORDER_SEED"] = 2873;
  return p;
}
std::array<double, 5> exact(int length, double coupling) {
  std::array<double, 5> sum{}; double z = 0;
  for (unsigned state = 0; state < (1u << length); ++state) {
    double energy = 0, mag = 0;
    for (int i = 0; i < length; ++i) {
      const int spin = 2 * int((state >> i) & 1) - 1;
      const int next = 2 * int((state >> ((i + 1) % length)) & 1) - 1;
      energy -= coupling * spin * next; mag += spin;
    }
    const double weight = std::exp(-energy / 2.);
    const std::array<double, 5> x{energy, energy * energy, mag, mag * mag, mag * mag * mag * mag};
    z += weight;
    for (size_t i = 0; i < sum.size(); ++i) sum[i] += weight * x[i];
  }
  for (auto& x : sum) x /= z;
  return sum;
}
void physical_checkpoint(mpi::communicator const& comm, parallel_ising_worker const& worker,
                         int length, double coupling) {
  memory_output out; worker.save(out); auto bytes = out.bytes(); memory_input in(bytes);
  alps::mc_steps steps; tabbed_vector<int> spins; double recorded_energy;
  in >> steps >> spins >> recorded_energy;
  bool valid = std::all_of(&spins[-1], &spins[int(spins.size())] + 1, [](int s) { return s == 0 || s == 1; });
  std::vector<int> local(&spins[0], &spins[int(spins.size())]), global;
  std::vector<std::vector<int>> blocks; mpi::all_gather(comm, local, blocks);
  for (auto const& block : blocks) global.insert(global.end(), block.begin(), block.end());
  valid = valid && global.size() == size_t(length);
  const auto& previous = blocks[(comm.rank() + comm.size() - 1) % comm.size()];
  const auto& next = blocks[(comm.rank() + 1) % comm.size()];
  valid = valid && spins[-1] == previous.back() && spins[spins.size()] == next.front();
  double energy = 0;
  for (size_t i = 0; i < global.size(); ++i)
    energy -= coupling * (2 * global[i] - 1) * (2 * global[(i + 1) % global.size()] - 1);
  if (comm.rank() == 0) valid = valid && energy == recorded_energy;
  require(comm, valid, "Saved physical configuration, halos or energy are inconsistent");
}
void thermodynamics(mpi::communicator const& comm, int length, double coupling) {
  auto p = parameters(comm, length, coupling);
  parallel_ising_worker worker(comm, p); alps::ObservableSet obs; worker.init_observables(p, obs);
  for (int i = 0; i < 102000; ++i) {
    worker.run(obs);
    if (i == 1999) obs.reset();
    if (comm.rank() == 0 && std::abs(worker.weight_parameter()) > std::abs(coupling) * length)
      throw std::runtime_error("Energy is outside physical bounds");
    if (i % 10000 == 0) physical_checkpoint(comm, worker, length, coupling);
  }
  bool valid = true;
  if (comm.rank() == 0) {
    const char* names[]{"Energy", "Energy^2", "Magnetization", "Magnetization^2", "Magnetization^4"};
    auto expected = exact(length, coupling);
    for (size_t i = 0; i < expected.size(); ++i) {
      alps::RealObsevaluator result = obs[names[i]];
      const auto count = dynamic_cast<alps::RealObservable const&>(obs[names[i]]).count();
      const double tolerance = std::max(.025 * std::max(1., std::abs(expected[i])), 6 * result.error());
      // Legacy evaluators discard incomplete bins; the live accumulator counts
      // every production update. Native ALEA migration removes that truncation.
      if (count != 100000 || std::abs(result.mean() - expected[i]) > tolerance) {
        std::cerr << "L=" << length << " J=" << coupling << " " << names[i] << ": "
                  << result.mean() << " expected " << expected[i] << " tolerance " << tolerance
                  << " count " << count << '\n';
        valid = false;
      }
    }
  }
  require(comm, valid, "Spatial-worker moments disagree with exact canonical enumeration");
}
std::vector<char> snapshot(parallel_ising_worker const& worker, alps::ObservableSet const& obs) {
  memory_output out; worker.save_worker(out); out << obs; return out.bytes();
}
void restart(mpi::communicator const& comm, int length, char const* rng) {
  auto p = parameters(comm, length); p["RNG"] = rng;
  parallel_ising_worker full(comm, p); alps::ObservableSet full_obs; full.init_observables(p, full_obs);
  for (int i = 0; i < 79; ++i) full.run(full_obs);
  parallel_ising_worker part(comm, p); alps::ObservableSet part_obs; part.init_observables(p, part_obs);
  for (int i = 0; i < 31; ++i) part.run(part_obs);
  auto checkpoint = snapshot(part, part_obs);
  parallel_ising_worker resumed(comm, p); alps::ObservableSet resumed_obs; resumed.init_observables(p, resumed_obs);
  memory_input in(checkpoint); resumed.load_worker(in); in >> resumed_obs;
  for (int i = 31; i < 79; ++i) resumed.run(resumed_obs);
  require(comm, snapshot(full, full_obs) == snapshot(resumed, resumed_obs), "Spatial checkpoint did not resume exactly");
}
void invalid_input(mpi::communicator const& comm, int length) {
  for (int fault = 0; fault < 7; ++fault) {
    auto p = parameters(comm, length);
    if (fault == 0) p["L"] = 2 * comm.size() - 1; // Old uneven-partition failure could hang rank zero.
    if (comm.rank() == 0) {
      if (fault == 1) p["T"] = 0;
      if (fault == 2) p["L"] = 2.5;
      if (fault == 3) p["J"] = "1/0";
      if (comm.size() > 1 && fault == 4) p["J"] = 2.;
      if (comm.size() > 1 && fault == 5) p["L"] = length + 1;
      if (comm.size() > 1 && fault == 6) p["T"] = 3.;
    }
    if (comm.size() == 1 && fault >= 4) continue;
    bool rejected = false;
    try { parallel_ising_worker worker(comm, p); } catch (std::invalid_argument const&) { rejected = true; }
    require(comm, rejected, "Invalid input was not rejected on every spatial rank");
  }
  auto p = parameters(comm, length); p.erase("T");
  parallel_ising_worker worker(comm, p); alps::ObservableSet obs; worker.init_observables(p, obs);
  if (comm.rank() == 0) worker.set_beta(.5);
  bool rejected = false;
  try { worker.run(obs); } catch (std::invalid_argument const&) { rejected = true; }
  // At one rank, beta was set and the update is valid; otherwise unset ranks force rejection everywhere.
  require(comm, rejected == (comm.size() > 1), "Unset beta did not reach consensus");
  for (double beta : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    rejected = false;
    worker.set_beta(.5);
    if (comm.rank() == 0) worker.set_beta(beta);
    try { worker.run(obs); } catch (std::invalid_argument const&) { rejected = true; }
    require(comm, rejected, "Invalid exchange beta was accepted");
  }
  if (comm.size() > 1) {
    worker.set_beta(comm.rank() == 0 ? 1 : .5); rejected = false;
    try { worker.run(obs); } catch (std::invalid_argument const&) { rejected = true; }
    require(comm, rejected, "Different valid inverse temperatures were accepted");
  }
  worker.set_beta(0); worker.run(obs); // Infinite-temperature exchange endpoint remains supported.
  p["J"] = 0.; parallel_ising_worker uncoupled(comm, p); alps::ObservableSet zero_obs;
  uncoupled.init_observables(p, zero_obs); uncoupled.set_beta(1e300);
  for (int i = 0; i < 19; ++i) uncoupled.run(zero_obs);
  require(comm, uncoupled.weight_parameter() == 0, "Uncoupled spins acquired energy at extreme beta");
}
}
int main(int argc, char** argv) {
  mpi::environment environment(argc, argv); mpi::communicator world;
  try {
    const int length = std::max(4, 2 * world.size()) + 1;
    invalid_input(world, length);
    restart(world, length, "mt19937"); restart(world, length, "lagged_fibonacci607");
    for (double coupling : {1., -1., 0.}) thermodynamics(world, length, coupling);
    if (world.rank() == 0) std::cout << "Spatial Ising checks passed at " << world.size() << " ranks\n";
  } catch (std::exception const& error) {
    std::cerr << "rank " << world.rank() << ": " << error.what() << '\n'; world.abort(1);
  }
}
