// SPDX-License-Identifier: MIT
#include <alps/parapack/parapack.h>
#include <alps/ngs/parapack/parapack.h>
#include <alps/hdf5.hpp>
#include <string>

class classic_worker : public alps::parapack::abstract_worker {
public:
  explicit classic_worker(alps::Parameters const&) {}
  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) override {
    obs << alps::RealObservable("Count");
  }
  void run(alps::ObservableSet& obs) override { obs["Count"] << double(++count_); }
  void load(alps::IDump& dump) override { dump >> count_; }
  void save(alps::ODump& dump) const override { dump << count_; }
  bool is_thermalized() const override { return true; }
  double progress() const override { return count_ / 4.; }
private:
  int count_ = 0;
};

class ngs_worker {
public:
  explicit ngs_worker(alps::params const&) {}
  bool run(boost::function<bool()> const&, boost::function<void(double)> const&) {
    ++count_;
    return fraction_completed() >= 1;
  }
  void load(alps::hdf5::archive& ar) { ar["/test/count"] >> count_; }
  void save(alps::hdf5::archive& ar) const { ar["/test/count"] << count_; }
  double fraction_completed() const { return count_ / 4.; }
private:
  int count_ = 0;
};

int main(int argc, char** argv) {
  if (argc < 2) return 2;
  bool ngs = std::string(argv[1]) == "ngs";
  --argc;
  ++argv;
  int result;
  if (ngs) {
    result = alps::ngs_parapack::start<ngs_worker>(argc, argv);
  } else {
    alps::parapack::worker_factory::instance()->register_worker<classic_worker>("test");
    result = alps::parapack::start(argc, argv);
  }
  // Check the API result before the OS maps negative process exit codes.
  if (argc == 2 && (std::string(argv[1]) == "--invalid-option" || std::string(argv[1]) == "missing.xml"))
    return result == (ngs ? -1 : 127) ? 0 : 1;
  return result;
}
