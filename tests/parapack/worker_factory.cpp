// SPDX-License-Identifier: MIT
#include <alps/parapack/worker_factory.h>
#include <sstream>
#include <stdexcept>

void require(bool condition) {
  if (!condition) throw std::runtime_error("worker factory selection failed");
}

template<int Id>
struct worker : alps::parapack::abstract_worker {
  explicit worker(alps::Parameters const&) {}
#ifdef ALPS_HAVE_MPI
  worker(boost::mpi::communicator const&, alps::Parameters const&) {}
#endif
  void load(alps::IDump&) override {}
  void save(alps::ODump&) const override {}
  bool is_thermalized() const override { return true; }
  double progress() const override { return Id; }
};

struct capture {
  std::ostringstream out, err;
  std::streambuf* old_out = std::cout.rdbuf(out.rdbuf());
  std::streambuf* old_err = std::cerr.rdbuf(err.rdbuf());
  ~capture() { std::cout.rdbuf(old_out); std::cerr.rdbuf(old_err); }
};

template<class Factory, class Create>
void check(Factory* factory, Create create, std::string const& name) {
  capture messages;
  alps::Parameters params;
  auto failure = [&](std::string const& expected) {
    messages.err.str("");
    try {
      create(params);
    } catch (std::runtime_error const& error) {
      require(error.what() == name + "::make_creator()");
      require(messages.err.str().find(expected) != std::string::npos);
      return;
    }
    throw std::runtime_error("missing factory selection error");
  };
  failure("no algorithm registered");
  require(factory->template register_worker<worker<1>>("one"));
  require(create(params)->progress() == 1);
  params["ALGORITHM"] = "unknown";
  require(create(params)->progress() == 1);
  require(messages.out.str().find("The only algorithm \"one\"") != std::string::npos);
  require(factory->template register_worker<worker<2>>("two"));
  failure("unknown algorithm: \"unknown\" (registered algorithms: \"one\", \"two\")");
  params = alps::Parameters();
  failure("no algorithm specified (registered algorithms: \"one\", \"two\")");
  params["WORKER"] = "two";
  require(create(params)->progress() == 2);
  require(messages.out.str().find("parameter WORKER is obsolete") != std::string::npos);
  messages.out.str("");
  params["ALGORITHM"] = "one";
  require(create(params)->progress() == 1 && messages.out.str().empty());
  require(!factory->template register_worker<worker<2>>("one"));
  require(create(params)->progress() == 2);
  require(factory->unregister_worker("one"));
  require(!factory->unregister_worker("one"));
  require(factory->unregister_worker("two"));
  failure("no algorithm registered");
}

int main(int argc, char** argv) {
  using namespace alps::parapack;
#ifdef ALPS_HAVE_MPI
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator comm;
  check(parallel_worker_factory::instance(), [&](alps::Parameters const& params) {
    return parallel_worker_factory::make_worker(comm, params);
  }, "parallel_worker_factory");
#endif
  check(worker_factory::instance(), [](alps::Parameters const& params) {
    return worker_factory::make_worker(params);
  }, "worker_factory");
}
