// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/mc/driver.hpp>
#include <csignal>
#include <iostream>

namespace {
void require(bool condition,char const* message) {
    int success=condition;
    MPI_Allreduce(MPI_IN_PLACE,&success,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if (!success) throw std::runtime_error(message);
}
void clear_signals() {
    alps::ngs::signal pending;
    while (!pending.empty()) pending.pop();
}

class simulation {
    size_t id_;
    uint64_t target_,steps_=0;
    boost::mpi::communicator const* communicator_;
public:
    simulation(alps::params const& p,size_t,size_t id,boost::mpi::communicator const* communicator=nullptr)
      :id_(id),target_(p["SWEEPS"].as<uint64_t>()),communicator_(communicator) {}
    static alps::params checkpoint_parameters(alps::params p) {return p;}
    uint64_t completed_sweeps() const {return steps_;}
    double fraction_completed() const {return double(steps_)/target_;}
    void update() {
        if (communicator_) communicator_->barrier();
        ++steps_;
        // A stop arriving during a sweep is local to rank 1. Every rank must
        // agree before the next collective update; root has no local signal.
        if (steps_==1 && (communicator_ ? communicator_->rank()==1 : id_==1))
            alps::ngs::signal::slot(SIGTERM);
    }
    void measure() {}
    void save(alps::hdf5::archive& ar) const {ar["steps"]<<steps_;}
    void load(alps::hdf5::archive& ar) {ar["steps"]>>steps_;}
};

struct collective_group : alps::mc::collective {
    template<class Simulation>
    auto make(alps::params const& p,size_t bins,size_t id) const {
        return std::make_unique<Simulation>(p,bins,id,&world);
    }
    // The test's physical counter is already identical on every participant.
    template<class Chains> void synchronize(Chains&) const {}
};

template<class Group> void contract(Group const& group,size_t count) {
    clear_signals();
    alps::run_configuration run;
    run.parameters["SWEEPS"]=64;run.parameters["THERMALIZATION"]=0;
    run.execution["seed"]=17;run.execution["rng"]="mt19937";
    run.execution["chains"]=count;run.execution["bins"]=8;
    run.execution["time_limit"]=0.;run.execution["checkpoint_interval"]=0.;
    run.execution["max_sweeps"]=0;
    auto prepare=[](alps::params&,alps::run_configuration const&) {};
    auto chains=alps::mc::prepare_chains<simulation>(run,prepare,group);
    bool published=false;
    auto publish=[&](auto const&,auto const& state,auto const&) {
        published=true;
        if (state.front()->completed_sweeps()!=(count==1 ? 1 : 32) ||
            (count==2 && state.back()->completed_sweeps()!=1))
            throw std::runtime_error("Driver stop changed physical sweep counts");
    };
    alps::mc::execute(run,chains,prepare,publish,alps::mc::snapshot_type<simulation>{},group);
    alps::ngs::signal pending;
    require(pending.empty()==(group.rank()==0),"Stop was not confined to rank 1");
    require(published==(group.rank()==0),"Driver published on the wrong rank");
    const auto owned=count==1 ? 0 : size_t(group.rank());
    require(chains[owned]->completed_sweeps()==(count==1 || group.rank()==1 ? 1 : 32),
            "Driver entered an extra collective sweep or changed independent scheduling");
    clear_signals();
}
}

int main(int argc,char** argv) {
    boost::mpi::environment environment(argc,argv);
    try {
        collective_group collective;
        require(collective.size()==2,"Driver stop contract requires two MPI ranks");
        contract(collective,1);
        contract(alps::mc::parallel{},2);
        if (!collective.rank()) std::cout<<"Native driver MPI stop contracts passed\n";
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
