/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#pragma once
#include "kernel.hpp"
#include "ising_schema.hpp"
#include <alps/lattice.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/mc/driver.hpp>
#include <alps/mc/physical_moments.hpp>
#include <alps/alea/transformer.hpp>

template<class Kernel> class basic_ising {
public:
    struct statistics {
        alps::alea::batch_acc<double> joint;
        alps::alea::autocorr_acc<double> diagnostics;
        alps::mc::physical_moments physical;
        explicit statistics(size_t bins):joint(6,bins),diagnostics(6),physical(1,bins){}
        void add(std::vector<double> const& x) {
            physical.add(joint,alps::alea::column<double>{x[1]});
            joint<<alps::alea::make_adapter(x); diagnostics<<alps::alea::make_adapter(x);
        }
    };
    template<class... Context>
    basic_ising(alps::params const& p,size_t bins,size_t chain,Context const&... context):parameters_(p),bins_(bins),chain_(chain),
        random_(p["SEED"].as<int>()+chain,p["RNG"].as<std::string>()),
        state_(alps::graph_helper<>(alps::make_deprecated_parameters(p)),p["J"].as<double>(),[&]{return random_();},context...) {
        const bool scan=p["ALGORITHM"].as<std::string>()=="ising; temperature scan";
        const size_t n=scan ? p["NUM_TEMPERATURES"].as<size_t>() : 1;
        uint64_t offset=0;
        for (size_t i=0;i<n;++i) {
            const double t=scan ? p["INITIAL_TEMPERATURE"].as<double>()+i*p["DIFF_TEMPERATURE"].as<double>() : p["T"].as<double>();
            if (!(t>0) || !std::isfinite(t) || !std::isfinite(1/t)) throw std::invalid_argument("Temperatures and inverse temperatures must be finite and positive");
            temperatures_.push_back(t);
            const uint64_t warm=i==0 && scan ? p.value_or<uint64_t>("INITIAL_THERMALIZATION",p["THERMALIZATION"].as<uint64_t>()) : p["THERMALIZATION"].as<uint64_t>();
            const uint64_t production=p["SWEEPS"].as<uint64_t>();
            if (warm>UINT64_MAX-offset || production>UINT64_MAX-offset-warm) throw std::overflow_error("Temperature-scan sweep count overflow");
            start_.push_back(offset); warm_.push_back(warm); offset+=warm+production;
            stats_.emplace_back(bins);
        }
        total_=offset;
    }
    static alps::params checkpoint_parameters(alps::params p) {
        if (p["ALGORITHM"].as<std::string>()=="ising") p.erase("SWEEPS");
        return p;
    }
    uint64_t completed_sweeps() const {return steps_;}
    double fraction_completed() const {return double(steps_)/total_;}
    size_t stages() const {return stats_.size();}
    double temperature(size_t i) const {return temperatures_.at(i);}
    auto const& stage_statistics(size_t i) const {return stats_.at(i);}
    void update() {
        if (steps_>=total_) throw std::logic_error("Run is complete");
        const size_t i=std::upper_bound(start_.begin(),start_.end(),steps_)-start_.begin()-1;
        state_.step(1/temperatures_[i],[&]{return random_();});
        ++steps_;
        if (steps_>start_[i]+warm_[i]) stats_[i].add(state_.sample());
    }
    void measure() {}
    void synchronize() {state_.synchronize();}
    void save(alps::hdf5::archive& ar) const {
        ar["/parameters"]<<parameters_; ar["checkpoint/engine"]<<random_;
        ar["checkpoint/sweeps"]<<steps_; ar["checkpoint/chain"]<<uint64_t(chain_);
        state_.save(ar);
        for (size_t i=0;i<stages();++i) {
            alps::alea::hdf5_serializer codec(ar,"stages/"+std::to_string(i));
            serialize(codec,"joint",stats_[i].joint); serialize(codec,"autocorrelation",stats_[i].diagnostics);
            stats_[i].physical.save(ar,"stages/"+std::to_string(i)+"/physical");
        }
    }
    void load(alps::hdf5::archive& ar) {
        auto restored=*this;
        alps::params saved; uint64_t chain;
        ar["/parameters"]>>saved; ar["checkpoint/chain"]>>chain;
        ar["checkpoint/sweeps"]>>restored.steps_; ar["checkpoint/engine"]>>restored.random_;
        if (checkpoint_parameters(saved)!=checkpoint_parameters(parameters_) || chain!=chain_
                || restored.random_.name()!=random_.name() || restored.steps_>total_
                || ar.list_children("stages").size()!=stages()) throw std::invalid_argument("Checkpoint does not match this scan");
        restored.state_.load(ar);
        for (size_t i=0;i<stages();++i) {
            auto& values=restored.stats_[i];
            alps::alea::hdf5_serializer codec(ar,"stages/"+std::to_string(i));
            deserialize(codec,"joint",values.joint); deserialize(codec,"autocorrelation",values.diagnostics);
            const uint64_t count=restored.steps_<=start_[i]+warm_[i] ? 0 : std::min(restored.steps_-start_[i]-warm_[i],parameters_["SWEEPS"].as<uint64_t>());
            if (values.joint.size()!=6 || values.joint.num_batches()!=bins_ || values.joint.count()!=count
                    || values.joint.current_batch_size()!=values.joint.cursor().factor()
                    || !values.joint.store().batch().allFinite() || values.diagnostics.size()!=6
                    || values.diagnostics.count()!=count || values.diagnostics.batch_size()!=1 || values.diagnostics.granularity()!=2)
                throw std::invalid_argument("Invalid temperature-stage measurements");
            values.physical.load(ar,"stages/"+std::to_string(i)+"/physical",values.joint);
        }
        *this=std::move(restored);
    }
private:
    alps::params parameters_;
    size_t bins_,chain_;
    alps::random01 random_;
    Kernel state_;
    std::vector<statistics> stats_;
    std::vector<double> temperatures_;
    std::vector<uint64_t> start_,warm_;
    uint64_t steps_=0,total_=0;
};

using single_ising=basic_ising<ising_kernel>;

template<class Simulation,class Group=alps::mc::parallel>
int ising_main(int argc,char** argv,char const* command) {
    auto prepare=[](alps::params& p,alps::run_configuration const&) {
        if (!p.exists("THERMALIZATION")) p["THERMALIZATION"]=p["SWEEPS"].as<uint64_t>()/8;
        if (p["ALGORITHM"].as<std::string>()=="ising") {
            if (!p.exists("T")) throw std::invalid_argument("ising requires T");
            for (auto key:{"NUM_TEMPERATURES","INITIAL_TEMPERATURE","DIFF_TEMPERATURE","INITIAL_THERMALIZATION"})
                if (p.exists(key)) throw std::invalid_argument("Scan parameters require the temperature-scan algorithm");
        } else {
            if (p.exists("T")) throw std::invalid_argument("Scans define temperatures through their initial value and increment");
            if (!p.exists("NUM_TEMPERATURES") || !p.exists("INITIAL_TEMPERATURE") || !p.exists("DIFF_TEMPERATURE"))
                throw std::invalid_argument("Temperature scans require count, initial temperature and increment");
        }
    };
    auto publish=[](alps::run_configuration const& run,auto const& chains,alps::params const& p) {
        alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(),[&](auto& ar) {
            ar["/parameters"]<<p; ar["/run_config"]<<run;
            const auto n=chains.front()->stages();
            for (size_t i=0;i<n;++i) {
                std::vector<alps::alea::batch_result<double>> batches;
                alps::mc::moment_results moments;
                for (auto const& chain:chains) {
                    auto const& stats=chain->stage_statistics(i);
                    batches.push_back(stats.joint.result());
                    auto bins=stats.physical.results(); moments.insert(moments.end(),bins.begin(),bins.end());
                }
                const auto joint=alps::alea::merge(batches);
                const std::string path=n==1 ? "/simulation" : "/simulation/replicas/"+std::to_string(i);
                auto stage=p; stage["T"]=chains.front()->temperature(i);
                ar[path+"/parameters"]<<stage;
                if (n==1) ar["/parameters"]<<stage;
                alps::alea::hdf5_serializer raw(ar,path),output(ar,path+"/results");
                serialize(raw,"joint",joint); ar.create_group(path+"/results");
                std::map<std::string,alps::alea::batch_result<double>> results;
                alps::mc::unavailable_results unavailable;
                if (joint.observations()>1) {
                    const char* names[]={"Number of Sites","Energy","Energy^2","Magnetization","Magnetization^2","Magnetization^4"};
                    for (size_t k=0;k<6;++k) {
                        Eigen::Matrix<double,1,6> select=Eigen::Matrix<double,1,6>::Zero(); select(k)=1;
                        results.emplace(names[k],alps::alea::transform(alps::alea::jackknife_prop(),alps::alea::linear_transformer<double>(select),joint));
                    }
                    alps::mc::estimate(results,&unavailable,"Binder Ratio of Magnetization",joint,[](auto const& x){return x(5)>0 ? x(4)*x(4)/x(5) : NAN;});
                    const auto centered=alps::mc::centered_batches(moments).first;
                    const double beta=1/chains.front()->temperature(i),sites=joint.mean()(0);
                    alps::mc::estimate(results,&unavailable,"Specific Heat",centered,[=](auto const& x){return (beta*((x(1)-x(0)*x(0))*beta))/sites;});
                } else unavailable["Statistics"]="At least two effective batches are required";
                for (auto const& [name,value]:results) serialize(output,ar.encode_segment(name),value);
                for (auto const& [name,reason]:unavailable) ar[path+"/unavailable/"+ar.encode_segment(name)]<<reason;
                for (size_t id=0;id<chains.size();++id) {
                    alps::alea::hdf5_serializer diag(ar,path+"/realizations/0/clones/"+std::to_string(id)+"/autocorrelation");
                    serialize(diag,"Moments",chains[id]->stage_statistics(i).diagnostics.result());
                }
            }
        });
    };
    return alps::mc::main<Simulation,Group>(argc,argv,command,ising_schema,{},prepare,publish);
}
