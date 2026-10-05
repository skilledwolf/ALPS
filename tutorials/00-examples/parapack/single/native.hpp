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

class single_ising {
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
    single_ising(alps::params const& p,size_t bins,size_t chain):parameters_(p),bins_(bins),chain_(chain),
        random_(p["SEED"].as<int>()+chain,p["RNG"].as<std::string>()),
        state_(alps::graph_helper<>(alps::make_deprecated_parameters(p)),p["J"].as<double>(),[&]{return random_();}) {
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
        single_ising restored(parameters_,bins_,chain_);
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
    ising_kernel state_;
    std::vector<statistics> stats_;
    std::vector<double> temperatures_;
    std::vector<uint64_t> start_,warm_;
    uint64_t steps_=0,total_=0;
};
