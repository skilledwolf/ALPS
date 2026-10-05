/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#pragma once
#include "ising_schema.hpp"
#include <alps/mc/driver.hpp>
#include <alps/alea/transform.hpp>
#include <alps/alea/transformer.hpp>
#include <alps/lattice.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/hdf5/vector.hpp>
#include <cmath>

// One sweep makes N random-site Metropolis proposals, with replacement.
// A symmetric proposal can keep the spin, so measurement strides remain ergodic.
class IsingSimulation : public alps::mcbase {
public:
    IsingSimulation(alps::params const& p, size_t bins, size_t chain)
        : mcbase(p,chain), chain_(chain), bins_(bins),
          warmup_(p["THERMALIZATION"].as<uint64_t>()), target_(p["SWEEPS"].as<uint64_t>()),
          beta_(1/p["T"].as<double>()), correlations_(p["CORRELATIONS"].as<bool>()) {
        if (!(p["T"].as<double>()>0) || target_>UINT64_MAX-warmup_)
            throw std::invalid_argument("Positive T and a representable sweep count are required");
        if (correlations_) {
            const auto n = p["L"].as<size_t>();
            neighbors_.resize(n);
            for (size_t i=0; i<n; ++i) neighbors_[i] = {(i+1)%n,(i+n-1)%n};
        } else {
            alps::graph_helper<> graph(alps::make_deprecated_parameters(p));
            if (graph.inhomogeneous()) throw std::invalid_argument("Disordered lattices are not supported");
            neighbors_.resize(graph.num_sites());
            for (size_t i=0; i<neighbors_.size(); ++i)
                for (auto [it,end]=graph.neighbors(i); it!=end; ++it) neighbors_[i].push_back(*it);
        }
        if (neighbors_.size()<2) throw std::invalid_argument("At least two sites are required");
        spins_.resize(neighbors_.size());
        for (auto& spin : spins_) spin = random()<.5 ? 1 : -1;
        const size_t dimension = correlations_ ? 4+spins_.size() : 2;
        measurements.emplace("Moments",std::make_shared<alps::alea::batch_acc<double>>(dimension,bins));
        measurements.emplace("Autocorrelation",std::make_shared<alps::alea::autocorr_acc<double>>(dimension));
    }
    static alps::params checkpoint_parameters(alps::params p) { p.erase("SWEEPS"); return p; }
    uint64_t completed_sweeps() const { return sweeps_; }
    uint64_t measurement_count() const { return measurement("Moments")->count(); }
    double fraction_completed() const override { return sweeps_<=warmup_ ? 0 : double(sweeps_-warmup_)/target_; }
    void update() override {
        if (sweeps_>=warmup_+target_ || measurement_count()!=count())
            throw std::logic_error("Complete each update/measure pair before advancing");
        for (size_t j=0; j<spins_.size(); ++j) {
            const size_t i = size_t(random()*spins_.size());
            if (random()<.5) continue; // propose the current spin with probability 1/2
            double delta = 0;
            for (auto neighbor : neighbors_[i])
                if (neighbor!=i) delta += 2.*spins_[i]*spins_[neighbor];
            if (delta<=0 || random()<std::exp(-beta_*delta)) spins_[i] = -spins_[i];
        }
        ++sweeps_;
    }
    void measure() override {
        if (sweeps_<=warmup_) return;
        if (measurement_count()!=count()-1) throw std::logic_error("Measure each update once");
        std::vector<double> values(correlations_ ? 4+spins_.size() : 2,0.);
        const double n = spins_.size();
        for (size_t i=0; i<spins_.size(); ++i) {
            values[1] += spins_[i]/n;
            for (auto neighbor : neighbors_[i]) values[0] -= spins_[i]*spins_[neighbor]/(2*n);
            if (correlations_) for (size_t d=0; d<spins_.size(); ++d)
                values[4+d] += spins_[i]*spins_[(i+d)%spins_.size()]/n;
        }
        if (correlations_) { values[2]=values[1]*values[1]; values[3]=values[2]*values[2]; }
        *measurement("Moments") << alps::alea::make_adapter(values);
        *measurement<alps::alea::autocorr_acc<double>>("Autocorrelation") << alps::alea::make_adapter(values);
    }
    alps::mc::batch_results collect_results() const { return {{"Moments",measurement("Moments")->result()}}; }
    void save(alps::hdf5::archive& ar) const override {
        if (measurement_count()!=count()) throw std::logic_error("Measure before checkpointing");
        mcbase::save(ar);
        ar["checkpoint/spins"] << spins_; ar["checkpoint/topology"] << topology();
        ar["checkpoint/sweeps"] << sweeps_; ar["checkpoint/chain"] << uint64_t(chain_);
    }
    void load(alps::hdf5::archive& ar) override {
        alps::params saved;
        uint64_t sweeps,chain;
        std::vector<int> spins;
        std::vector<uint64_t> graph;
        ar["/parameters"] >> saved;
        ar["checkpoint/spins"] >> spins; ar["checkpoint/topology"] >> graph;
        ar["checkpoint/sweeps"] >> sweeps; ar["checkpoint/chain"] >> chain;
        if (checkpoint_parameters(saved)!=checkpoint_parameters(parameters) || chain!=chain_
                || sweeps>warmup_+target_ || graph!=topology() || spins.size()!=spins_.size()
                || !std::all_of(spins.begin(),spins.end(),[](int s){return s==1 || s==-1;}))
            throw std::invalid_argument("Checkpoint physical state does not match this run");
        alps::mc::validate_measurements(measurements,ar,sweeps>warmup_ ? sweeps-warmup_ : 0,bins_);
        auto invocation = parameters;
        mcbase::load(ar); parameters=std::move(invocation);
        spins_=std::move(spins); sweeps_=sweeps;
    }
private:
    uint64_t count() const { return sweeps_>warmup_ ? sweeps_-warmup_ : 0; }
    std::vector<uint64_t> topology() const {
        std::vector<uint64_t> data;
        for (auto const& neighbors : neighbors_) {
            data.push_back(neighbors.size()); data.insert(data.end(),neighbors.begin(),neighbors.end());
        }
        return data;
    }
    size_t chain_,bins_;
    uint64_t warmup_,target_,sweeps_=0;
    double beta_;
    bool correlations_;
    std::vector<int> spins_;
    std::vector<std::vector<size_t>> neighbors_;
};
