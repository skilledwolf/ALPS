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
#include <alps/hdf5/vector.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// One physical heat-bath kernel shared by the native command and the remaining
// legacy exchange callers. Colors contain no connected pair of distinct sites.
class ising_kernel {
public:
    template<class Graph, class Random>
    ising_kernel(Graph const& graph, double coupling, Random&& random) : coupling_(coupling) {
        if (!graph.num_sites()) throw std::invalid_argument("The lattice must contain sites");
        neighbors_.resize(graph.num_sites());
        std::vector<size_t> color(neighbors_.size(),neighbors_.size());
        for (size_t i=0; i<neighbors_.size(); ++i) {
            std::vector<bool> used(colors_.size(),false);
            for (auto [it,end]=graph.neighbors(i); it!=end; ++it) {
                const size_t neighbor=*it;
                neighbors_[i].push_back(neighbor);
                if (neighbor!=i && color[neighbor]<used.size()) used[color[neighbor]]=true;
            }
            const size_t c=std::find(used.begin(),used.end(),false)-used.begin();
            if (c==colors_.size()) colors_.emplace_back();
            color[i]=c; colors_[c].push_back(i);
        }
        spins_.resize(neighbors_.size()); proposals_.resize(spins_.size());
        for (auto& spin:spins_) spin=random()<.5 ? 1 : -1;
    }
    template<class Random> void step(double beta, Random&& random) {
        for (auto const& sites:colors_) {
            // Consume RNG serially in a fixed order, independently of thread count.
            for (auto i:sites) proposals_[i]=random();
#ifdef _OPENMP
#pragma omp parallel for if(sites.size()>=256)
#endif
            for (size_t k=0; k<sites.size(); ++k) {
                const size_t i=sites[k];
                double field=0;
                for (auto j:neighbors_[i]) if (j!=i) field+=coupling_*spins_[j];
                const double x=(2*spins_[i]*field)*beta;
                // Logistic heat-bath flip probability, stable at both extremes.
                const double e=std::exp(-std::abs(x));
                if (proposals_[i]<(x>0 ? e/(1+e) : 1/(1+e))) spins_[i]=-spins_[i];
            }
        }
    }
    size_t size() const { return spins_.size(); }
    double energy() const {
        double result=0;
        for (size_t i=0; i<size(); ++i) for (auto j:neighbors_[i]) result-=coupling_*spins_[i]*spins_[j]/2;
        return result;
    }
    std::vector<double> sample() const {
        const double e=energy(),m=std::accumulate(spins_.begin(),spins_.end(),0.);
        return {double(size()),e,e*e,m,m*m,m*m*m*m};
    }
    auto const& spins() const { return spins_; }
    void restore(std::vector<int> spins) {
        if (spins.size()!=size() || !std::all_of(spins.begin(),spins.end(),[](int s){return s==1 || s==-1;}))
            throw std::invalid_argument("Invalid Ising checkpoint spins");
        spins_=std::move(spins);
    }
    std::vector<uint64_t> topology() const {
        std::vector<uint64_t> out;
        for (auto const& neighbors:neighbors_) {
            out.push_back(neighbors.size()); out.insert(out.end(),neighbors.begin(),neighbors.end());
        }
        return out;
    }
    void save(alps::hdf5::archive& ar) const {
        ar["checkpoint/spins"]<<spins_; ar["checkpoint/topology"]<<topology();
    }
    void load(alps::hdf5::archive& ar) {
        std::vector<int> spins; std::vector<uint64_t> graph;
        ar["checkpoint/spins"]>>spins; ar["checkpoint/topology"]>>graph;
        if (graph!=topology()) throw std::invalid_argument("Checkpoint lattice changed");
        restore(std::move(spins));
    }
private:
    double coupling_;
    std::vector<int> spins_;
    std::vector<double> proposals_;
    std::vector<std::vector<size_t>> neighbors_,colors_;
};
