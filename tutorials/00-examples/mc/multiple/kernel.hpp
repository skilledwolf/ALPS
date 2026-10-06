// SPDX-License-Identifier: MIT
#pragma once
#include "../single/kernel.hpp"
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/communicator.hpp>
#endif

// One periodic chain, split into contiguous blocks. The serial fallback keeps
// arbitrary graph support; spatial MPI requires a naturally ordered ring.
class spatial_ising_kernel {
public:
    template<class Graph,class Random>
    spatial_ising_kernel(Graph const& graph,double coupling,Random&& random):coupling_(coupling),size_(graph.num_sites()) {
        if (!size_ || size_>size_t(std::numeric_limits<int>::max()-2) || !std::isfinite(coupling_))
            throw std::invalid_argument("Spatial Ising requires a nonempty integer-sized lattice and finite J");
        bool ring=size_>=2;
        for (size_t i=0;i<size_;++i) {
            std::vector<size_t> neighbors;
            for (auto [it,end]=graph.neighbors(i);it!=end;++it) neighbors.push_back(*it);
            topology_.push_back(neighbors.size()); topology_.insert(topology_.end(),neighbors.begin(),neighbors.end());
            const size_t previous=(i+size_-1)%size_,next=(i+1)%size_;
            ring=ring && neighbors.size()==2 &&
                ((neighbors[0]==previous && neighbors[1]==next) || (neighbors[1]==previous && neighbors[0]==next));
        }
        if (!ring) {
            fallback_.emplace(graph,coupling_,std::forward<Random>(random));
            topology_=std::vector<uint64_t>{};
            return;
        }
        // Every rank consumes this same initialization stream. Subsequent
        // proposals come exclusively from the worker communicator's root.
        global_.resize(size_);
        for (auto& spin:global_) spin=random()<.5 ? 1 : -1;
        partition(1);
    }
    void validate_partition(int ranks) const {
        if (ranks<1) throw std::invalid_argument("Spatial Ising requires a positive rank count");
        if (ranks>1 && (fallback_ || size_/size_t(ranks)<2))
            throw std::invalid_argument("Spatial MPI Ising requires a periodic nearest-neighbor ring with at least two sites per rank");
    }
#ifdef ALPS_HAVE_MPI
    template<class Graph,class Random>
    spatial_ising_kernel(Graph const& graph,double coupling,Random&& random,boost::mpi::communicator const& comm)
        :spatial_ising_kernel(graph,coupling,std::forward<Random>(random)) {
        distribute(comm);
    }
    // Repartition the local full state without communication or random draws.
    // Nested runners can construct and load before binding a worker subgroup.
    void distribute(boost::mpi::communicator const& comm) {
        validate_partition(comm.size());
        if (communicator_) {
            if (MPI_Comm(*communicator_)!=MPI_Comm(comm))
                throw std::logic_error("Spatial Ising is already bound to another communicator");
            return;
        }
        if (!fallback_) global_.assign(spins_.begin()+1,spins_.end()-1);
        rank_=comm.rank(); communicator_=comm; // Share the handle; no communicator duplication.
        if (!fallback_) {
            partition(comm.size());
            if (rank_!=0) global_=std::vector<int>{};
        }
        synchronized_=true;
    }
#endif
    template<class Random> void step(double beta,Random&& random) {
        if (fallback_) {fallback_->step(beta,std::forward<Random>(random));return;}
        for (auto const& phase:phases_) {
            if (rank_==0)
                for (int i=0;i<phase.offsets.back()+phase.counts.back();++i) global_proposals_[i]=random();
#ifdef ALPS_HAVE_MPI
            if (counts_.size()>1) {
                BOOST_MPI_CHECK_RESULT(MPI_Scatterv,(rank_==0 ? global_proposals_.data() : nullptr,
                    phase.counts.data(),phase.offsets.data(),MPI_DOUBLE,proposals_.data(),int(phase.sites.size()),MPI_DOUBLE,0,*communicator_));
            } else
#endif
                std::copy_n(global_proposals_.begin(),phase.sites.size(),proposals_.begin());
#ifdef _OPENMP
#pragma omp parallel for if(phase.sites.size()>=256)
#endif
            for (size_t k=0;k<phase.sites.size();++k) {
                const size_t i=phase.sites[k]+1;
                const double difference=2*spins_[i]*(coupling_*(spins_[i-1]+spins_[i+1]));
                if (ising_kernel::heatbath_flip(difference,beta,proposals_[k])) spins_[i]=-spins_[i];
            }
            halos();
        }
        synchronized_=false;
    }
    size_t size() const {return size_;}
    std::vector<double> sample() const {
        if (fallback_) return fallback_->sample();
        std::int64_t local[2]{},global[2]{};
        for (size_t i=1;i+1<spins_.size();++i) {
            local[0]+=spins_[i]*spins_[i+1]; local[1]+=spins_[i];
        }
#ifdef ALPS_HAVE_MPI
        if (counts_.size()>1) {
            BOOST_MPI_CHECK_RESULT(MPI_Allreduce,(local,global,2,MPI_INT64_T,MPI_SUM,*communicator_));
        } else
#endif
            std::copy(local,local+2,global);
        const double energy=-coupling_*global[0],magnetization=global[1];
        return {double(size_),energy,energy*energy,magnetization,magnetization*magnetization,
                magnetization*magnetization*magnetization*magnetization};
    }
    void synchronize() {
        if (fallback_) return;
#ifdef ALPS_HAVE_MPI
        if (counts_.size()>1) {
            BOOST_MPI_CHECK_RESULT(MPI_Gatherv,(spins_.data()+1,counts_[rank_],MPI_INT,
                rank_==0 ? global_.data() : nullptr,counts_.data(),offsets_.data(),MPI_INT,0,*communicator_));
        } else
#endif
            global_.assign(spins_.begin()+1,spins_.end()-1);
        synchronized_=true;
    }
    void save(alps::hdf5::archive& ar) const {
        if (fallback_) {fallback_->save(ar);return;}
        if (rank_!=0 || (counts_.size()>1 && !synchronized_))
            throw std::logic_error("Gather spatial Ising state before saving on its root");
        ar["checkpoint/spins"]<<(counts_.size()==1 ? std::vector<int>(spins_.begin()+1,spins_.end()-1) : global_);
        ar["checkpoint/topology"]<<topology_;
    }
    void load(alps::hdf5::archive& ar) {
        if (fallback_) {fallback_->load(ar);return;}
        std::vector<int> spins; std::vector<uint64_t> graph;
        ar["checkpoint/spins"]>>spins; ar["checkpoint/topology"]>>graph;
        if (graph!=topology_ || spins.size()!=size_ ||
            !std::all_of(spins.begin(),spins.end(),[](int s){return s==1 || s==-1;}))
            throw std::invalid_argument("Invalid spatial Ising checkpoint spins or lattice");
        restore_local(spins);
        if (rank_==0) global_=std::move(spins);
        synchronized_=true;
    }
private:
    struct phase {
        std::vector<size_t> sites;
        std::vector<int> counts,offsets;
    };
    size_t color(size_t i) const {return size_%2 && i==size_-1 ? 2 : i%2;}
    void partition(int ranks) {
        counts_.assign(ranks,int(size_/ranks)); counts_[0]+=size_%ranks;
        offsets_.resize(ranks);
        for (int r=1;r<ranks;++r) offsets_[r]=offsets_[r-1]+counts_[r-1];
        phases_=std::vector<phase>(size_%2 ? 3 : 2);
        for (auto& phase:phases_) {phase.counts.resize(ranks);phase.offsets.resize(ranks);}
        for (int r=0;r<ranks;++r) {
            for (int i=0;i<counts_[r];++i) {
                auto& phase=phases_[color(offsets_[r]+i)]; ++phase.counts[r];
                if (r==rank_) phase.sites.push_back(i);
            }
            if (r>0) for (auto& phase:phases_) phase.offsets[r]=phase.offsets[r-1]+phase.counts[r-1];
        }
        size_t local_max=0,global_max=0;
        for (auto const& phase:phases_) {
            local_max=std::max(local_max,phase.sites.size());
            global_max=std::max(global_max,size_t(phase.offsets.back()+phase.counts.back()));
        }
        proposals_=std::vector<double>(local_max);
        global_proposals_=std::vector<double>(rank_==0 ? global_max : 0);
        spins_=std::vector<int>(counts_[rank_]+2);
        restore_local(global_);
    }
    void restore_local(std::vector<int> const& global) {
        const size_t begin=offsets_[rank_],count=counts_[rank_];
        std::copy(global.begin()+begin,global.begin()+begin+count,spins_.begin()+1);
        spins_.front()=global[(begin+size_-1)%size_]; spins_.back()=global[(begin+count)%size_];
    }
    void halos() {
#ifdef ALPS_HAVE_MPI
        if (counts_.size()>1) {
            const int previous=(rank_+counts_.size()-1)%counts_.size(),next=(rank_+1)%counts_.size();
            communicator_->sendrecv(next,0,spins_[spins_.size()-2],previous,0,spins_.front());
            communicator_->sendrecv(previous,1,spins_[1],next,1,spins_.back());
            return;
        }
#endif
        spins_.front()=spins_[spins_.size()-2]; spins_.back()=spins_[1];
    }
    double coupling_;
    size_t size_;
    int rank_=0;
    bool synchronized_=true;
    std::optional<ising_kernel> fallback_;
    std::vector<uint64_t> topology_;
    std::vector<int> spins_,global_,counts_,offsets_;
    std::vector<double> proposals_,global_proposals_;
    std::vector<phase> phases_;
#ifdef ALPS_HAVE_MPI
    std::optional<boost::mpi::communicator> communicator_;
#endif
};
