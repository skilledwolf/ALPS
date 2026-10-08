// Copyright (C) 1997-2012 Synge Todo; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include "single/native.hpp"
#include "classical_schema.hpp"
#include <alps/mc/replica_exchange.hpp>
#include <alps/mc/replica_parallel.hpp>
#include <algorithm>
#include <array>
#include <optional>

// The same bond Hamiltonian for Ising heat-bath and Heisenberg Metropolis
// updates. Self-bonds contribute constant energy and never a local field.
template<size_t Dimension> class classical_walker {
    static_assert(Dimension==1 || Dimension==3);
    alps::random01 random_;
    std::vector<std::vector<double>> spins_;
    std::vector<std::vector<std::pair<size_t,double>>> neighbors_;
    std::vector<uint64_t> topology_;
    std::vector<double> proposal() {
        if constexpr (Dimension==1) return {random_()<.5 ? -1. : 1.};
        else {
            const double z=2*random_()-1,phi=2*std::acos(-1.)*random_(),radius=std::sqrt(std::max(0.,1-z*z));
            return {radius*std::cos(phi),radius*std::sin(phi),z};
        }
    }
public:
    classical_walker(alps::params const& p,size_t offset):random_(p["SEED"].as<int>()+offset,p["RNG"].as<std::string>()) {
        alps::graph_helper<> graph(p);
        if (!graph.num_sites()) throw std::invalid_argument("The lattice must contain sites");
        neighbors_.resize(graph.num_sites());
        for (auto [it,end]=graph.bonds();it!=end;++it) {
            auto a=graph.source(*it),b=graph.target(*it);
            auto type=graph.bond_type(*it);
            double coupling=p.value_or<double>("J"+std::to_string(type),p["J"].as<double>());
            if (!std::isfinite(coupling)) throw std::invalid_argument("Couplings must be finite");
            neighbors_[a].emplace_back(b,coupling); neighbors_[b].emplace_back(a,coupling);
            topology_.insert(topology_.end(),{uint64_t(a),uint64_t(b),uint64_t(type)});
        }
        for (size_t i=0;i<neighbors_.size();++i) spins_.push_back(proposal());
    }
    static auto names() {
        std::vector<std::string> result{"Number of Sites","Energy","Energy^2",
            Dimension==1 ? "Magnetization" : "Magnetization Z","Magnetization^2","Magnetization^4"};
        if constexpr (Dimension==3) {result.push_back("Magnetization Z^2"); result.push_back("Magnetization Z^4");}
        return result;
    }
    void step(double beta) {
        for (size_t i=0;i<spins_.size();++i) {
            std::array<double,Dimension> field{};
            for (auto [j,coupling]:neighbors_[i]) if (i!=j)
                for (size_t k=0;k<Dimension;++k) field[k]+=coupling*spins_[j][k];
            if constexpr (Dimension==1) {
                if (ising_kernel::heatbath_flip(2*spins_[i][0]*field[0],beta,random_())) spins_[i][0]*=-1;
            } else {
                auto next=proposal(); double change=0;
                for (size_t k=0;k<Dimension;++k) change+=(spins_[i][k]-next[k])*field[k];
                if (std::log(random_())<-beta*change) spins_[i]=std::move(next);
            }
        }
    }
    void synchronize() {}
    double energy() const {
        double energy=0;
        for (size_t i=0;i<spins_.size();++i) for (auto [j,coupling]:neighbors_[i])
            for (size_t k=0;k<Dimension;++k) energy-=.5*coupling*spins_[i][k]*spins_[j][k];
        return energy;
    }
    auto sample() const {
        const double e=energy(),n=spins_.size();
        std::array<double,Dimension> magnetization{};
        for (auto const& spin:spins_) for (size_t k=0;k<Dimension;++k) magnetization[k]+=spin[k]/n;
        double m2=0;for (auto value:magnetization) m2+=value*value;
        const double z=magnetization.back();
        std::vector<double> result{n,e,e*e,z,m2,m2*m2};
        if constexpr (Dimension==3) {result.push_back(z*z);result.push_back(z*z*z*z);}
        return result;
    }
    void save(alps::hdf5::archive& ar) const {
        ar["spins"]<<spins_; ar["topology"]<<topology_; ar["rng"]<<random_;
    }
    void load(alps::hdf5::archive& ar) {
        std::vector<uint64_t> topology;
        ar["topology"]>>topology;
        auto spins=spins_;auto random=random_;
        ar["spins"]>>spins; ar["rng"]>>random;
        if (topology!=topology_ || spins.size()!=spins_.size() || random.name()!=random_.name())
            throw std::invalid_argument("Checkpoint graph or RNG changed");
        for (auto const& spin:spins) {
            double norm=std::inner_product(spin.begin(),spin.end(),spin.begin(),0.);
            if (spin.size()!=Dimension || !std::isfinite(norm) || std::abs(norm-1)>1.e-12)
                throw std::invalid_argument("Invalid classical checkpoint spin");
            if constexpr (Dimension==1) if (std::abs(spin[0])!=1) throw std::invalid_argument("Invalid Ising spin");
        }
        spins_=std::move(spins);random_=std::move(random);
    }
};

// Groups that spread one walker over several ranks check owned walkers
// locally on construction, bind them to their team before collective work,
// and store each walker from one rank.
template<class Group,class Walkers> auto check_walkers(Group const& group,Walkers const& walkers,int)
    -> decltype(group.check(walkers),void()) {group.check(walkers);}
template<class Group,class Walkers> void check_walkers(Group const&,Walkers const&,long) {}
template<class Group,class Walkers> auto bind_walkers(Group& group,Walkers& walkers,int)
    -> decltype(group.bind(walkers),void()) {group.bind(walkers);}
template<class Group,class Walkers> void bind_walkers(Group&,Walkers&,long) {}
template<class Group> auto stores_walker(Group const& group,size_t i,int)
    -> decltype(group.stores_walker(i)) {return group.stores_walker(i);}
template<class Group> bool stores_walker(Group const& group,size_t i,long) {return group.owns_walker(i);}

template<class Walker,class Group=alps::mc::replica_parallel> class classical_application {
    alps::params parameters_;
    size_t bins_,chain_;
    uint64_t steps_=0;
    std::vector<std::optional<Walker>> walkers_;
    std::unique_ptr<alps::mc::replica_checkpoint> checkpoint_;
    Group group_;
    std::vector<spin_statistics> stats_;
    std::vector<std::map<std::string,alps::mc::batch>> diagnostics_;
    std::optional<alps::mc::replica_exchange<double>> exchange_;
    uint64_t production() const {
        auto warm=parameters_["THERMALIZATION"].template as<uint64_t>();
        return exchange_ ? exchange_->production_sweeps() : steps_>warm ? steps_-warm : 0;
    }
    uint64_t diagnostic_count(size_t i,std::string const& name) const {
        if (name=="EXMC: Temperature" || name=="EXMC: Inverse Temperature") return production();
        auto interval=parameters_["EXCHANGE_INTERVAL"].template as<uint64_t>();
        auto first=(completed_sweeps()-production())/interval,last=completed_sweeps()/interval;
        if (name!="EXMC: Acceptance Rate" || parameters_["RANDOM_EXCHANGE"].template as<bool>()) return last-first;
        auto even=last/2-first/2;
        return i%2 ? last-first-even : even;
    }
public:
    classical_application(alps::params const& p,size_t bins,size_t chain,Group group):parameters_(p),bins_(bins),chain_(chain),group_(std::move(group)) {
        if (p["ALGORITHM"].as<std::string>().find("exchange")!=std::string::npos) exchange_.emplace(p,chain,0.,true);
        if (group_.distributed && !exchange_) throw std::invalid_argument("parallel = replicas requires an exchange algorithm");
        const size_t n=exchange_ ? exchange_->size() : 1;
        walkers_.resize(n);
        for (size_t i=0;i<n;++i) {
            if (group_.owns_walker(i)) walkers_[i].emplace(p,exchange_ ? chain*(n+1)+i : chain);
            stats_.emplace_back(bins,names().size());
        }
        // Idle ranks still validate the model before collective preflight.
        if (std::none_of(walkers_.begin(),walkers_.end(),[](auto const& w){return bool(w);})) Walker(p,chain*(n+1));
        check_walkers(group_,walkers_,0);
        diagnostics_.resize(n);
        if (exchange_) exchange_->init_diagnostics([&](size_t i,char const* name){diagnostics_[i].emplace(name,alps::mc::batch(1,bins));});
    }
    static auto names() {return Walker::names();}
    static alps::params checkpoint_parameters(alps::params p) {p.erase("SWEEPS");return p;}
    size_t stages() const {return walkers_.size();}
    double inverse_temperature(size_t i) const {return exchange_ ? exchange_->beta(i) : parameters_["BETA"].template as<double>();}
    double temperature(size_t i) const {return 1/inverse_temperature(i);}
    auto const& stage_statistics(size_t i) const {return stats_.at(i);}
    auto const& diagnostics(size_t i) const {return diagnostics_.at(i);}
    uint64_t completed_sweeps() const {return exchange_ ? exchange_->completed_sweeps() : steps_;}
    double fraction_completed() const {return double(production())/parameters_["SWEEPS"].template as<uint64_t>();}
    void synchronize() {
        bind_walkers(group_,walkers_,0);
        for (auto& walker:walkers_) if (walker) walker->synchronize();
        group_.synchronize_physical([&](auto& ar) {
            for (size_t i=0;i<stages();++i) if (walkers_[i] && stores_walker(group_,i,0))
                ar["/walkers/"+std::to_string(i)]<<*walkers_[i];
        },checkpoint_,"/walkers");
    }
    void update() {
        bind_walkers(group_,walkers_,0);
        if (!exchange_) {
            walkers_[0]->step(parameters_["BETA"].template as<double>());++steps_;
            if (production()) stats_[0].add(walkers_[0]->sample());
            return;
        }
        bool sampling=false;
        std::vector<double> energies(stages());
        exchange_->step([&](auto const& walkers,auto const& betas,bool record) {
            sampling=record;
            const auto components=names().size();
            std::vector<double> samples(stages()*components);
            group_.update([&] {
                for (size_t i=0;i<stages();++i) if (group_.owns_walker(walkers[i])) {
                    auto& walker=*walkers_[walkers[i]];walker.step(betas[i]);
                    auto sample=walker.sample();
                    std::copy(sample.begin(),sample.end(),samples.begin()+i*components);
                }
            });
            group_.collect(samples);
            for (size_t i=0;i<stages();++i) {
                energies[walkers[i]]=samples[i*components+1];
                if (record) stats_[i].add(std::vector<double>(samples.begin()+i*components,samples.begin()+(i+1)*components));
            }
        },[&] {return energies;
        },[](double energy,double beta){return -beta*energy;},
        [&](size_t i,char const* name,double value) {if (sampling) diagnostics_[i].at(name)<<alps::alea::column<double>{value};});
    }
    void measure() {}
    void save(alps::hdf5::archive& ar) const {
        ar["/parameters"]<<parameters_;ar["chain"]<<uint64_t(chain_);
        if (exchange_) exchange_->save(ar); else ar["steps"]<<steps_;
        for (size_t i=0;i<stages();++i) {
            auto path="walkers/"+std::to_string(i);
            if (group_.distributed) checkpoint_->copy("/"+path,ar,path);
            else ar[path]<<*walkers_[i];
            stats_[i].save(ar,"stages/"+std::to_string(i));
            alps::alea::hdf5_serializer codec(ar,"stages/"+std::to_string(i)+"/exchange");
            for (auto const& [name,value]:diagnostics_[i]) serialize(codec,ar.encode_segment(name),value);
        }
    }
    void load(alps::hdf5::archive& ar) {
        classical_application restored(parameters_,bins_,chain_,group_);alps::params p;uint64_t chain;
        ar["/parameters"]>>p;ar["chain"]>>chain;
        if (checkpoint_parameters(p)!=checkpoint_parameters(parameters_) || chain!=chain_ ||
                ar.list_children("walkers").size()!=stages() || ar.list_children("stages").size()!=stages())
            throw std::invalid_argument("Checkpoint does not match the classical model or ladder");
        if (exchange_) restored.exchange_->load(ar);
        else {
            ar["steps"]>>restored.steps_;
            if (restored.production()>parameters_["SWEEPS"].template as<uint64_t>()) throw std::invalid_argument("Invalid classical checkpoint counter");
        }
        for (size_t i=0;i<stages();++i) {
            if (restored.walkers_[i]) ar["walkers/"+std::to_string(i)]>>*restored.walkers_[i];
            restored.stats_[i].load(ar,"stages/"+std::to_string(i),restored.production());
            if (!exchange_) continue;
            const auto path="stages/"+std::to_string(i)+"/exchange";
            if (ar.list_children(path).size()!=diagnostics_[i].size()) throw std::invalid_argument("Invalid exchange diagnostic names");
            alps::alea::hdf5_serializer codec(ar,path);
            for (auto& [name,value]:restored.diagnostics_[i]) {
                deserialize(codec,ar.encode_segment(name),value);
                if (value.size()!=1 || value.num_batches()!=bins_ || value.count()!=restored.diagnostic_count(i,name)
                        || value.current_batch_size()!=value.cursor().factor() || !value.store().batch().allFinite())
                    throw std::invalid_argument("Invalid exchange diagnostic measurements");
            }
        }
        *this=std::move(restored);
    }
};

// The first algorithm is the default; names containing "exchange" run a ladder.
template<class Walker,class Group=alps::mc::replica_parallel>
int classical_main(int argc,char** argv,char const* command,std::vector<std::string> algorithms,std::string const& group_schema="") {
    using simulation=classical_application<Walker,Group>;
    auto schema=std::string(classical_schema)+spin_common_schema+alps::mc::replica_exchange_schema+alps::mc::replica_parallel_schema+group_schema;
    schema+="[parameters.ALGORITHM]\ntype = \"string\"\ndefault = \""+algorithms.front()+"\"\nchoices = [";
    for (size_t i=0;i<algorithms.size();++i) schema+=(i ? ", \"" : "\"")+algorithms[i]+"\"";
    schema+="]\n";
    return alps::mc::main<simulation,Group>(argc,argv,command,schema.c_str(),
        [](std::string const& key)->char const* {
            return key.size()>1 && key[0]=='J' && key.find_first_not_of("0123456789",1)==std::string::npos ? "float64" : nullptr;
        },[=](alps::params& p,alps::run_configuration const& run) {
            auto algorithm=p.value_or<std::string>("ALGORITHM",algorithms.front());p["ALGORITHM"]=algorithm;
            const bool exchange=algorithm.find("exchange")!=std::string::npos;
            if (!p.exists("THERMALIZATION")) p["THERMALIZATION"]=p["SWEEPS"].as<uint64_t>()/8;
            if (exchange) {
                if (p.exists("T") || p.exists("BETA")) throw std::invalid_argument("Exchange uses a temperature ladder, not T/BETA");
                if ((p["OPTIMIZE_TEMPERATURE"].as<bool>() || p["TEMPERATURE_OPTIMIZATION"].as<bool>()) && run.execution["chains"].as<size_t>()!=1)
                    throw std::invalid_argument("Temperature optimization requires execution.chains = 1");
            } else {
                for (auto key:{"NUM_REPLICAS","TEMPERATURE_SET","INVERSE_TEMPERATURE_SET","T_MIN","T_MAX","BETA_MIN","BETA_MAX","INITIAL_BLOCK_SWEEPS","OPTIMIZATION_ITERATIONS","BLOCK_SWEEP_FACTOR"})
                    if (p.exists(key)) throw std::invalid_argument("Replica ladder parameters require the exchange algorithm");
                if (p["OPTIMIZE_TEMPERATURE"].as<bool>() || p["TEMPERATURE_OPTIMIZATION"].as<bool>())
                    throw std::invalid_argument("Temperature optimization requires the exchange algorithm");
                if (p.exists("T") && p.exists("BETA")) throw std::invalid_argument("Specify only T or BETA");
                double beta=p.exists("T") ? 1/p["T"].as<double>() : p.value_or<double>("BETA",0.);
                if (!std::isfinite(beta) || beta<0) throw std::invalid_argument("BETA must be finite and nonnegative; T must be positive");
                p["BETA"]=beta;
                if (beta>0 && !p.exists("T")) p["T"]=1/beta;
            }
        },[](auto const& run,auto const& chains,alps::params const& p) {
            alps::hdf5::save_checkpoint(run.output["results"].template as<std::string>(),[&](auto& ar) {
                write_spin_moments(ar,run,chains,p,simulation::names());
                for (size_t i=0;i<chains.front()->stages();++i) {
                    auto path="/simulation/replicas/"+std::to_string(i)+"/results";
                    alps::alea::hdf5_serializer codec(ar,path);
                    for (auto const& [name,unused]:chains.front()->diagnostics(i)) {
                        std::vector<alps::alea::batch_result<double>> results;
                        for (auto const& chain:chains) results.push_back(chain->diagnostics(i).at(name).result());
                        serialize(codec,ar.encode_segment(name),alps::alea::merge(results));
                    }
                }
            });
        });
}
