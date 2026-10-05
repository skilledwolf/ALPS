// Copyright (C) 1999–2009 Matthias Troyer, Fabian Stoeckli;
// modifications (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "../moment_difference.hpp"
#include <alps/mcbase.hpp>
#include <alps/lattice.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/alea/checkpoint.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/alea/transform.hpp>
#include <alps/hdf5/stdarray.hpp>
#include <alps/hdf5/vector.hpp>
#include <cmath>
#include <iostream>
#include <limits>

namespace spinmc {
using results_type = std::map<std::string,alps::alea::batch_result<double>>;
inline std::vector<double> values(alps::params const& p, std::string const& key, std::vector<double> fallback) {
    if (!p.exists(key)) return fallback;
    if (p[key].isType<std::vector<double>>()) return p[key].as<std::vector<double>>();
    return {p[key].as<double>()};
}
inline double inverse_temperature(alps::params const& p) {
    if (p.exists("T") == p.exists("beta")) throw std::invalid_argument("Specify exactly one spinmc T or beta");
    double beta = p.exists("beta") ? p["beta"].as<double>() : 1. / p["T"].as<double>();
    if (!std::isfinite(beta) || beta < 0 || (p.exists("T") && !(p["T"].as<double>() > 0)))
        throw std::invalid_argument("spinmc requires a finite nonnegative inverse temperature");
    return beta;
}
inline Eigen::MatrixXd matrix(std::vector<double> const& values, std::size_t n) {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n,n);
    if (values.size() == 1) result.diagonal().setConstant(values[0]);
    else if (values.size() == n)
        for (std::size_t i=0; i<n; ++i) result(i,i)=values[i];
    else if (values.size() == n*(n+1)/2) {
        std::size_t k=0;
        for (std::size_t i=0; i<n; ++i)
            for (std::size_t j=0; j<=i; ++j) result(i,j)=result(j,i)=values[k++];
    } else if (values.size() == n*n)
        for (std::size_t i=0; i<n; ++i)
            for (std::size_t j=0; j<n; ++j) result(i,j)=values[i*n+j];
    else throw std::invalid_argument("spinmc matrix needs scalar, diagonal, packed symmetric or full values");
    if (!result.allFinite()) throw std::invalid_argument("nonfinite spinmc matrix");
    return result;
}
inline bool isotropic(Eigen::MatrixXd const& m) {
    return m == m(0,0) * Eigen::MatrixXd::Identity(m.rows(),m.cols());
}

// Classical H = -sum_ij S_i S_j s_i^T J_ij s_j
//               -sum_i S_i^2 s_i^T D_i s_i - g sum_i S_i h.s_i.
// Quantum convention reverses J and uses sqrt(S(S+1)) for bond/field moments.
// Potts replaces the bond dot product with equal-color indicators and uses
// color0 occupation for its field and magnetization measurements.
class simulation : public alps::mcbase, private alps::graph_helper<> {
    struct bond { std::size_t source, target, type; Eigen::MatrixXd coupling; };
public:
    using results_type = std::map<std::string,alps::alea::batch_result<double>>;
    results_type collect_results(result_names_type const& names={}) const {
        return collect_results_as<alps::alea::batch_result<double>>(names);
    }
    simulation(alps::params const& p, std::size_t bins=128, std::size_t chain=0)
        : mcbase(p,chain), graph_helper<>(graph_parameters(p)), bins_(bins), chain_(chain),
          model_(p["MODEL"].as<std::string>()), dim_(model_dimension(model_)), potts_(model_=="Potts"),
          q_(potts_ ? p["q"].as<unsigned>() : 0), beta_(inverse_temperature(p)),
          production_(p["SWEEPS"].as<uint64_t>()),
          thermalization_(p.value_or<uint64_t>("THERMALIZATION", production_/10)),
          neighbors_(num_sites()), factors_(num_sites()), onsite_(num_sites()) {
        if (!num_sites() || !production_ || (potts_ && q_!=3 && q_!=4 && q_!=10))
            throw std::invalid_argument("spinmc requires sites, positive SWEEPS and Potts q=3,4 or10");
        if (thermalization_ > (UINT64_MAX-production_)/num_sites()
                || thermalization_ * num_sites() > UINT64_MAX-(num_sites()-1))
            throw std::invalid_argument("spinmc thermalization counters would overflow");
        if (inhomogeneous()) throw std::invalid_argument("spinmc does not support inhomogeneous lattice types");
        auto convention=p.value_or("CONVENTION","classical");
        if (convention!="classical" && convention!="quantum") throw std::invalid_argument("unknown spinmc convention");
        bool quantum=convention=="quantum", cluster_legal=true;
        double g=p.value_or("g",1.);
        auto h=values(p,"h",{0.});
        field_=Eigen::VectorXd::Zero(dim_);
        if (h.size()==1) field_(dim_-1)=h[0];
        else if (h.size()==dim_) for (std::size_t i=0;i<dim_;++i) field_(i)=h[i];
        else throw std::invalid_argument("spinmc h needs one value or one per component");
        if (potts_ && field_(0)!=0)
            throw std::invalid_argument("Potts h requires a scalar or [0,h]; its first component is unused");
        if (!field_.allFinite() || !std::isfinite(g) || (quantum && !field_.isZero(0)))
            throw std::invalid_argument("invalid spinmc field or quantum convention with nonzero field");
        auto field_norm=field_.stableNorm();
        direction_=field_norm ? Eigen::VectorXd(field_/field_norm) : Eigen::VectorXd::Zero(dim_);
        field_*=g;
        if (!field_.allFinite() || !std::isfinite(field_norm)) throw std::invalid_argument("spinmc field would overflow");
        cluster_legal=field_.isZero(0);
        double energy_bound=0, magnetization_bound=0;
        for (std::size_t site=0;site<num_sites();++site) {
            auto type=std::to_string(site_type(site));
            double s=p.value_or("S"+type,p.value_or("S",quantum?.5:1.));
            if (!std::isfinite(s) || s<0) throw std::invalid_argument("spinmc S must be finite and nonnegative");
            factors_[site]=std::sqrt(s*(s+(quantum?1.:0.)));
            onsite_[site]=matrix(values(p,"D"+type,values(p,"D",{0.})),dim_)*s*s;
            if (!std::isfinite(factors_[site]) || !onsite_[site].allFinite()) throw std::invalid_argument("spinmc spin or onsite interaction would overflow");
            if (potts_ && !onsite_[site].isZero(0)) throw std::invalid_argument("Potts does not define a D interaction");
            if (onsite_[site]!=onsite_[site].transpose()) throw std::invalid_argument("spinmc D must be symmetric");
            cluster_legal=cluster_legal && isotropic(onsite_[site]);
            energy_bound+=onsite_[site].cwiseAbs().sum()+factors_[site]*field_.stableNorm();
            magnetization_bound+=factors_[site]/num_sites();
        }
        for (auto [it,end]=bonds();it!=end;++it) {
            auto u=std::size_t(source(*it)), v=std::size_t(target(*it)), type=std::size_t(bond_type(*it));
            auto coupling=matrix(values(p,"J"+std::to_string(type),values(p,"J",{1.})),dim_);
            if (potts_ && !isotropic(coupling)) throw std::invalid_argument("Potts requires isotropic scalar J");
            coupling*= (quantum?-1.:1.)*factors_[u]*factors_[v];
            if (!coupling.allFinite()) throw std::invalid_argument("spinmc bond interaction would overflow");
            cluster_legal=cluster_legal && isotropic(coupling);
            if (u!=v) {
                ferro_=ferro_ && coupling(0,0)>=0;
                antiferro_=antiferro_ && coupling(0,0)<=0;
                if (potts_ && coupling(0,0)<0) cluster_legal=false;
            }
            energy_bound+=coupling.cwiseAbs().sum();
            neighbors_[u].push_back(edges_.size());
            if (u!=v) neighbors_[v].push_back(edges_.size());
            edges_.push_back({u,v,type,std::move(coupling)});
            bond_types_=std::max(bond_types_,type+1);
        }
        auto update=p.value_or("UPDATE","auto");
        if (update!="local" && update!="cluster" && update!="auto") throw std::invalid_argument("unknown spinmc update");
        if (cluster_legal) cluster_legal=unfrustrated();
        if (update=="cluster" && !cluster_legal)
            throw std::invalid_argument("spinmc cluster updates require zero field, isotropic matrices and unfrustrated couplings; Potts requires ferromagnetic couplings");
        cluster_=update!="local" && cluster_legal;
        double largest=std::max({1.,double(num_sites()),energy_bound*energy_bound,
            energy_bound*std::pow(magnetization_bound,4),std::pow(magnetization_bound,4),
            beta_*num_sites()*magnetization_bound*magnetization_bound});
        if (!std::isfinite(energy_bound) || !std::isfinite(magnetization_bound)
                || !(largest<=std::sqrt(std::numeric_limits<double>::max())/(4*double(production_))))
            throw std::invalid_argument("spinmc measurement moments would overflow");
        spins_=Eigen::MatrixXd::Constant(potts_?1:dim_,num_sites(),potts_?0.:1./std::sqrt(double(dim_)));
        for (auto const& name : {"Number of Sites","Energy","Energy Density","Energy^2","|Magnetization|",
             "Magnetization along Field","Magnetization^2","Magnetization^4","E.Magnetization^2","E.Magnetization^4","Susceptibility"}) add(name);
        if (bond_types_) { add("Bond-type Energy",bond_types_); add("Bond-type Energy Density",bond_types_); }
        if (is_bipartite()) { add("|Staggered Magnetization|"); add("Staggered Magnetization^2"); }
        if (cluster_) {
            add("Cluster size");
            if (!potts_ && ferro_) { add("Improved Magnetization^2"); add("Improved Susceptibility"); }
            if (!potts_ && antiferro_ && is_bipartite()) { add("Improved Staggered Magnetization^2"); add("Improved Staggered Susceptibility"); }
        }
        if (p.exists("ERROR_VARIABLE") != p.exists("ERROR_LIMIT")) throw std::invalid_argument("spinmc error stopping requires variable and limit");
        if (p.exists("ERROR_VARIABLE")) {
            error_variable_=p["ERROR_VARIABLE"].as<std::string>(); error_limit_=p["ERROR_LIMIT"].as<double>();
            if (!(error_limit_>0) || !std::isfinite(error_limit_) || !measurements.count(error_variable_)
                    || measurement(error_variable_)->size()!=1) throw std::invalid_argument("invalid spinmc scalar error stopping criterion");
        }
    }
    simulation(simulation const&)=delete;
    simulation& operator=(simulation const&)=delete;
    using mcbase::save;
    using mcbase::load;
    std::string effective_update() const { return cluster_?"cluster":"local"; }
    uint64_t completed_sweeps() const { return updates_; }
    uint64_t measurement_count() const { return measurement("Energy")->count(); }
    std::size_t chain_id() const { return chain_; }
    static alps::params checkpoint_parameters(alps::params const& p) {
        auto identity=p;
        for (auto const* key:{"ERROR_VARIABLE","ERROR_LIMIT","PRINT_SWEEPS"}) identity.erase(key);
        return identity;
    }
    double fraction_completed() const override {
        double progress=double(measurement_count())/production_;
        if (!error_variable_.empty()) {
            auto result=measurement(error_variable_)->result();
            if ((result.store().count().array()>0).count()>1) {
                double error=result.stderror()(0);
                if (std::isfinite(error) && error<=error_limit_) progress=1.;
            }
        }
        return progress;
    }
    void update() override {
        if (fraction_completed()>=1.) throw std::logic_error("spinmc run is complete");
        if (measurement_count()!=updates_-warmup_updates_) throw std::logic_error("measure spinmc before updating again");
        cluster_size_=num_sites(); cluster_projection_=staggered_projection_=0;
        if (cluster_) cluster_update();
        else for (std::size_t i=0;i<num_sites();++i) local_update(std::size_t(random()*num_sites()));
        ++updates_;
        if (warmup_sites_<thermalization_*num_sites()) { warmup_sites_+=cluster_size_; ++warmup_updates_; }
    }
    void measure() override {
        if (updates_==warmup_updates_) return;
        if (measurement_count()+1!=updates_-warmup_updates_) throw std::logic_error("measure each spinmc update exactly once");
        double energy=0;
        Eigen::VectorXd bond_energy=Eigen::VectorXd::Zero(bond_types_);
        for (auto const& edge:edges_) {
            double e=edge_energy(edge,spins_.col(edge.source),spins_.col(edge.target));
            energy+=e; bond_energy(edge.type)+=e;
        }
        Eigen::VectorXd magnetization=Eigen::VectorXd::Zero(potts_?1:dim_), staggered=magnetization;
        double field_projection=0;
        for (std::size_t site=0;site<num_sites();++site) {
            energy+=site_energy(site,spins_.col(site));
            Eigen::VectorXd moment=potts_?Eigen::VectorXd::Constant(1,spins_(0,site)==0):Eigen::VectorXd(spins_.col(site));
            moment*=factors_[site]/num_sites();
            magnetization+=moment; staggered+=parity(site)*moment;
            field_projection+=potts_ ? moment(0)*direction_(dim_-1) : moment.dot(direction_);
        }
        double m2=magnetization.squaredNorm(), m4=m2*m2;
        record("Number of Sites",double(num_sites())); record("Energy",energy); record("Energy Density",energy/num_sites()); record("Energy^2",energy*energy);
        if (bond_types_) { record("Bond-type Energy",bond_energy); record("Bond-type Energy Density",Eigen::VectorXd(bond_energy/num_sites())); }
        record("|Magnetization|",magnetization.norm()); record("Magnetization along Field",field_projection);
        record("Magnetization^2",m2); record("Magnetization^4",m4);
        record("E.Magnetization^2",energy*m2); record("E.Magnetization^4",energy*m4);
        record("Susceptibility",beta_*num_sites()*m2/(potts_?1:dim_));
        if (is_bipartite()) { record("|Staggered Magnetization|",staggered.norm()); record("Staggered Magnetization^2",staggered.squaredNorm()); }
        if (cluster_) {
            record("Cluster size",double(cluster_size_)/num_sites());
            auto improved=[&](std::string const& prefix,double projected) {
                record(prefix+"Magnetization^2",dim_*projected*projected/(double(cluster_size_)*num_sites()));
                record(prefix+"Susceptibility",beta_*projected*projected/cluster_size_);
            };
            if (!potts_ && ferro_) improved("Improved ",cluster_projection_);
            if (!potts_ && antiferro_ && is_bipartite()) improved("Improved Staggered ",staggered_projection_);
        }
        auto print=parameters.value_or<uint64_t>("PRINT_SWEEPS",0);
        if (print && measurement_count()%print==0)
            for (std::size_t site=0;site<num_sites();++site) {
                std::cout<<site<<' '; for (double c:coordinate(site)) std::cout<<c<<' ';
                std::cout<<spins_.col(site).transpose()<<'\n';
            }
    }
    void save(alps::hdf5::archive& ar) const override {
        if (measurement_count()!=updates_-warmup_updates_) throw std::logic_error("measure spinmc before checkpointing");
        mcbase::save(ar);
        ar["checkpoint/chain_id"]<<uint64_t(chain_); ar["checkpoint/updates"]<<updates_;
        ar["checkpoint/warmup_updates"]<<warmup_updates_; ar["checkpoint/warmup_sites"]<<warmup_sites_;
        ar["checkpoint/topology"]<<topology(); ar["checkpoint/hamiltonian"]<<hamiltonian();
        ar["checkpoint/update"]<<effective_update();
        alps::alea::hdf5_serializer codec(ar,"checkpoint");
        alps::serialization::serialize(codec,"spins",spins_);
    }
    void load(alps::hdf5::archive& ar) override {
        alps::params p; uint64_t chain,updates,warmup_updates,warmup_sites;
        std::string update; Eigen::MatrixXd spins(spins_.rows(),spins_.cols());
        std::vector<std::array<uint64_t,3>> graph; std::vector<double> coefficients;
        ar["/parameters"]>>p; ar["checkpoint/chain_id"]>>chain; ar["checkpoint/updates"]>>updates;
        ar["checkpoint/warmup_updates"]>>warmup_updates; ar["checkpoint/warmup_sites"]>>warmup_sites;
        ar["checkpoint/topology"]>>graph; ar["checkpoint/hamiltonian"]>>coefficients; ar["checkpoint/update"]>>update;
        alps::alea::hdf5_serializer codec(ar,"checkpoint");
        alps::serialization::deserialize(codec,"spins",spins);
        uint64_t threshold=thermalization_*num_sites();
        if (checkpoint_parameters(p)!=checkpoint_parameters(parameters) || chain!=chain_ || graph!=topology() || coefficients!=hamiltonian() || update!=effective_update()
                || warmup_updates>updates || updates-warmup_updates>production_ || warmup_sites<warmup_updates
                || warmup_sites/num_sites()+bool(warmup_sites%num_sites())>warmup_updates
                || warmup_sites>threshold+num_sites()-1 || (!threshold && (warmup_updates || warmup_sites))
                || (warmup_sites<threshold && updates!=warmup_updates)
                || (!cluster_ && (warmup_updates!=std::min(updates,thermalization_) || warmup_sites!=warmup_updates*num_sites()))
                || spins.rows()!=spins_.rows() || spins.cols()!=spins_.cols() || !spins.allFinite())
            throw std::invalid_argument("spinmc checkpoint does not match this run");
        for (Eigen::Index site=0;site<spins.cols();++site) {
            bool valid=potts_ ? spins(0,site)>=0 && spins(0,site)<q_ && spins(0,site)==std::floor(spins(0,site))
                             : std::abs(spins.col(site).squaredNorm()-1.)<=1e-10;
            if (!valid || (!potts_ && dim_==1 && std::abs(spins(0,site))!=1.)) throw std::invalid_argument("invalid spinmc checkpoint state");
        }
        if (ar.list_children("measurements").size()!=measurements.size()) throw std::invalid_argument("unexpected spinmc checkpoint observables");
        alps::alea::hdf5_serializer measurements_codec(ar,"measurements");
        for (auto const& entry:measurements) {
            alps::alea::batch_acc<double> value;
            alps::alea::deserialize(measurements_codec,ar.encode_segment(entry.first),value);
            if (value.size()!=measurement(entry.first)->size() || value.num_batches()!=bins_
                    || value.current_batch_size()!=value.cursor().factor() || value.count()!=updates-warmup_updates
                    || !value.store().batch().allFinite()) throw std::invalid_argument("invalid spinmc checkpoint measurement state");
        }
        auto current_parameters=parameters;
        mcbase::load(ar); parameters=std::move(current_parameters);
        spins_=std::move(spins); updates_=updates; warmup_updates_=warmup_updates; warmup_sites_=warmup_sites;
    }
private:
    static std::size_t model_dimension(std::string const& model) {
        if (model=="Ising") return 1; if (model=="XY" || model=="Potts") return 2;
        if (model=="Heisenberg") return 3; if (model=="O(4)") return 4;
        throw std::invalid_argument("unknown spinmc model");
    }
    static alps::Parameters graph_parameters(alps::params const& p) {
        alps::Disorder::seed(p.value_or<uint32_t>("DISORDER_SEED",uint32_t(p.value_or("SEED",42))));
        return alps::make_deprecated_parameters(p);
    }
    Eigen::VectorXd axis() {
        Eigen::VectorXd direction(dim_);
        if (dim_==1) { direction(0)=1.; return direction; }
        if (dim_==2) { double phi=6.2831853071795864769*random(); direction<<std::cos(phi),std::sin(phi); return direction; }
        if (dim_==3) {
            double phi=6.2831853071795864769*random(), z=2*random()-1, r=std::sqrt(std::max(0.,1-z*z));
            direction<<r*std::cos(phi),r*std::sin(phi),z; return direction;
        }
        // Only O(4) reaches this paired Box–Muller draw; no cached RNG state.
        do {
            for (std::size_t i=0;i<dim_;i+=2) {
                double radius=std::sqrt(-2*std::log(1-random())), phi=6.2831853071795864769*random();
                direction(i)=radius*std::cos(phi); direction(i+1)=radius*std::sin(phi);
            }
        } while (direction.squaredNorm()==0);
        return direction.normalized();
    }
    Eigen::VectorXd proposed(std::size_t site,Eigen::VectorXd const& direction) {
        if (potts_) {
            unsigned color=unsigned((q_-1)*random()); if (color>=spins_(0,site)) ++color;
            return Eigen::VectorXd::Constant(1,color);
        }
        return Eigen::VectorXd(spins_.col(site)-2*direction.dot(spins_.col(site))*direction).normalized();
    }
    double site_energy(std::size_t site,Eigen::VectorXd const& s) const {
        return potts_ ? -(s(0)==0)*field_(dim_-1)*factors_[site]
                      : -s.dot(onsite_[site]*s)-factors_[site]*field_.dot(s);
    }
    double edge_energy(bond const& edge,Eigen::VectorXd const& u,Eigen::VectorXd const& v) const {
        return potts_ ? -edge.coupling(0,0)*(u(0)==v(0)) : -u.dot(edge.coupling*v);
    }
    void local_update(std::size_t site) {
        auto candidate=proposed(site,potts_?Eigen::VectorXd():axis());
        double delta=site_energy(site,candidate)-site_energy(site,spins_.col(site));
        for (auto index:neighbors_[site]) {
            auto const& edge=edges_[index];
            delta+=edge_energy(edge,edge.source==site?candidate:Eigen::VectorXd(spins_.col(edge.source)),
                                   edge.target==site?candidate:Eigen::VectorXd(spins_.col(edge.target)))
                  -edge_energy(edge,spins_.col(edge.source),spins_.col(edge.target));
        }
        // Heat-bath Ising flips avoid the even-sweep parity trap at beta=0.
        bool accept=!potts_ && dim_==1 ? random()<.5*(1+std::tanh(-.5*beta_*delta))
                                      : delta<=0 || random()<std::exp(-beta_*delta);
        if (accept) spins_.col(site)=candidate;
    }
    void cluster_update() {
        std::size_t seed=std::size_t(random()*num_sites());
        auto direction=potts_?Eigen::VectorXd():axis();
        auto replacement=proposed(seed,direction);
        std::vector<std::size_t> stack{seed}; std::vector<bool> member(num_sites(),false); member[seed]=true;
        cluster_size_=0;
        while (!stack.empty()) {
            auto site=stack.back(); stack.pop_back(); ++cluster_size_;
            double projected=potts_?0.:spins_.col(site).dot(direction)*factors_[site];
            cluster_projection_+=projected; staggered_projection_+=parity(site)*projected;
            for (auto index:neighbors_[site]) {
                auto const& edge=edges_[index]; auto neighbor=edge.source==site?edge.target:edge.source;
                if (member[neighbor]) continue;
                double delta=potts_ ? edge.coupling(0,0)*(spins_(0,site)==spins_(0,neighbor))
                    : 2*edge.coupling(0,0)*spins_.col(site).dot(direction)*spins_.col(neighbor).dot(direction);
                if (delta>0 && random() < -std::expm1(-beta_*delta)) { member[neighbor]=true; stack.push_back(neighbor); }
            }
            spins_.col(site)=potts_?replacement:proposed(site,direction);
        }
    }
    bool unfrustrated() const {
        std::vector<int> gauge(num_sites(),0);
        for (std::size_t root=0;root<num_sites();++root) if (!gauge[root]) {
            gauge[root]=1; std::vector<std::size_t> stack{root};
            while (!stack.empty()) {
                auto site=stack.back(); stack.pop_back();
                for (auto index:neighbors_[site]) {
                    auto const& edge=edges_[index]; auto neighbor=edge.source==site?edge.target:edge.source;
                    if (neighbor==site || edge.coupling(0,0)==0) continue;
                    int sign=edge.coupling(0,0)>0?1:-1;
                    if (!gauge[neighbor]) { gauge[neighbor]=gauge[site]*sign; stack.push_back(neighbor); }
                    else if (gauge[neighbor]!=gauge[site]*sign) return false;
                }
            }
        }
        return true;
    }
    void add(std::string const& name,std::size_t components=1) {
        measurements.emplace(name,std::make_shared<alps::alea::batch_acc<double>>(components,bins_));
    }
    template<class T> void record(std::string const& name,T const& value) { *measurement(name)<<alps::alea::make_adapter(value); }
    std::vector<std::array<uint64_t,3>> topology() const {
        std::vector<std::array<uint64_t,3>> result;
        for (auto const& edge:edges_) result.push_back({edge.source,edge.target,edge.type});
        return result;
    }
    std::vector<double> hamiltonian() const {
        std::vector<double> result=factors_;
        for (auto const& m:onsite_) result.insert(result.end(),m.data(),m.data()+m.size());
        for (auto const& edge:edges_) result.insert(result.end(),edge.coupling.data(),edge.coupling.data()+edge.coupling.size());
        return result;
    }
    std::size_t bins_,chain_;
    std::string model_,error_variable_;
    std::size_t dim_,bond_types_=0;
    bool potts_,cluster_=false,ferro_=true,antiferro_=true;
    unsigned q_;
    double beta_,error_limit_=0,cluster_projection_=0,staggered_projection_=0;
    uint64_t production_,thermalization_,updates_=0,warmup_updates_=0,warmup_sites_=0,cluster_size_=0;
    Eigen::VectorXd field_,direction_;
    Eigen::MatrixXd spins_;
    std::vector<bond> edges_;
    std::vector<std::vector<std::size_t>> neighbors_;
    std::vector<double> factors_;
    std::vector<Eigen::MatrixXd> onsite_;
};

// All nonlinear estimates use aligned direct physical moments and weighted
// native jackknife propagation, including their cross-observable covariance.
inline results_type derive(results_type const& raw,alps::params const& p) {
    auto results=raw;
    for (auto const& entry:raw) {
        auto const& r=entry.second;
        if (!r.valid() || !r.store().batch().allFinite() || (r.count() && !r.mean().allFinite())
                || ((r.store().count().array()>0).count()>1 && !r.stderror().allFinite()))
            throw std::overflow_error("spinmc moments or uncertainties are not representable");
    }
    if (raw.empty() || (raw.at("Energy").store().count().array()>0).count()<2) return results;
    double beta=inverse_temperature(p), n=raw.at("Number of Sites").mean()(0);
    auto difference=native_mc::moment_difference(raw.at("Energy").count(),raw.at("Energy").num_batches());
    struct function : alps::alea::transformer<double> {
        std::size_t size; std::function<double(alps::alea::column<double> const&)> evaluate;
        function(std::size_t n,decltype(evaluate) f):size(n),evaluate(std::move(f)){}
        std::size_t in_size()const override{return size;} std::size_t out_size()const override{return 1;}
        alps::alea::column<double> operator()(alps::alea::column<double> const& x)const override{return {evaluate(x)};}
    };
    auto transform=[&](std::string const& name,std::initializer_list<const char*> keys,auto fn) {
        auto it=keys.begin(); auto inputs=raw.at(*it++);
        for (;it!=keys.end();++it) inputs=alps::alea::join(inputs,raw.at(*it));
        function tf(inputs.size(),fn);
        auto result=alps::alea::transform(alps::alea::jackknife_prop{},tf,inputs);
        if (result.store().batch().allFinite()) {
            if (!result.mean().allFinite() || !result.stderror().allFinite()) throw std::overflow_error("spinmc derived uncertainty is not representable: "+name);
            results.emplace(name,std::move(result));
        }
    };
    transform("Specific Heat",{"Energy","Energy^2"},[=](auto const& x){return beta==0 ? 0. : beta*(beta*difference(x(1),x(0)*x(0)))/n;});
    transform("Binder Cumulant U2",{"|Magnetization|","Magnetization^2"},[](auto const& x){return x(1)/(x(0)*x(0));});
    transform("Binder Cumulant",{"Magnetization^2","Magnetization^4"},[](auto const& x){return x(1)/(x(0)*x(0));});
    transform("Connected Susceptibility",{"|Magnetization|","Magnetization^2"},[=](auto const& x){return beta==0 ? 0. : beta*n*difference(x(1),x(0)*x(0));});
    transform("Magnetization^2 slope",{"Energy","Magnetization^2","E.Magnetization^2"},[=](auto const& x){return beta==0 ? 0. : beta*(beta*difference(x(2),x(0)*x(1)));});
    transform("Magnetization^4 slope",{"Energy","Magnetization^4","E.Magnetization^4"},[=](auto const& x){return beta==0 ? 0. : beta*(beta*difference(x(2),x(0)*x(1)));});
    transform("Binder Cumulant slope",{"Energy","Magnetization^2","Magnetization^4","E.Magnetization^2","E.Magnetization^4"},
        [=](auto const& x){return beta==0 && x(1)!=0 ? 0. : beta*(beta*difference(difference(x(4),x(0)*x(2))/(x(1)*x(1)),2*x(2)*difference(x(3),x(0)*x(1))/(x(1)*x(1)*x(1))));});
    return results;
}
}
