// Copyright (C) 1997-2011 Synge Todo; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#pragma once
#include <alps/lattice.h>
#include <alps/mc/driver.hpp>
#include <numeric>

namespace wl {
using interval=std::vector<int64_t>;
inline size_t width(interval const& r) {
    if (r.size()!=2 || r[0]>r[1] || r[0]<INT_MIN || r[1]>INT_MAX || r[1]-r[0]>=1000000)
        throw std::invalid_argument("Energy ranges require two ordered integer endpoints and fewer than one million bins");
    return size_t(r[1]-r[0]+1);
}
inline bool contains(interval const& r,int64_t e) {return e>=r[0] && e<=r[1];}

// A measured window distinguishes unvisited energies from arbitrary log(g)
// offsets. Zero or negative log(g) is a valid weight, not a missing-bin flag.
struct density {
    interval range;
    uint64_t sites=0;
    int64_t coupling=0;
    std::vector<uint64_t> bonds;
    std::vector<double> logg;
    std::vector<bool> visited;
    bool complete=false;
    void same_model(density const& other) const {
        if (sites!=other.sites || coupling!=other.coupling || bonds!=other.bonds)
            throw std::invalid_argument("Wang-Landau input Hamiltonians differ");
    }
    void save(alps::hdf5::archive& ar) const {
        ar["range"]<<range;ar["sites"]<<sites;ar["coupling"]<<coupling;ar["bonds"]<<bonds;
        ar["logg"]<<logg;ar["visited"]<<visited;ar["complete"]<<complete;
    }
    void load(alps::hdf5::archive& ar) {
        ar["range"]>>range;ar["sites"]>>sites;ar["coupling"]>>coupling;ar["bonds"]>>bonds;
        ar["logg"]>>logg;ar["visited"]>>visited;ar["complete"]>>complete;
        if (!sites || bonds.size()%2 || logg.size()!=width(range) || visited.size()!=logg.size()
                || !std::all_of(logg.begin(),logg.end(),[](double x){return std::isfinite(x);})
                || (complete && std::none_of(visited.begin(),visited.end(),[](bool x){return x;}))
                || std::any_of(bonds.begin(),bonds.end(),[&](uint64_t x){return x>=sites;}))
            throw std::invalid_argument("Invalid Wang-Landau density of states");
    }
    density crop(interval const& target) const {
        width(target);
        if (target[0]<range[0] || target[1]>range[1]) throw std::invalid_argument("Weights do not cover the measurement range");
        auto out=*this;out.range=target;
        size_t a=target[0]-range[0],b=a+width(target);
        out.logg={logg.begin()+a,logg.begin()+b};out.visited={visited.begin()+a,visited.begin()+b};
        return out;
    }
};

inline density stitch(std::vector<density> windows) {
    if (windows.empty()) throw std::invalid_argument("At least one weight window is required");
    std::stable_sort(windows.begin(),windows.end(),[](auto const& a,auto const& b){return a.range<b.range;});
    auto result=windows.front();
    std::vector<size_t> counts(result.logg.size(),1);
    for (size_t w=1;w<windows.size();++w) {
        auto const& next=windows[w];result.same_model(next);
        int64_t low=next.range[0],high=std::min(result.range[1],next.range[1]);
        double shift=0;size_t overlap=0;
        for (int64_t e=low;e<=high;++e) {
            size_t a=e-result.range[0],b=e-next.range[0];
            if (result.visited[a]!=next.visited[b]) throw std::invalid_argument("Weight windows disagree on accessible energies");
            if (next.visited[b]) {shift+=result.logg[a]-next.logg[b];++overlap;}
        }
        if (!overlap) throw std::invalid_argument("Weight windows require an overlapping visited energy");
        shift/=overlap;
        result.range[1]=std::max(result.range[1],next.range[1]);
        auto n=width(result.range);result.logg.resize(n);result.visited.resize(n);counts.resize(n,0);
        for (size_t b=0;b<next.logg.size();++b) if (next.visited[b]) {
            size_t a=next.range[0]+b-result.range[0];
            result.logg[a]+=(next.logg[b]+shift-result.logg[a])/++counts[a];
            result.visited[a]=true;
        }
        result.complete=result.complete && next.complete;
    }
    return result;
}

inline void validate_joint(alps::alea::batch_result<double> const& data,size_t energies) {
    if (data.size()!=4*energies || !data.store().batch().allFinite()) throw std::invalid_argument("Invalid microcanonical moments");
    for (size_t b=0;b<data.num_batches();++b) {
        const double count=data.store().count()(b),tolerance=1.e-10*count;
        double total=0;
        for (size_t i=0;i<energies;++i) {
            auto x=data.store().batch().col(b).segment(4*i,4);total+=x(0);
            if (x(0)<0 || std::abs(x(1))>x(0)+tolerance || x(2)<0 || x(2)>x(0)+tolerance || x(3)<0 || x(3)>x(2)+tolerance)
                throw std::invalid_argument("Microcanonical moments exceed physical bounds");
        }
        if (total>count+tolerance) throw std::invalid_argument("Microcanonical visits exceed sample count");
    }
}

class simulation {
    alps::params p_;
    size_t bins_,chain_;
    alps::random01 random_;
    std::string mode_;
    interval walk_,measure_;
    density weights_;
    std::vector<std::vector<size_t>> neighbors_;
    std::vector<int> spins_;
    std::vector<uint64_t> histogram_,overall_,ascending_,measured_;
    alps::mc::batch joint_,roundtrips_;
    uint64_t steps_=0,stage_=0;
    int64_t energy_=0;
    int direction_=0,last_direction_=0;
    bool done_=false;
    std::vector<alps::alea::batch_result<double>> input_moments_;
    uint64_t production() const {auto warm=p_["THERMALIZATION"].as<uint64_t>();return steps_>warm ? steps_-warm : 0;}
    double log_factor() const {return done_ ? 0 : std::ldexp(std::log(p_["INITIAL_UPDATE_FACTOR"].as<double>()),-int(stage_));}
    int64_t spin_energy() const {
        int64_t e=0;
        for (size_t b=0;b<weights_.bonds.size();b+=2) e-=weights_.coupling*spins_[weights_.bonds[b]]*spins_[weights_.bonds[b+1]];
        return e;
    }
    void refine() {
        uint64_t minimum=UINT64_MAX;long double total=0;size_t count=0;
        for (int64_t e=measure_[0];e<=measure_[1];++e) {
            size_t i=e-walk_[0];
            if (overall_[i]) {minimum=std::min(minimum,histogram_[i]);total+=histogram_[i];++count;}
        }
        if (!total || minimum<p_["FLATNESS_THRESHOLD"].as<double>()*total/count) return;
        if (log_factor()<=std::log(p_["FINAL_UPDATE_FACTOR"].as<double>())) done_=true;
        else {++stage_;std::fill(histogram_.begin(),histogram_.end(),0);}
    }
public:
    simulation(alps::params const& p,size_t bins,size_t chain):p_(p),bins_(bins),chain_(chain),
        random_(p["SEED"].as<int>()+chain,p["RNG"].as<std::string>()),mode_(p["MODE"].as<std::string>()),
        walk_(p["ENERGY_WALK_RANGE"].as<interval>()),measure_(p["ENERGY_MEASURE_RANGE"].as<interval>()),
        joint_(4*width(measure_),bins),roundtrips_(1,bins) {
        auto size=width(walk_);
        if (measure_[0]<walk_[0] || measure_[1]>walk_[1]) throw std::invalid_argument("Measurement range must be inside the walk range");
        alps::graph_helper<> graph(p);
        weights_.sites=graph.num_sites();weights_.coupling=p["COUPLING"].as<int64_t>();weights_.range=walk_;
        if (!weights_.sites) throw std::invalid_argument("The lattice must contain sites");
        neighbors_.resize(weights_.sites);
        for (auto [it,end]=graph.bonds();it!=end;++it) {
            size_t a=graph.source(*it),b=graph.target(*it);
            weights_.bonds.insert(weights_.bonds.end(),{a,b});
            if (a!=b) {neighbors_[a].push_back(b);neighbors_[b].push_back(a);}
        }
        const long double bound=std::abs(double(weights_.coupling))*graph.num_bonds();
        if (bound>INT_MAX || walk_[0]>bound || walk_[1]<-bound) throw std::invalid_argument("Energy range is outside the supported Hamiltonian bounds");
        for (size_t i=0;i<weights_.sites;++i) spins_.push_back(random_()<.5 ? 1 : -1);
        energy_=spin_energy();weights_.logg.resize(size);weights_.visited.resize(size);
        histogram_.resize(size);overall_.resize(size);ascending_.resize(size);measured_.resize(size);
        if (mode_=="learn") return;
        std::vector<density> windows;
        for (auto const& file:p["INPUT_FILES"].as<std::vector<std::string>>()) {
            alps::hdf5::archive ar(file);density window;ar["/weights"]>>window;weights_.same_model(window);
            if (!window.complete) throw std::invalid_argument("Input weights are not complete");
            windows.push_back(window);
            if (mode_=="reweight") {
                bool complete;ar["/microcanonical/complete"]>>complete;
                if (!complete || window.range!=measure_) throw std::invalid_argument("Reweighting requires complete measurements on the same energy range");
                alps::alea::batch_result<double> moments;
                alps::alea::hdf5_serializer codec(ar,"/microcanonical");deserialize(codec,"joint",moments);
                validate_joint(moments,width(measure_));
                for (size_t i=0;i<window.visited.size();++i)
                    if (!window.visited[i] && moments.store().batch().row(4*i).sum()!=0) throw std::invalid_argument("Measured energy is absent from the density of states");
                input_moments_.push_back(std::move(moments));
                if (windows.size()>1) {
                    auto const& first=windows.front();double shift=0;bool anchored=false;
                    for (size_t i=0;i<window.logg.size();++i) {
                        if (first.visited[i]!=window.visited[i]) throw std::invalid_argument("Measurement weights differ");
                        if (!window.visited[i]) continue;
                        double delta=window.logg[i]-first.logg[i];
                        if (!anchored) {shift=delta;anchored=true;}
                        if (std::abs(delta-shift)>1.e-10*(1+std::abs(delta)+std::abs(shift))) throw std::invalid_argument("Measurement weights differ");
                    }
                }
            }
        }
        auto loaded=stitch(std::move(windows)).crop(measure_);
        if (p.exists("REFERENCE_BIN") && !loaded.visited[p["REFERENCE_BIN"].as<int64_t>()-measure_[0]]) throw std::invalid_argument("Reference normalization requires an accessible energy");
        for (size_t i=0;i<loaded.logg.size();++i) {
            size_t j=measure_[0]+i-walk_[0];weights_.logg[j]=loaded.logg[i];weights_.visited[j]=loaded.visited[i];
        }
        weights_.complete=true;
    }
    static alps::params checkpoint_parameters(alps::params p) {p.erase("SWEEPS");return p;}
    uint64_t completed_sweeps() const {return steps_;}
    double fraction_completed() const {return mode_=="reweight" ? 1 : mode_=="learn" ? double(done_) : double(production())/p_["SWEEPS"].as<uint64_t>();}
    auto moments() const {return mode_=="reweight" ? alps::alea::merge(input_moments_) : joint_.result();}
    auto roundtrips() const {return roundtrips_.result();}
    auto const& histogram() const {return histogram_;}
    auto const& overall() const {return overall_;}
    auto const& ascending() const {return ascending_;}
    auto const& measured() const {return measured_;}
    auto const& walk_range() const {return walk_;}
    auto density_of_states() const {
        auto result=weights_.crop(measure_);
        if (mode_=="learn") {
            result.complete=done_;
            for (size_t i=0;i<result.visited.size();++i) result.visited[i]=overall_[measure_[0]+i-walk_[0]]>0;
        }
        return result;
    }
    void update() {
        ++steps_;
        for (size_t s=0;s<spins_.size();++s) {
            // A symmetric keep/flip proposal avoids the all-accepted sweep's
            // even-parity trap, including zero coupling and flat weights.
            if (random_()<.5) {
                int64_t next=energy_;
                for (auto j:neighbors_[s]) next+=2*weights_.coupling*spins_[s]*spins_[j];
                bool accept=energy_<walk_[0] ? next>=energy_ && next<=walk_[1]
                    : energy_>walk_[1] ? next<=energy_ && next>=walk_[0]
                    : contains(walk_,next) && std::log(random_())<weights_.logg[energy_-walk_[0]]-weights_.logg[next-walk_[0]];
                if (accept) {spins_[s]*=-1;energy_=next;}
            }
            if (!contains(walk_,energy_)) continue;
            size_t i=energy_-walk_[0];++histogram_[i];++overall_[i];
            if (energy_==walk_[0]) direction_=1;else if (energy_==walk_[1]) direction_=2;
            if (direction_==1) ++ascending_[i];
            if (mode_=="learn") weights_.logg[i]+=log_factor()*(contains(measure_,energy_) ? 1 : p_["VISIT_PENALTY"].as<double>());
            else if (contains(measure_,energy_) && !weights_.visited[i]) throw std::invalid_argument("Encountered an energy missing from input weights");
        }
        if (mode_=="learn" || production()) {
            if (contains(walk_,energy_)) {
                roundtrips_<<alps::alea::column<double>{double(last_direction_==1 && direction_==2)};
                last_direction_=direction_;
            }
            if (contains(measure_,energy_)) ++measured_[energy_-walk_[0]];
        }
        if (mode_=="learn") {
            if (steps_%p_["CHECK_INTERVAL"].as<uint64_t>()==0) refine();
        } else if (production()) {
            auto sample=alps::alea::column<double>::Zero(joint_.size()).eval();
            if (contains(measure_,energy_)) {
                double m=double(std::accumulate(spins_.begin(),spins_.end(),int64_t(0)))/spins_.size();
                sample.segment(4*(energy_-measure_[0]),4)<<1,m,m*m,m*m*m*m;
            }
            joint_<<sample;
        }
    }
    void measure() {}
    void save(alps::hdf5::archive& ar) const {
        ar["parameters"]<<p_;ar["chain"]<<uint64_t(chain_);ar["rng"]<<random_;
        ar["spins"]<<spins_;ar["energy"]<<energy_;ar["weights"]<<weights_;
        ar["steps"]<<steps_;ar["stage"]<<stage_;ar["done"]<<done_;
        ar["direction"]<<direction_;ar["last_direction"]<<last_direction_;
        ar["histogram"]<<histogram_;ar["overall"]<<overall_;ar["ascending"]<<ascending_;ar["measured"]<<measured_;
        alps::alea::hdf5_serializer codec(ar,"measurements");serialize(codec,"joint",joint_);serialize(codec,"roundtrips",roundtrips_);
    }
    void load(alps::hdf5::archive& ar) {
        auto restored=*this;alps::params parameters;uint64_t chain;
        ar["parameters"]>>parameters;ar["chain"]>>chain;
        if (checkpoint_parameters(parameters)!=checkpoint_parameters(p_) || chain!=chain_) throw std::invalid_argument("Wang-Landau checkpoint parameters differ");
        ar["rng"]>>restored.random_;ar["spins"]>>restored.spins_;ar["energy"]>>restored.energy_;
        auto& w=restored.weights_;ar["weights"]>>w;weights_.same_model(w);
        ar["steps"]>>restored.steps_;ar["stage"]>>restored.stage_;ar["done"]>>restored.done_;
        ar["direction"]>>restored.direction_;ar["last_direction"]>>restored.last_direction_;
        ar["histogram"]>>restored.histogram_;ar["overall"]>>restored.overall_;ar["ascending"]>>restored.ascending_;ar["measured"]>>restored.measured_;
        const size_t n=width(walk_);
        if (w.range!=walk_ || w.complete!=weights_.complete || w.visited!=weights_.visited || restored.spins_.size()!=spins_.size()
                || restored.random_.name()!=random_.name() || restored.stage_>64 || restored.direction_<0 || restored.direction_>2
                || restored.last_direction_<0 || restored.last_direction_>2
                || !std::all_of(restored.spins_.begin(),restored.spins_.end(),[](int x){return x==1 || x==-1;})
                || restored.histogram_.size()!=n || restored.overall_.size()!=n || restored.ascending_.size()!=n || restored.measured_.size()!=n)
            throw std::invalid_argument("Invalid Wang-Landau checkpoint state");
        uint64_t final_stage=0;
        if (mode_=="learn") while (std::ldexp(std::log(p_["INITIAL_UPDATE_FACTOR"].as<double>()),-int(final_stage))>std::log(p_["FINAL_UPDATE_FACTOR"].as<double>())) ++final_stage;
        if (restored.energy_!=restored.spin_energy()
                || (mode_!="learn" && (w.logg!=weights_.logg || restored.stage_ || restored.done_ || restored.production()>p_["SWEEPS"].as<uint64_t>()))
                || (mode_=="learn" && (restored.stage_>final_stage || restored.stage_>restored.steps_/p_["CHECK_INTERVAL"].as<uint64_t>()
                    || (restored.done_ && (restored.stage_!=final_stage || restored.steps_%p_["CHECK_INTERVAL"].as<uint64_t>())))))
            throw std::invalid_argument("Invalid Wang-Landau checkpoint energy, weights or refinement stage");
        long double visits=0,measurements=0;
        for (size_t i=0;i<n;++i) {
            visits+=restored.overall_[i];measurements+=restored.measured_[i];
            if (restored.histogram_[i]>restored.overall_[i] || restored.ascending_[i]>restored.overall_[i]
                    || (!contains(measure_,walk_[0]+i) && restored.measured_[i])) throw std::invalid_argument("Invalid Wang-Landau histograms");
        }
        const uint64_t samples=mode_=="learn" ? restored.steps_ : restored.production();
        if (visits>static_cast<long double>(restored.steps_)*spins_.size() || measurements>samples) throw std::invalid_argument("Invalid Wang-Landau visit count");
        alps::alea::hdf5_serializer codec(ar,"measurements");deserialize(codec,"joint",restored.joint_);deserialize(codec,"roundtrips",restored.roundtrips_);
        for (auto const* value:{&restored.joint_,&restored.roundtrips_})
            if (value->num_batches()!=bins_ || value->current_batch_size()!=value->cursor().factor() || !value->store().batch().allFinite()) throw std::invalid_argument("Invalid Wang-Landau statistical batches");
        if (restored.joint_.count()!=(mode_=="learn" ? 0 : samples) || restored.roundtrips_.size()!=1 || restored.roundtrips_.count()>samples)
            throw std::invalid_argument("Invalid Wang-Landau statistical counts");
        validate_joint(restored.joint_.result(),width(measure_));
        if (mode_=="measure") for (size_t i=0;i<width(measure_);++i)
            if (restored.joint_.store().batch().row(4*i).sum()!=restored.measured_[measure_[0]+i-walk_[0]]) throw std::invalid_argument("Measurement histogram and batches disagree");
        if ((restored.roundtrips_.store().batch().array()<0).any() ||
                (restored.roundtrips_.store().batch().row(0).array()>restored.roundtrips_.store().count().cast<double>().array()).any()) throw std::invalid_argument("Invalid round-trip counts");
        *this=std::move(restored);
    }
};

inline void write_estimates(alps::hdf5::archive& ar,std::string const& path,alps::mc::batch_results const& results,alps::mc::unavailable_results const& unavailable) {
    alps::alea::hdf5_serializer codec(ar,path+"/results");ar.create_group(path+"/results");
    for (auto const& [name,value]:results) serialize(codec,ar.encode_segment(name),value);
    for (auto const& [name,reason]:unavailable) ar[path+"/unavailable/"+ar.encode_segment(name)]<<reason;
}

inline void reweight(alps::hdf5::archive& ar,density const& weights,alps::alea::batch_result<double> const& joint,alps::params const& p) {
    auto temperatures=p["TEMPERATURE_SET"].as<std::vector<double>>();
    for (size_t t=0;t<temperatures.size();++t) {
        const long double beta=1/temperatures[t],origin=weights.range[0];
        long double scale=-std::numeric_limits<long double>::infinity();
        for (size_t i=0;i<weights.logg.size();++i) if (weights.visited[i]) scale=std::max(scale,weights.logg[i]-beta*i);
        std::vector<long double> w(weights.logg.size());
        for (size_t i=0;i<w.size();++i) if (weights.visited[i]) w[i]=std::exp(weights.logg[i]-beta*i-scale);
        auto evaluate=[&](auto const& x,std::string const& name)->double {
            long double z=0,e=0,e2=0,m=0,m2=0,m4=0;
            for (size_t i=0;i<w.size();++i) {
                auto q=w[i]*x(4*i);z+=q;e+=q*i;e2+=q*i*i;
                m+=w[i]*x(4*i+1);m2+=w[i]*x(4*i+2);m4+=w[i]*x(4*i+3);
            }
            if (z<=0) return NAN;
            e/=z;e2/=z;m/=z;m2/=z;m4/=z;
            if (name=="Energy") return origin+e;
            if (name=="Energy^2") return origin*origin+2*origin*e+e2;
            if (name=="Energy Density") return (origin+e)/weights.sites;
            if (name=="Specific Heat") return beta*beta*std::max(0.L,e2-e*e)/weights.sites;
            if (name=="Magnetization") return m;
            if (name=="Magnetization^2") return m2;
            if (name=="Magnetization^4") return m4;
            if (name=="Binder Ratio of Magnetization") return m4>0 ? m2*m2/m4 : NAN;
            size_t ref=p["REFERENCE_BIN"].as<int64_t>()-weights.range[0];
            if (x(4*ref)<=0) return NAN;
            const long double logz=std::log(z)+scale-beta*origin+p["REFERENCE_LOGG"].as<double>()-weights.logg[ref]-std::log(x(4*ref));
            return name=="Free-Energy" ? -logz/beta : name=="Free-Energy Density" ? -logz/beta/weights.sites
                : name=="Entropy" ? logz+beta*(origin+e) : (logz+beta*(origin+e))/weights.sites;
        };
        alps::mc::batch_results results;alps::mc::unavailable_results unavailable;
        alps::mc::estimate(results,&unavailable,"Number of Sites",joint,[&](auto const&){return double(weights.sites);});
        alps::mc::estimate(results,&unavailable,"Temperature",joint,[&](auto const&){return temperatures[t];});
        alps::mc::estimate(results,&unavailable,"Inverse Temperature",joint,[&](auto const&){return double(beta);});
        for (std::string name:{"Energy","Energy^2","Energy Density","Specific Heat","Magnetization","Magnetization^2","Magnetization^4","Binder Ratio of Magnetization",
                              "Free-Energy","Free-Energy Density","Entropy","Entropy Density"}) {
            if ((name.find("Free-Energy")==0 || name.find("Entropy")==0) && !p.exists("REFERENCE_BIN")) unavailable[name]="Absolute normalization requires REFERENCE_BIN and REFERENCE_LOGG";
            else alps::mc::estimate(results,&unavailable,name,joint,[&](auto const& x){return evaluate(x,name);});
        }
        const auto path="/simulation/replicas/"+std::to_string(t);auto parameters=p;parameters["T"]=temperatures[t];parameters["BETA"]=double(beta);
        ar[path+"/parameters"]<<parameters;write_estimates(ar,path,results,unavailable);
    }
}

template<class Chains> void publish(alps::run_configuration const& run,Chains const& chains,alps::params const& p) {
    alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(),[&](auto& ar) {
        ar["/run_config"]<<run;ar["/parameters"]<<p;
        std::vector<density> windows;std::vector<alps::alea::batch_result<double>> moments,roundtrips;
        bool complete=true;
        for (size_t i=0;i<chains.size();++i) {
            auto const& chain=*chains[i];windows.push_back(chain.density_of_states());moments.push_back(chain.moments());roundtrips.push_back(chain.roundtrips());
            complete=complete && chain.fraction_completed()>=1;
            const auto path="/simulation/realizations/0/clones/"+std::to_string(i);
            ar[path+"/completed_sweeps"]<<chain.completed_sweeps();ar[path+"/weights"]<<windows.back();
            ar[path+"/walk_range"]<<chain.walk_range();ar[path+"/Final Histogram"]<<chain.histogram();ar[path+"/Overall Histogram"]<<chain.overall();
            ar[path+"/Histogram of Ascending Walker"]<<chain.ascending();ar[path+"/Measurement Histogram"]<<chain.measured();
        }
        // Unfinished learning windows may have different coverage. Retain each
        // independently; only a completed estimate may be used as new weights.
        auto weights=complete || p["MODE"].as<std::string>()!="learn" ? stitch(windows) : windows.front();weights.complete=weights.complete && complete;
        ar["/weights"]<<weights;
        if (p["MODE"].as<std::string>()=="reweight") {reweight(ar,weights,alps::alea::merge(moments),p);return;}
        alps::mc::batch_results diagnostics;alps::mc::unavailable_results unavailable;
        auto trips=alps::alea::merge(roundtrips);
        alps::mc::estimate(diagnostics,&unavailable,"Inverse Round-trip Time",trips,[](auto const& x){return x(0);});
        alps::mc::estimate(diagnostics,&unavailable,"Round-trip Time",trips,[](auto const& x){return x(0)>0 ? 1/x(0) : NAN;});
        write_estimates(ar,"/simulation",diagnostics,unavailable);
        if (p["MODE"].as<std::string>()=="learn") return;
        auto joint=alps::alea::merge(moments);alps::alea::hdf5_serializer codec(ar,"/microcanonical");serialize(codec,"joint",joint);ar["/microcanonical/complete"]<<complete;
        for (size_t i=0;i<weights.logg.size();++i) if (weights.visited[i]) {
            alps::mc::batch_results results;alps::mc::unavailable_results missing;
            alps::mc::estimate(results,&missing,"Number of Sites",joint,[&](auto const&){return double(weights.sites);});
            alps::mc::estimate(results,&missing,"Energy",joint,[&](auto const&){return double(weights.range[0]+int64_t(i));});
            for (size_t k=1;k<4;++k) alps::mc::estimate(results,&missing,std::string("Magnetization")+(k==1 ? "" : "^"+std::to_string(k==2 ? 2 : 4)),joint,
                [=](auto const& x){return x(4*i)>0 ? x(4*i+k)/x(4*i) : NAN;});
            const auto path="/simulation/replicas/"+std::to_string(i);auto parameters=p;parameters["Energy"]=weights.range[0]+int64_t(i);parameters["Number of Sites"]=weights.sites;
            ar[path+"/parameters"]<<parameters;write_estimates(ar,path,results,missing);
        }
    });
}
} // namespace wl
