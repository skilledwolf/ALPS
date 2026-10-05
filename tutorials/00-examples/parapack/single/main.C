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

#include "native.hpp"

int main(int argc,char** argv) {
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
    return alps::mc::main<single_ising>(argc,argv,"ising_single",ising_schema,{},prepare,publish);
}
