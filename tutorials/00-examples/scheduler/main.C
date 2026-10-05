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

#include "ising.h"

int main(int argc, char** argv) {
    auto prepare = [](alps::params& p, alps::run_configuration const&) {
        p["CORRELATIONS"] = bool(ISING_VARIANT!=2);
        if constexpr (ISING_VARIANT==2) {
            if (p.exists("LATTICE")==p.exists("GRAPH"))
                throw std::invalid_argument("Specify exactly one LATTICE or GRAPH");
        } else if (p.exists("LATTICE") || p.exists("GRAPH"))
            throw std::invalid_argument("Use ising2 for a lattice-library graph");
    };
    auto publish = [](alps::run_configuration const& run, auto const& chains, alps::params const& p) {
        std::vector<alps::mc::batch_results> results;
        for (auto const& chain : chains) results.push_back(chain->collect_results());
        const auto joint = alps::mc::pool(results).at("Moments");
        alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(),[&](auto& ar) {
            ar["/parameters"] << p; ar["/run_config"] << run;
            alps::alea::hdf5_serializer raw(ar,"/simulation"), output(ar,"/simulation/results");
            serialize(raw,"joint",joint);
            ar.create_group("/simulation/results");
            if (joint.observations()>1) {
                std::vector<std::string> names{"Energy","Magnetization"};
                if constexpr (ISING_VARIANT!=2) names.insert(names.end(),{"Magnetization^2","Magnetization^4","Correlations"});
                for (size_t i=0; i<names.size(); ++i) {
                    const size_t width = i==4 ? joint.size()-4 : 1;
                    Eigen::MatrixXd select = Eigen::MatrixXd::Zero(width,joint.size());
                    select.block(0,i,width,width).setIdentity();
                    const auto result = alps::alea::transform(alps::alea::jackknife_prop(),
                        alps::alea::linear_transformer<double>(select),joint);
                    serialize(output,ar.encode_segment(names[i]),result);
                    std::cout << names[i] << ": " << result.mean().transpose() << " +/- " << result.stderror().transpose() << '\n';
                }
            } else ar["/simulation/unavailable/Statistics"] << std::string("At least two effective batches are required");
            for (size_t id=0; id<chains.size(); ++id) {
                const auto path = "/simulation/realizations/0/clones/"+std::to_string(id);
                ar[path+"/completed_sweeps"] << chains[id]->completed_sweeps();
                alps::alea::hdf5_serializer diagnostics(ar,path+"/autocorrelation");
                serialize(diagnostics,"Moments",chains[id]->template measurement<alps::alea::autocorr_acc<double>>("Autocorrelation")->result());
            }
        });
    };
    return alps::mc::main<IsingSimulation>(argc,argv,ISING_COMMAND,ising_schema,{},prepare,publish);
}
