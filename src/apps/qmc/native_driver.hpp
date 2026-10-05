// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "../mc/native_driver.hpp"
#include "analysis.hpp"

namespace native_qmc {

struct no_extra_results {
    template<class... T> void operator()(T const&...) const {}
};

template<class Chains,class Finalize=no_extra_results>
void publish_into(alps::hdf5::archive& ar,alps::run_configuration const& run,Chains const& chains,
                  std::string const& root,Finalize finalize={}) {
    using result=alps::alea::batch_result<double>;
    std::map<std::string,result> pooled;
    double reference=0;
    for (auto const& chain:chains) if (std::isfinite(chain->density_reference())) { reference=chain->density_reference(); break; }
    for (auto const& name:native_mc::batch_names(chains.front()->get_measurements())) {
        std::vector<result> values;
        for (auto const& chain:chains) {
            auto value=chain->measurement(name)->result();
            if (name=="Centered Density Moments" && value.count()) {
                double delta=chain->density_reference()-reference;
                auto& data=value.store().batch();
                data.row(1)+=2*delta*data.row(0)+delta*delta*data.row(2);
                data.row(0)+=delta*data.row(2);
            }
            values.push_back(std::move(value));
        }
        pooled.emplace(name,alps::alea::merge(values));
    }
    auto results=pooled;
    std::map<std::string,std::string> unavailable;
    for (auto const& name:chains.front()->signed_measurements()) {
        auto const& raw=pooled.at(name);
        if ((raw.store().count().array()>0).count()<2) {
            unavailable[name]="At least two occupied batches are required for sign reweighting";
            results.erase(name);
            continue;
        }
        auto value=alps::alea::transform(alps::alea::jackknife_prop{},divide_sign(raw.size()),raw);
        if (!value.mean().allFinite() || !value.stderror().allFinite()) {
            unavailable[name]="Sign reweighting has an undefined value or uncertainty";
            results.erase(name);
        } else results[name]=std::move(value);
    }
    derive(results,run.parameters,chains.front()->num_sites(),unavailable);
    finalize(pooled,results,unavailable);
    alps::save_results(results,run.parameters,ar,root+"/results");
    ar[root+"/density_reference"] << reference;
    auto const& signed_names=chains.front()->signed_measurements();
    ar[root+"/signed_measurements"] << std::vector<std::string>(signed_names.begin(),signed_names.end());
    ar[root+"/number_of_sites"] << uint64_t(chains.front()->num_sites());
    for (auto const& [name,labels]:chains.front()->measurement_labels())
        if (results.count(name)) ar[root+"/results/"+ar.encode_segment(name)+"/labels"] << labels;
    for (auto const& [name,reason]:unavailable) ar[root+"/unavailable/"+ar.encode_segment(name)] << reason;
    for (size_t i=0;i<chains.size();++i) {
        auto path=root+"/realizations/0/clones/"+std::to_string(i);
        ar[path+"/sampling_parameters"] << chains[i]->sampling_parameters();
        ar[path+"/density_reference"] << chains[i]->density_reference();
        ar[path+"/completed_sweeps"] << chains[i]->completed_sweeps();
        ar[path+"/complete"] << (chains[i]->fraction_completed()>=1.);
        auto names=native_mc::batch_names(chains[i]->get_measurements());
        alps::save_results(chains[i]->template collect_results_as<result>(names),run.parameters,ar,path+"/results");
        native_mc::save_diagnostics(*chains[i],ar,path,names);
    }
}

template<class Chains,class Finalize=no_extra_results>
void publish(alps::run_configuration const& run,Chains const& chains,alps::params const&,Finalize finalize={}) {
    alps::hdf5::save_checkpoint(run.output["results"].as<std::string>(),[&](alps::hdf5::archive& ar) {
        publish_into(ar,run,chains,"/simulation",finalize);
        ar["/run_config"] << run;
    });
}

inline std::string schema(char const* application,char const* common,char const* specific) {
    auto schema=toml::parse(common), extra=toml::parse(specific);
    schema.insert("application",application);
    for (auto const& [name,value]:*extra["parameters"].as_table())
        schema["parameters"].as_table()->insert(name,value);
    std::ostringstream text; text<<schema;
    return text.str();
}

template<class Simulation>
int main(int argc,char** argv,char const* application,char const* common,char const* specific) {
    auto text=schema(application,common,specific);
    return native_mc::main<Simulation>(argc,argv,application,text.c_str(),
        [](std::string const&,toml::node const&)->char const*{return nullptr;},
        [](alps::params& p,alps::run_configuration const&){
            if (!p.exists("THERMALIZATION")) p["THERMALIZATION"]=p["SWEEPS"].as<uint64_t>()/10;
        }, [](auto const& run,auto const& chains,auto const& p){ publish(run,chains,p); });
}
}
