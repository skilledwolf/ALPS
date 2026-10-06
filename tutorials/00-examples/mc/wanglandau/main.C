// Copyright (C) 1997-2011 Synge Todo; 2026 ALPS Collaboration.
// SPDX-License-Identifier: MIT
#include "native.hpp"
#include "schema.hpp"

int main(int argc,char** argv) {
    return alps::mc::main<wl::simulation>(argc,argv,"wanglandau",wl_schema,{},
        [](alps::params& p,alps::run_configuration const& run) {
            for (auto key:{"ENERGY_WALK_RANGE","ENERGY_MEASURE_RANGE"})
                if (!p.exists(key)) p[key]=p["ENERGY_RANGE"];
            for (auto key:{"ENERGY_WALK_RANGE","ENERGY_MEASURE_RANGE"}) wl::width(p[key].as<wl::interval>());
            auto mode=p["MODE"].as<std::string>();
            if (mode=="learn") {
                if (run.input.exists("weights") || run.input.exists("measurements")) throw std::invalid_argument("Learning starts with flat weights; resume with input.checkpoint");
                if (p["INITIAL_UPDATE_FACTOR"].as<double>()<p["FINAL_UPDATE_FACTOR"].as<double>() || p["FINAL_UPDATE_FACTOR"].as<double>()<=1 || p["FLATNESS_THRESHOLD"].as<double>()<=0)
                    throw std::invalid_argument("Learning requires INITIAL_UPDATE_FACTOR >= FINAL_UPDATE_FACTOR > 1 and positive flatness");
            } else {
                const std::string key=mode=="measure" ? "weights" : "measurements";
                if (run.input.exists(mode=="measure" ? "measurements" : "weights") || !run.input.exists(key)) throw std::invalid_argument("This mode requires only input."+key);
                p["INPUT_FILES"]=run.input[key];
                auto files=p["INPUT_FILES"].as<std::vector<std::string>>();std::set<std::filesystem::path> unique;
                for (auto const& file:files) unique.insert(std::filesystem::weakly_canonical(file));
                if (files.empty() || unique.size()!=files.size()) throw std::invalid_argument("Input file list must be nonempty and contain distinct files");
            }
            if (mode=="reweight") {
                if (run.execution["chains"].as<size_t>()!=1 || run.input.exists("checkpoint") || run.output.exists("checkpoint")) throw std::invalid_argument("Reweighting is one deterministic analysis, without simulation checkpoints");
                auto temperatures=p["TEMPERATURE_SET"].as<std::vector<double>>();
                if (temperatures.empty() || std::any_of(temperatures.begin(),temperatures.end(),[](double t){return !std::isfinite(t) || t<=0 || !std::isfinite(1/t);})) throw std::invalid_argument("Reweighting temperatures and inverse temperatures must be positive and finite");
                if (p.exists("REFERENCE_BIN")!=p.exists("REFERENCE_LOGG")) throw std::invalid_argument("Reference normalization requires both REFERENCE_BIN and REFERENCE_LOGG");
                if (p.exists("REFERENCE_BIN") && !wl::contains(p["ENERGY_MEASURE_RANGE"].as<wl::interval>(),p["REFERENCE_BIN"].as<int64_t>())) throw std::invalid_argument("Reference bin is outside the measured range");
            } else if (p.exists("TEMPERATURE_SET") || p.exists("REFERENCE_BIN") || p.exists("REFERENCE_LOGG")) throw std::invalid_argument("Temperatures and reference normalization belong to reweight mode");
        },[](auto const& run,auto const& chains,auto const& p){wl::publish(run,chains,p);});
}
