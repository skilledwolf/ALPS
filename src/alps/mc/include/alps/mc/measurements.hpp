// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/mcbase.hpp>
#include <alps/alea/checkpoint.hpp>
#include <alps/alea/hdf5.hpp>

namespace alps::mc {
using batch = alps::alea::batch_acc<double>;
using autocorr = alps::alea::autocorr_acc<double>;
inline std::string diagnostic(std::string const& name) { return "Autocorrelation: " + name; }
inline void add_measurement(alps::mcbase& sim, std::string const& name, size_t size, size_t bins) {
    auto& values=sim.get_measurements();
    values.emplace(name,std::make_shared<batch>(size,bins));
    values.emplace(diagnostic(name),std::make_shared<autocorr>(size));
}
template<class T> void record(alps::mcbase& sim, std::string const& name, T const& value) {
    auto sample=alps::alea::make_adapter(value);
    *sim.measurement(name)<<sample;
    *sim.template measurement<autocorr>(diagnostic(name))<<sample;
}
inline auto batch_names(alps::mcbase::observable_collection_type const& values) {
    alps::mcbase::result_names_type names;
    for (auto const& [name,value]:values)
        if (std::holds_alternative<std::shared_ptr<batch>>(value)) names.push_back(name);
    return names;
}
inline batch validate_measurements(alps::mcbase::observable_collection_type const& values,
        alps::hdf5::archive& ar, uint64_t samples, size_t bins) {
    if (ar.list_children("measurements").size()!=values.size())
        throw std::invalid_argument("Unexpected checkpoint measurements");
    alps::alea::hdf5_serializer codec(ar,"measurements");
    batch energy;
    for (auto const& [name,handle]:values) {
        if (auto original=std::get_if<std::shared_ptr<batch>>(&handle)) {
            batch value;
            alps::alea::deserialize(codec,ar.encode_segment(name),value);
            if (value.size()!=(*original)->size() || value.num_batches()!=bins
                    || value.current_batch_size()!=value.cursor().factor()
                    || value.count()!=samples || !value.store().batch().allFinite())
                throw std::invalid_argument("Invalid checkpoint batch measurements");
            if (name=="Energy") energy=std::move(value);
        } else {
            autocorr value;
            alps::alea::deserialize(codec,ar.encode_segment(name),value);
            if (value.size()!=std::get<std::shared_ptr<autocorr>>(handle)->size()
                    || value.batch_size()!=1 || value.granularity()!=2 || value.count()!=samples)
                throw std::invalid_argument("Invalid checkpoint autocorrelation measurements");
        }
    }
    return energy;
}
// Retain each chain's chronology independently; pooling bins creates no chronology.
template<class Simulation>
void save_diagnostics(Simulation const& sim, alps::hdf5::archive& ar, std::string const& path,
                             alps::mcbase::result_names_type const& names) {
    alps::alea::hdf5_serializer series(ar,path+"/series"), analysis(ar,path+"/autocorrelation");
    for (auto const& name:names) {
        auto key=ar.encode_segment(name);
        alps::alea::serialize(series,key,*sim.measurement(name));
        alps::alea::serialize(analysis,key,sim.template measurement<autocorr>(diagnostic(name))->result());
    }
}
}
