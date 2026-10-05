/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/mcbase.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/alea/checkpoint.hpp>

namespace alps {

    mcbase::mcbase(parameters_type const & parms, std::size_t seed_offset)
        : parameters(parms)
        , random((parameters.value_or("SEED", 42)) + seed_offset, parameters.value_or<std::string>("RNG", "mt19937"))
    {
    }

    mcbase::~mcbase() = default;

    void mcbase::save(boost::filesystem::path const & filename) const {
        alps::hdf5::save_checkpoint(filename, [this](alps::hdf5::archive& ar) {
            ar["/simulation/realizations/0/clones/0"] << *this;
        });
    }

    void mcbase::load(boost::filesystem::path const & filename) {
        alps::hdf5::archive ar(filename);
        try {
            ar["/simulation/realizations/0/clones/0"] >> *this;
        } catch (...) {
            // Borrowed callback views must not keep the file open after this
            // operation, including when a Python load hook raises.
            try { ar.close(); } catch (...) {}
            throw;
        }
        ar.close();
    }

    bool mcbase::run(std::function<bool ()> const & stop_callback) {
        bool stopped = false;
        while(!(stopped = stop_callback()) && fraction_completed() < 1.) {
            update();
            measure();
        }
        return !stopped;
    }

    // implement a nice keys(m) function
    mcbase::result_names_type mcbase::result_names() const {
        result_names_type names;
        for(observable_collection_type::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
            names.push_back(it->first);
        return names;
    }

    mcbase::results_type mcbase::collect_results() const {
        return collect_results(result_names());
    }

    mcbase::results_type mcbase::collect_results(result_names_type const & names) const {
        results_type partial_results;
        for (auto const& name : names) {
            auto const& value = measurements.at(name);
            partial_results.emplace(name,std::visit([&](auto const& handle) -> alea::result::variant_type {
                if (!handle) throw std::invalid_argument("null MC measurement: " + name);
                return handle->result();
            },value));
        }
        return partial_results;
    }

    void mcbase::save(alps::hdf5::archive & ar) const {
        ar["/parameters"] << parameters;
        ar.create_group("measurements");
        for (auto const& child : ar.list_children("measurements")) {
            auto const path = "measurements/" + child;
            if (ar.is_group(path)) ar.delete_group(path);
            else ar.delete_data(path);
        }
        alps::alea::hdf5_serializer serializer(ar, "measurements");
        for (auto const& entry : measurements) {
            std::visit([&](auto const& handle) {
                if (!handle) throw std::invalid_argument("null MC measurement: " + entry.first);
                alps::alea::serialize(serializer,ar.encode_segment(entry.first),*handle);
            },entry.second);
        }
        ar["checkpoint/engine"] << random;
    }

    namespace {
    template<class A> mcbase::observable_type read_accumulator(alea::deserializer& codec,std::string const& name) {
        auto value=std::make_shared<A>();
        alea::deserialize(codec,name,*value);
        return value;
    }
    mcbase::observable_type read_measurement(hdf5::archive& ar,alea::deserializer& codec,std::string const& name) {
        auto path="measurements/"+name;
        uint32_t kind; ar[path+"/@kind"] >> kind;
        auto field=kind==6 ? "/batch/sum" : kind==10 ? "/levels/0/value" : "/value";
        bool complex=ar.is_complex(path+field);
        using C=std::complex<double>;
        switch (kind) {
        case 6: return complex ? read_accumulator<alea::batch_acc<C>>(codec,name) : read_accumulator<alea::batch_acc<double>>(codec,name);
        case 7: return complex ? read_accumulator<alea::mean_acc<C>>(codec,name) : read_accumulator<alea::mean_acc<double>>(codec,name);
        case 8:
            if (complex && ar.extent(path+"/centered_moment").size()==3)
                return read_accumulator<alea::var_acc<C,alea::elliptic_var>>(codec,name);
            return complex ? read_accumulator<alea::var_acc<C>>(codec,name) : read_accumulator<alea::var_acc<double>>(codec,name);
        case 9:
            if (complex && ar.extent(path+"/centered_moment").size()==4)
                return read_accumulator<alea::cov_acc<C,alea::elliptic_var>>(codec,name);
            return complex ? read_accumulator<alea::cov_acc<C>>(codec,name) : read_accumulator<alea::cov_acc<double>>(codec,name);
        case 10: return complex ? read_accumulator<alea::autocorr_acc<C>>(codec,name) : read_accumulator<alea::autocorr_acc<double>>(codec,name);
        default: throw std::invalid_argument("Unsupported MC accumulator checkpoint kind");
        }
    }
    }

    void mcbase::load(alps::hdf5::archive & ar) {
        parameters_type restored_parameters;
        observable_collection_type restored_measurements;
        auto restored_random = random;
        ar["/parameters"] >> restored_parameters;
        alps::alea::hdf5_serializer serializer(ar, "measurements");
        for (auto const& child : ar.list_children("measurements")) {
            auto value = read_measurement(ar,serializer,child);
            auto const name = ar.decode_segment(child);
            auto known = measurements.find(name);
            if (known != measurements.end()) {
                auto size=[](auto const& handle) { return handle ? handle->size() : 0; };
                if (known->second.index()!=value.index()
                    || std::visit(size,known->second)!=std::visit(size,value)) throw alps::alea::size_mismatch();
            }
            if (!restored_measurements.emplace(name, std::move(value)).second)
                throw std::invalid_argument("duplicate MC checkpoint measurement: " + name);
        }
        for (auto const& known : measurements)
            if (restored_measurements.find(known.first) == restored_measurements.end())
                throw std::invalid_argument("missing MC checkpoint measurement: " + known.first);
        ar["checkpoint/engine"] >> restored_random;
        parameters = std::move(restored_parameters);
        measurements = std::move(restored_measurements);
        random = std::move(restored_random);
    }

}
