// Copyright (C) 2011 Synge Todo; 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/fortran/fortran_wrapper.h>
#include <alps/fortran/fwrapper_impl.h>
#include <alps/mc/driver.hpp>
#include <algorithm>
#include <cstring>

namespace alps {
fortran_wrapper::fortran_wrapper(params const& p, std::size_t bins, std::size_t chain)
    : mcbase(p, chain), bins_(bins), chain_(chain) {
    try { alps_init(this); check(); alps_init_observables(this); check(); }
    catch (...) { alps_finalize(this); throw; }
}
fortran_wrapper::~fortran_wrapper() { alps_finalize(this); }
void fortran_wrapper::update() {
    int thermalized = 0;
    alps_is_thermalized(&thermalized, this); check();
    collecting_ = thermalized != 0;
    alps_run(this); check();
    collecting_ = false;
    ++steps_;
}
double fortran_wrapper::fraction_completed() const {
    double value = 0;
    alps_progress(&value, const_cast<fortran_wrapper*>(this)); check();
    if (!std::isfinite(value) || value < 0) throw std::invalid_argument("Invalid Fortran progress");
    return value;
}
std::string fortran_wrapper::field(int type, std::size_t count, bool writing) {
    if (!archive_ || writing != writing_) throw std::logic_error("Fortran checkpoint call outside save/load");
    auto path = "checkpoint/fortran/fields/" + std::to_string(field_++);
    if (writing) { archive_->create_group(path); (*archive_)[path+"/@type"] << type; (*archive_)[path+"/@count"] << count; }
    else {
        int saved_type; std::size_t saved_count;
        (*archive_)[path+"/@type"] >> saved_type; (*archive_)[path+"/@count"] >> saved_count;
        if (saved_type != type || saved_count != count) throw std::invalid_argument("Fortran checkpoint field mismatch");
    }
    return path+"/value";
}
void fortran_wrapper::save(hdf5::archive& ar) const {
    auto& self = *const_cast<fortran_wrapper*>(this);
    mcbase::save(ar);
    ar["checkpoint/fortran/version"] << 1;
    ar["checkpoint/fortran/chain"] << chain_;
    ar["checkpoint/fortran/steps"] << steps_;
    self.archive_ = &ar; self.field_ = 0; self.writing_ = true;
    alps_save(&self);
    self.archive_ = nullptr; check();
    ar["checkpoint/fortran/fields_count"] << self.field_;
}
void fortran_wrapper::load(hdf5::archive& ar) {
    fortran_wrapper staged(parameters, bins_, chain_);
    int version; std::size_t chain, fields;
    params saved;
    ar["checkpoint/fortran/version"] >> version;
    ar["checkpoint/fortran/chain"] >> chain;
    ar["/parameters"] >> saved;
    if (version != 1 || chain != chain_ || checkpoint_parameters(saved) != checkpoint_parameters(parameters))
        throw std::invalid_argument("Incompatible Fortran checkpoint");
    staged.mcbase::load(ar);
    staged.parameters = parameters; // Retain the requested production length.
    ar["checkpoint/fortran/steps"] >> staged.steps_;
    ar["checkpoint/fortran/fields_count"] >> fields;
    if (ar.list_children("measurements").size() != measurements.size())
        throw std::invalid_argument("Unexpected Fortran measurements");
    for (auto const& name : mc::batch_names(staged.measurements)) {
        auto const& batch = *staged.measurement(name);
        if (batch.num_batches()!=bins_ || batch.count()!=staged.measurement<mc::autocorr>(mc::diagnostic(name))->count())
            throw std::invalid_argument("Inconsistent Fortran checkpoint measurements");
    }
    staged.archive_ = &ar;
    alps_load(&staged);
    staged.archive_ = nullptr; staged.check();
    if (staged.field_ != fields) throw std::invalid_argument("Unconsumed Fortran checkpoint fields");
    staged.fraction_completed();
    std::swap(context, staged.context);
    std::swap(random, staged.random);
    measurements.swap(staged.measurements);
    steps_ = staged.steps_;
}
int fortran_wrapper::main(int argc, char** argv, char const* application, char const* schema) {
    return mc::main<fortran_wrapper>(argc, argv, application, schema, {},
        [](auto&, auto const&) {}, [](auto const& run, auto const& chains, auto const&) {
        mc::batch_results results;
        auto names = mc::batch_names(chains.front()->get_measurements());
        // User observables may have different sampling intervals.
        for (auto const& name : names) {
            std::vector<alea::batch_result<double>> parts;
            for (auto const& chain : chains) parts.push_back(chain->measurement(name)->result());
            results.emplace(name, alea::merge(parts));
        }
        hdf5::save_checkpoint(run.output["results"].template as<std::string>(), [&](hdf5::archive& ar) {
            save_results(results, run.parameters, ar, "/simulation/results");
            ar["/run_config"] << run;
            for (std::size_t i=0; i<chains.size(); ++i) {
                auto path = "/simulation/realizations/0/clones/" + std::to_string(i);
                ar[path+"/completed_sweeps"] << chains[i]->completed_sweeps();
                mc::save_diagnostics(*chains[i], ar, path, names);
            }
        });
    });
}
}

namespace {
using worker = alps::fortran_wrapper;
std::string text(char const* data, std::size_t width) {
    std::string value(data, width);
    value.erase(value.find_last_not_of(' ')+1);
    return value;
}
void put(char* data, std::size_t width, std::string const& value) {
    if (value.size()>width) throw std::invalid_argument("Fortran character buffer too small");
    std::fill_n(data, width, ' '); std::copy(value.begin(), value.end(), data);
}
template<class F> void numeric(int type, F f) {
    switch (type) {
    case 1: f(int{}); break;
    case 2: f(std::int64_t{}); break;
    case 3: f(float{}); break;
    case 4: f(double{}); break;
    default: throw std::invalid_argument("Invalid Fortran numeric type");
    }
}
void transfer(worker& w, void* data, std::size_t count, int type, std::size_t width, bool writing) {
    auto path = w.field(type, count, writing);
    auto values = [&](auto values) {
        if (writing) w.archive()[path] << values;
        else w.archive()[path] >> values;
        if (values.size()!=count) throw std::invalid_argument("Invalid Fortran checkpoint field size");
        return values;
    };
    if (type == 0) {
        auto* p = static_cast<char*>(data);
        std::vector<std::string> strings;
        for (std::size_t i=0; i<count; ++i) strings.push_back(writing ? text(p+i*width,width) : "");
        strings = values(std::move(strings));
        if (!writing) for (std::size_t i=0; i<count; ++i) put(p+i*width,width,strings[i]);
    } else numeric(type, [&](auto tag) {
        using T = decltype(tag);
        auto* p = static_cast<T*>(data);
        std::vector<T> v(count);
        if (writing) std::copy_n(p,count,v.begin());
        v = values(std::move(v));
        if (!writing) std::copy(v.begin(),v.end(),p);
    });
}
}
// No C++ exception crosses a Fortran frame. Callbacks return after alps_failed()
// and the driver rethrows the saved exception at the C++ boundary.
extern "C" {
void* alps_get_context(worker* w) noexcept { return w->context; }
void alps_set_context(worker* w, void* context) noexcept { w->context = context; }
std::int64_t alps_completed_sweeps(worker* w) noexcept { return w->completed_sweeps(); }
bool alps_failed(worker* w) noexcept { return w->failed(); }
void alps_fail(worker* w, char const* message) noexcept {
    w->guard([&] { throw std::invalid_argument(message); });
}
double alps_random(worker* w) noexcept {
    double value = 0; w->guard([&] { value = w->get_random()(); }); return value;
}
bool alps_parameter_defined(worker* w, char const* name) noexcept {
    bool value = false; w->guard([&] { value = w->get_parameters().exists(name); }); return value;
}
void alps_get_parameter(worker* w, void* data, char const* name, int type, std::size_t width) noexcept {
    w->guard([&] {
        auto value = w->get_parameters()[name];
        if (type==0) put(static_cast<char*>(data),width,value.as<std::string>());
        else numeric(type,[&](auto tag) { *static_cast<decltype(tag)*>(data) = value.as<decltype(tag)>(); });
    });
}
void alps_dump(worker* w, void* data, std::size_t count, int type, std::size_t width) noexcept {
    w->guard([&] { transfer(*w,data,count,type,width,true); });
}
void alps_restore(worker* w, void* data, std::size_t count, int type, std::size_t width) noexcept {
    w->guard([&] { transfer(*w,data,count,type,width,false); });
}
void alps_init_observable(worker* w, std::size_t count, char const* name) noexcept {
    w->guard([&] {
        if (!count || !*name || w->get_measurements().count(name) || std::string(name).find("Autocorrelation: ")==0)
            throw std::invalid_argument("Invalid or duplicate Fortran observable");
        alps::mc::add_measurement(*w,name,count,w->bins());
    });
}
void alps_accumulate_observable(worker* w, void* data, std::size_t count, int type, char const* name) noexcept {
    w->guard([&] {
        if (w->measurement(name)->size()!=count) throw std::invalid_argument("Fortran observable size mismatch");
        numeric(type,[&](auto tag) {
            auto* p = static_cast<decltype(tag)*>(data);
            std::vector<double> v(p,p+count);
            if (!std::all_of(v.begin(),v.end(),[](double x) { return std::isfinite(x); }))
                throw std::invalid_argument("Nonfinite Fortran observation");
            if (w->collecting()) alps::mc::record(*w,name,v);
        });
    });
}
}
