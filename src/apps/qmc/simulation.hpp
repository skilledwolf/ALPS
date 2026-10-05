// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include "../mc/measurements.hpp"
#include <set>
#include <valarray>

namespace native_qmc {
// Scientific measurements carry their own aligned sign denominator. This also
// handles improved estimators, whose sign can differ from a configuration sign.
class simulation : public alps::mcbase {
public:
  simulation(alps::params const& p,size_t bins,size_t chain):alps::mcbase(p,chain),bins_(bins) {}
  static alps::params checkpoint_parameters(alps::params p) { p.erase("SWEEPS"); return p; }
  alps::params sampling_parameters() const { return parameters; }
  void record_measurements(bool enabled) { recording_=enabled; }
  auto const& measurement_labels() const { return labels_; }
  auto const& signed_measurements() const { return signed_names_; }
  void add_measurement(std::string const& name,size_t size=1,bool sign=true) {
    if (is_signed_ && sign) signed_names_.insert(name);
    native_mc::add_measurement(*this,name,size+signed_names_.count(name),bins_);
  }
  void add_measurement(std::string const& name,std::vector<std::string> const& labels,bool sign=true) {
    if (labels.empty()) throw std::invalid_argument("Empty QMC measurement labels: "+name);
    labels_[name]=labels;
    add_measurement(name,labels.size(),sign);
  }
  template<class T> void record(std::string const& name,T const& value,double sign) {
    if (!recording_) return;
    auto adapter=alps::alea::make_adapter(value);
    alps::alea::column<double> sample(adapter.size()+signed_names_.count(name));
    sample.setZero();
    adapter.add_to(alps::alea::view<double>(sample.data(),adapter.size()));
    if (signed_names_.count(name)) sample[sample.size()-1]=sign;
    native_mc::record(*this,name,sample);
  }
  void record(std::string const& name,std::valarray<double> const& value,double sign) {
    record(name,std::vector<double>(std::begin(value),std::end(value)),sign);
  }
  void validate_measurements(alps::hdf5::archive& ar,uint64_t samples) const {
    if (ar.list_children("measurements").size()!=measurements.size())
      throw std::invalid_argument("Unexpected QMC checkpoint measurements");
    alps::alea::hdf5_serializer codec(ar,"measurements");
    for (auto const& [name,handle]:measurements) {
      if (!std::holds_alternative<std::shared_ptr<native_mc::batch>>(handle)) continue;
      native_mc::batch batch;
      native_mc::autocorr diagnostic;
      alps::alea::deserialize(codec,ar.encode_segment(name),batch);
      alps::alea::deserialize(codec,ar.encode_segment(native_mc::diagnostic(name)),diagnostic);
      if (batch.size()!=std::get<std::shared_ptr<native_mc::batch>>(handle)->size() ||
          batch.num_batches()!=bins_ || batch.current_batch_size()!=batch.cursor().factor() ||
          batch.count()>samples || !batch.store().batch().allFinite() ||
          diagnostic.size()!=batch.size() || diagnostic.count()!=batch.count() ||
          diagnostic.batch_size()!=1 || diagnostic.granularity()!=2)
        throw std::invalid_argument("Invalid QMC checkpoint measurement state");
    }
  }
protected:
  size_t bins_;
  bool is_signed_=false;
  std::map<std::string,std::vector<std::string>> labels_;
  std::set<std::string> signed_names_;
  bool recording_=true;
};
} // namespace native_qmc
