// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/mc/measurements.hpp>
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
    alps::mc::add_measurement(*this,name,size+signed_names_.count(name),bins_);
  }
  void add_measurement(std::string const& name,std::vector<std::string> const& labels,bool sign=true) {
    if (labels.empty()) throw std::invalid_argument("Empty QMC measurement labels: "+name);
    labels_[name]=labels;
    add_measurement(name,labels.size(),sign);
  }
  using samples_type=std::vector<std::pair<std::string,std::vector<double>>>;
  // Capture raw, signed observations for routing between replica ranks.
  // Replay uses the same native accumulators in the original sample order.
  template<class F> samples_type sample(F const& advance) {
    samples_type samples;
    pending_=&samples;
    try {advance();} catch (...) {pending_=nullptr;throw;}
    pending_=nullptr;
    return samples;
  }
  void record(samples_type const& samples) {
    for (auto const& [name,value]:samples) alps::mc::record(*this,name,value);
  }
  template<class T> void record(std::string const& name,T const& value,double sign) {
    if (!recording_) return;
    auto adapter=alps::alea::make_adapter(value);
    alps::alea::column<double> sample(adapter.size()+signed_names_.count(name));
    sample.setZero();
    adapter.add_to(alps::alea::view<double>(sample.data(),adapter.size()));
    if (signed_names_.count(name)) sample[sample.size()-1]=sign;
    if (pending_) pending_->emplace_back(name,std::vector<double>(sample.data(),sample.data()+sample.size()));
    else alps::mc::record(*this,name,sample);
  }
  void record(std::string const& name,std::valarray<double> const& value,double sign) {
    record(name,std::vector<double>(std::begin(value),std::end(value)),sign);
  }
  void validate_measurements(alps::hdf5::archive& ar,uint64_t samples) const {
    if (ar.list_children("measurements").size()!=measurements.size())
      throw std::invalid_argument("Unexpected QMC checkpoint measurements");
    alps::alea::hdf5_serializer codec(ar,"measurements");
    for (auto const& [name,handle]:measurements) {
      if (!std::holds_alternative<std::shared_ptr<alps::mc::batch>>(handle)) continue;
      alps::mc::batch batch;
      alps::mc::autocorr diagnostic;
      alps::alea::deserialize(codec,ar.encode_segment(name),batch);
      alps::alea::deserialize(codec,ar.encode_segment(alps::mc::diagnostic(name)),diagnostic);
      if (batch.size()!=std::get<std::shared_ptr<alps::mc::batch>>(handle)->size() ||
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
  samples_type* pending_=nullptr;
};
} // namespace native_qmc
