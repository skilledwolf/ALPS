// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/alea/checkpoint.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/alea/convert.hpp>
#include <alps/alea/transform.hpp>
#include <map>

namespace alps::mc {
using moment_results=std::vector<alps::alea::cov_result<double>>;
// Reuse the batch accumulator's cursor and native centered covariance algebra.
// Moments describe individual physical samples within each bin, not bin means.
class physical_moments {
    std::vector<alps::alea::cov_acc<double>> bins_;
public:
    physical_moments(size_t size,size_t bins):bins_(bins,alps::alea::cov_acc<double>(size)){}
    void add(alps::alea::batch_acc<double> const& layout, alps::alea::column<double> const& sample) {
        auto next=layout.cursor();
        if (layout.store().count()(next.current())>=layout.current_batch_size()) {
            ++next;
            if (next.merge_mode()) {
                bins_[next.merge_into()] << bins_[next.current()].result();
                bins_[next.current()].reset();
            }
        }
        bins_[next.current()] << alps::alea::make_adapter(sample);
    }
    moment_results results() const {
        moment_results out;
        for (auto const& bin:bins_) out.push_back(bin.result());
        return out;
    }
    void save(alps::hdf5::archive& ar,std::string const& path) const {
        alps::alea::hdf5_serializer codec(ar,path);
        for (size_t i=0;i<bins_.size();++i) alps::alea::serialize(codec,std::to_string(i),bins_[i]);
    }
    void load(alps::hdf5::archive& ar,std::string const& path,alps::alea::batch_acc<double> const& layout) {
        if (ar.list_children(path).size()!=bins_.size()) throw std::invalid_argument("Wrong physical moment bin count");
        alps::alea::hdf5_serializer codec(ar,path);
        auto staged=bins_;
        for (size_t i=0;i<staged.size();++i) {
            auto& bin=staged[i];
            alps::alea::deserialize(codec,std::to_string(i),bin);
            if (bin.size()!=bins_[i].size() || bin.batch_size()!=1 || bin.current().count()!=0
                    || bin.count()!=layout.store().count()(i) || bin.store().count2()!=double(bin.count())
                    || !bin.store().data().allFinite() || !bin.store().data2().allFinite())
                throw std::invalid_argument("Invalid physical moment state");
        }
        bins_=std::move(staged);
    }
};

// Express each bin around one pooled reference before native jackknife analysis.
// The output keeps the original weights and all within-bin physical covariance.
inline auto centered_batches(moment_results const& bins) {
    auto pooled=alps::alea::merge(bins);
    auto n=pooled.size();
    alps::alea::column<double> reference=alps::alea::column<double>::Zero(n);
    if (pooled.count()) reference=pooled.mean();
    alps::alea::batch_data<double> data(n+n*n,bins.size());
    for (size_t i=0;i<bins.size();++i) {
        auto count=bins[i].count();
        data.count()(i)=count;
        if (!count) continue;
        alps::alea::column<double> delta=bins[i].mean()-reference;
        Eigen::MatrixXd centered=count*(delta*delta.transpose());
        if (count>1) centered+=(count-1)*bins[i].cov();
        data.batch().col(i).head(n)=count*delta;
        data.batch().col(i).tail(n*n)=Eigen::Map<Eigen::VectorXd>(centered.data(),n*n);
    }
    return std::make_pair(alps::alea::batch_result<double>(data),reference);
}

// One common publication rule for every optional thermodynamic estimate.
using unavailable_results=std::map<std::string,std::string>;
template<class F> void estimate(std::map<std::string,alps::alea::batch_result<double>>& results,
        unavailable_results* unavailable, std::string const& name,
        alps::alea::batch_result<double> const& input, F evaluate) {
    std::string reason;
    if ((input.store().count().array()>0).count()<2) reason="At least two occupied batches are required";
    else {
        struct function : alps::alea::transformer<double> {
            size_t size; F evaluate;
            function(size_t n,F f):size(n),evaluate(std::move(f)){}
            size_t in_size()const override{return size;} size_t out_size()const override{return 1;}
            alps::alea::column<double> operator()(alps::alea::column<double> const& x)const override{return {evaluate(x)};}
        } transform(input.size(),std::move(evaluate));
        auto value=alps::alea::transform(alps::alea::jackknife_prop{},transform,input);
        if (value.store().batch().allFinite() && value.mean().allFinite() && value.stderror().allFinite()) {
            results.emplace(name,std::move(value));
            return;
        }
        reason="Value or uncertainty is undefined or exceeds floating-point range";
    }
    if (unavailable) (*unavailable)[name]=std::move(reason);
}
}
