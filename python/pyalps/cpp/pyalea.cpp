// Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
//               Lukas Gamper <gamperl@gmail.com>,
//               Matthias Troyer <troyer@itp.phys.ethz.ch>,
//               Maximilian Poprawe <poprawem@ethz.ch>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/alea/checkpoint.hpp>
#include <alps/alea.hpp>
#include <alps/alea/testing.hpp>
#include <nanobind/stl/function.h>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5.hpp>
#include "numpy_compat.hpp"
#include "archive_savable.hpp"
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
namespace nb = nanobind;
namespace {
// __repr__ helper for any ALPS type with an ostream operator.
template <typename T>
std::string stream_repr(T const & x) {
    std::ostringstream ss;
    ss << x;
    return ss.str();
}
template<class T> void bind_statistics_io(nb::class_<T>& cls) {
    cls.def("save", [](T const& value, nb::handle object, std::string const& path) {
        pyalps::with_native_archive(object, [&](auto& archive) {
            alps::alea::hdf5_serializer bridge(archive, path);
            serialize(bridge, "", value);
        });
    }, nb::arg("archive"), nb::arg("path") = "")
    .def("load", [](T& value, nb::handle object, std::string const& path) {
        pyalps::with_native_archive(object, [&](auto& archive) {
            alps::alea::hdf5_serializer bridge(archive, path);
            deserialize(bridge, "", value);
        });
    }, nb::arg("archive"), nb::arg("path") = "")
    .def_static("read", [](nb::handle object, std::string const& path) {
        T value;
        pyalps::with_native_archive(object, [&](auto& archive) {
            alps::alea::hdf5_serializer bridge(archive, path);
            deserialize(bridge, "", value);
        });
        return value;
    }, nb::arg("archive"), nb::arg("path") = "");
    pyalps::mark_archive_savable(cls);
}

namespace aa = alps::alea;

template<class Derived> nb::object estimate_array(Eigen::MatrixBase<Derived> const& values) {
    if constexpr (std::is_same_v<typename Derived::Scalar,aa::complex_op<double>>) {
        std::vector<double> data;
        for (Eigen::Index i=0;i<values.rows();++i) for (Eigen::Index j=0;j<values.cols();++j) {
            auto x=values(i,j);
            data.insert(data.end(),{x.rere(),x.reim(),x.imre(),x.imim()});
        }
        std::vector<size_t> shape{size_t(values.rows())};
        if constexpr (Derived::ColsAtCompileTime != 1) shape.push_back(values.cols());
        shape.insert(shape.end(),{2,2});
        return alps::python::make_numpy_array(std::move(data),shape);
    } else return nb::cast(values.eval());
}

template<class T> struct python_transform : aa::transformer<T> {
    std::function<aa::column<T>(aa::column<T> const&)> function;
    size_t input, output;
    python_transform(decltype(function) f, size_t in, size_t out)
        : function(std::move(f)), input(in), output(out) {
        if (!output) throw nb::value_error("A transform must return at least one component");
    }
    size_t in_size() const override { return input; }
    size_t out_size() const override { return output; }
    aa::column<T> operator()(aa::column<T> const& x) const override {
        auto y=function(x);
        if (size_t(y.size())!=output) throw nb::value_error("Transform output size changed");
        return y;
    }
};

template<class R> void bind_estimate(nb::class_<R>& result, nb::module_& module) {
    using T = typename aa::traits<R>::value_type;
    result.def_prop_ro("count", &R::count)
        .def_prop_ro("size", &R::size)
        .def("__copy__", [](R const& value) { return R(value); })
        .def("__deepcopy__", [](R const& value, nb::handle) { return R(value); })
        .def_prop_ro("mean", [](R const& value) { return value.mean().eval(); })
        .def("__repr__", &stream_repr<R>)
        .def("join", [](R const& a, R const& b) { return aa::join(a,b); }, nb::arg("other"),
             "Concatenate aligned batches, or assume independence for summary results.");
    module.def("merge", [](std::vector<R> const& values) { return aa::merge(values); }, nb::arg("results"));
    if constexpr (aa::traits<R>::HAVE_VAR) {
        result.def_prop_ro("error", [](R const& value) { return estimate_array(value.stderror()); })
            .def_prop_ro("variance", [](R const& value) { return estimate_array(value.var()); })
            .def_prop_ro("count2", &R::count2)
            .def_prop_ro("observations", &R::observations)
            .def("test_mean", [](R const& value, aa::column<T> const& expected, double atol) {
                return aa::test_mean(value,expected,atol);
            }, nb::arg("expected"), nb::arg("atol")=1e-14)
            .def("test_mean", [](R const& value, R const& expected, double atol) {
                return aa::test_mean(value,expected,atol);
            }, nb::arg("expected"), nb::arg("atol")=1e-14);
    }
    if constexpr (aa::traits<R>::HAVE_COV)
        result.def_prop_ro("covariance", [](R const& value) { return estimate_array(value.cov()); });
    if constexpr (aa::traits<R>::HAVE_TAU) {
        result.def_prop_ro("tau", [](R const& value) { return value.tau().eval(); })
            .def_prop_ro("tau_available", &R::tau_available)
            .def_prop_ro("converged_errors", &R::converged_errors)
            .def_prop_ro("levels", &R::nlevel)
            .def("level", [](R const& value, size_t i) {
                if (i>=value.nlevel()) throw nb::index_error();
                return value.level(i);
            }, nb::arg("index"));
    }
    if constexpr (!std::is_same_v<T,double>)
        result.def("real_components", [](R const& value) { return aa::real_components(value); });
    // Elliptic operators are real 2x2 covariance blocks, not complex scalars.
    // Keep them explicit; native ratio_real_imag handles signed estimators.
    if constexpr (!std::is_same_v<R,aa::var_result<std::complex<double>,aa::elliptic_var>>
                  && !std::is_same_v<R,aa::cov_result<std::complex<double>,aa::elliptic_var>>) {
        result.def("transform", [](R const& value,
                    std::function<aa::column<T>(aa::column<T> const&)> function,
                    size_t output_size, std::string const& method, double dx) -> nb::object {
            python_transform<T> f(std::move(function),value.size(),output_size);
            if (method=="none") return nb::cast(aa::transform(aa::no_prop(),f,value));
            if constexpr (aa::traits<R>::HAVE_VAR)
                if (method=="linear") return nb::cast(aa::transform(aa::linear_prop(dx),f,value));
            if constexpr (aa::traits<R>::HAVE_BATCH)
                if (method=="jackknife") return nb::cast(aa::transform(aa::jackknife_prop(),f,value));
            throw nb::value_error("Propagation method is unavailable for this estimator");
        }, nb::arg("function"),nb::arg("output_size")=1,
           nb::arg("method")=aa::traits<R>::HAVE_BATCH ? "jackknife" : aa::traits<R>::HAVE_VAR ? "linear" : "none",
           nb::arg("dx")=0.);
    }
    bind_statistics_io(result);
}

template<class A> void bind_accumulator(nb::class_<A>& accumulator) {
    using T = typename A::value_type;
    accumulator.def_prop_ro("size", &A::size)
        .def_prop_ro("count", &A::count)
        .def("__copy__", [](A const& value) { return A(value); })
        .def("__deepcopy__", [](A const& value, nb::handle) { return A(value); })
        .def("reset", &A::reset)
        .def("result", &A::result)
        .def("__lshift__", [](A& self, nb::handle sample) -> A& {
            auto array = alps::python::numpy_module().attr("asarray")(sample);
            auto kind = nb::cast<std::string>(array.attr("dtype").attr("kind"));
            if (std::string("biufc").find(kind) == std::string::npos
                    || (std::is_same_v<T, double> && kind == "c"))
                throw nb::type_error("ALEA expects a numeric sample of the accumulator's value type");
            auto values = alps::python::as_contiguous<T>(array);
            if (values.ndim() != 1) throw nb::value_error("ALEA samples must be scalars or one-dimensional arrays");
            typename aa::eigen<T>::const_col_map vector(values.data(), values.shape(0));
            self << aa::make_adapter(vector);
            return self;
        }, nb::rv_policy::none);
    bind_statistics_io(accumulator);
}

template<class A> void bind_moments(nb::module_& module, std::string const& prefix) {
    using R = typename aa::traits<A>::result_type;
    nb::class_<R> result(module,(prefix+"Result").c_str());
    bind_estimate(result,module);
    nb::class_<A> accumulator(module,(prefix+"Accumulator").c_str());
    accumulator.def(nb::init<size_t,uint64_t>(),nb::arg("size")=1,nb::arg("batch_size")=1);
    bind_accumulator(accumulator);
}

template<class T> void bind_native_statistics(nb::module_& module, std::string const& prefix) {
    using M = aa::mean_acc<T>; using MR = aa::mean_result<T>;
    nb::class_<MR> mean(module,(prefix+"MeanResult").c_str()); bind_estimate(mean,module);
    nb::class_<M> ma(module,(prefix+"MeanAccumulator").c_str());
    ma.def(nb::init<size_t>(),nb::arg("size")=1); bind_accumulator(ma);
    bind_moments<aa::var_acc<T>>(module,prefix+"Variance");
    bind_moments<aa::cov_acc<T>>(module,prefix+"Covariance");
    using C = aa::autocorr_acc<T>; using CR = aa::autocorr_result<T>;
    nb::class_<CR> corr(module,(prefix+"AutocorrelationResult").c_str()); bind_estimate(corr,module);
    nb::class_<C> ca(module,(prefix+"AutocorrelationAccumulator").c_str());
    ca.def(nb::init<size_t,uint64_t,size_t>(),nb::arg("size")=1,nb::arg("batch_size")=1,nb::arg("granularity")=2);
    bind_accumulator(ca);
    using A = aa::batch_acc<T>; using R = aa::batch_result<T>;
    nb::class_<R> result(module,(prefix+"BatchResult").c_str()); bind_estimate(result,module);
    result.def_prop_ro("batch_sums", [](R const& value) { return value.store().batch().transpose().eval(); })
        .def_prop_ro("batch_counts", [](R const& value) { return value.store().count(); });
    nb::class_<A> accumulator(module,(prefix+"BatchAccumulator").c_str());
    accumulator.def(nb::init<size_t,size_t,uint64_t>(),nb::arg("size")=1,nb::arg("num_batches")=64,nb::arg("base_size")=1)
        .def_prop_ro("batch_offsets", [](A const& value) { return value.offset().eval(); });
    bind_accumulator(accumulator);
}
} // namespace
NB_MODULE(pyalea_c, m) {
    m.doc() = "ALPS alea bindings (nanobind)";
    bind_native_statistics<double>(m, "");
    bind_native_statistics<std::complex<double>>(m, "Complex");
    bind_moments<aa::var_acc<std::complex<double>,aa::elliptic_var>>(m,"EllipticVariance");
    bind_moments<aa::cov_acc<std::complex<double>,aa::elliptic_var>>(m,"EllipticCovariance");
    m.def("ratio_real_imag", &aa::ratio_real_imag, nb::arg("result"));
    nb::class_<aa::t2_result>(m,"MeanTest")
        .def_prop_ro("score", &aa::t2_result::score)
        .def_prop_ro("pvalue", &aa::t2_result::pvalue)
        .def_prop_ro("pvalue_lower", &aa::t2_result::pvalue_lower)
        .def_prop_ro("pvalue_upper", &aa::t2_result::pvalue_upper)
        .def_prop_ro("has_plower", &aa::t2_result::has_plower);
}
