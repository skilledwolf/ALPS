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
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/alea/mcdata.hpp>
#include <alps/alea/mcanalyze.hpp>
#include <alps/alea/value_with_error.hpp>
#include <alps/alea/checkpoint.hpp>
#include <alps/alea.hpp>
#include <alps/alea/testing.hpp>
#include <nanobind/stl/function.h>
#include <alps/alea/hdf5.hpp>
#include <alps/hdf5.hpp>
#include <alps/numeric/vector_functions.hpp>
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
// Build a 1-D numpy.ndarray (dtype=float64) from any sequence-like
// alps container (std::vector<double>, std::valarray<double>).
template <typename Container>
nb::object seq_to_numpy(Container const & v) {
    std::size_t n = static_cast<std::size_t>(v.size());
    std::vector<double> tmp(n);
    for (std::size_t i = 0; i < n; ++i)
        tmp[i] = static_cast<double>(v[i]);
    return alps::python::make_numpy_array<double>(std::move(tmp), {n});
}
// Copy a numpy array into a std::vector<ValueType>. Used when the
// caller still instantiates mctimeseries<ValueType> from a
// Python array.
template <typename ValueType>
std::vector<ValueType>
numpy_to_vector(nb::handle arr) {
    auto view = alps::python::as_contiguous<double>(arr);
    if (view.ndim() != 1)
        throw std::invalid_argument(
            "mctimeseries ctor: expected 1-D array");
    std::size_t n = static_cast<std::size_t>(view.shape(0));
    std::vector<ValueType> out(n);
    double const * data = view.data();
    for (std::size_t i = 0; i < n; ++i)
        out[i] = static_cast<ValueType>(data[i]);
    return out;
}
// __repr__ helper for any ALPS type with an ostream operator.
template <typename T>
std::string stream_repr(T const & x) {
    std::ostringstream ss;
    ss << x;
    return ss.str();
}
template <typename T>
std::string value_with_error_repr(alps::alea::value_with_error<T> const & v) {
    std::ostringstream ss;
    ss << v.mean() << " +/- " << v.error();
    return ss.str();
}
// Numpy-returning wrappers for vector-valued alps::alea free
// functions (mean, variance, uncorrelated_error, binning_error) —
// the scalar overloads are bound directly and nanobind casts
// their `double` return to a Python float automatically.
template <typename T>
nb::object mean_vector(T const & x) {
    return seq_to_numpy(alps::alea::mean(x));
}
template <typename T>
nb::object variance_vector(T const & x) {
    return seq_to_numpy(alps::alea::variance(x));
}
template <typename T>
nb::object uncorrelated_error_vector(T const & x) {
    return seq_to_numpy(alps::alea::uncorrelated_error(x));
}
template <typename T>
nb::object binning_error_vector(T const & x) {
    return seq_to_numpy(alps::alea::binning_error(x));
}
// mctimeseries.timeseries() returns std::vector<ValueType>; hand
// back to Python as numpy. For scalar ValueType we pack 1-D; for
// vector<double> ValueType we pack 2-D. mctimeseries_view has the
// same surface.
template <typename ValueType, typename TS>
nb::object ts_to_numpy_scalar(TS const & ts) {
    auto const & v = ts.timeseries();
    return seq_to_numpy(v);
}
template <typename TS>
nb::object ts_to_numpy_vector_rows(TS const & ts) {
    auto const & rows = ts.timeseries();
    if (rows.empty())
        return alps::python::make_numpy_array<double>(
            std::vector<double>(), {std::size_t{0}, std::size_t{0}});
    std::size_t nrows = rows.size();
    std::size_t ncols = rows.front().size();
    for (auto const & row : rows)
        if (row.size() != ncols)
            throw std::runtime_error("mctimeseries has ragged rows; cannot shape as numpy 2-D");
    std::vector<double> flat(nrows * ncols);
    double * dst = flat.data();
    for (auto const & row : rows) {
        for (std::size_t j = 0; j < ncols; ++j)
            *dst++ = static_cast<double>(row[j]);
    }
    return alps::python::make_numpy_array<double>(std::move(flat), {nrows, ncols});
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
    // ─── value_with_error ────────────────────────────────────────────
    nb::class_<alps::alea::value_with_error<double>>(m, "ValueWithError")
        .def(nb::init<double, double>(),
             nb::arg("mean") = 0.0, nb::arg("error") = 0.0)
        .def_prop_ro("mean",  &alps::alea::value_with_error<double>::mean)
        .def_prop_ro("error", &alps::alea::value_with_error<double>::error)
        .def("__repr__", &value_with_error_repr<double>);
    // ─── StdPairDouble ─────────────────────────────────────────────
    // nanobind's STL caster already registered std::pair<double, double> as a
    // Python tuple converter via the stl.h header; nanobind takes the
    // same path via <nanobind/stl/pair.h>. Binding it again as a
    // class_ would fight that, so use a thin attribute-access
    // wrapper instead, and give integrated_autocorrelation_time a
    // Python-side signature that accepts either StdPairDouble or a
    // plain (float, float) tuple.
    struct StdPairDouble {
        double first{0.0};
        double second{0.0};
        StdPairDouble() = default;
        StdPairDouble(double f, double s) : first(f), second(s) {}
        operator std::pair<double, double>() const { return {first, second}; }
    };
    nb::class_<StdPairDouble>(m, "StdPairDouble",
        "Pair of (fit slope, fit intercept) returned by the autocorrelation fit helpers.")
        .def(nb::init<>())
        .def(nb::init<double, double>(), nb::arg("first"), nb::arg("second"))
        .def_rw("first",  &StdPairDouble::first)
        .def_rw("second", &StdPairDouble::second)
        .def("__repr__", [](StdPairDouble const & p) {
            std::ostringstream ss;
            ss << "StdPairDouble(" << p.first << ", " << p.second << ")";
            return ss.str();
        });
    // ─── mctimeseries<T> / mctimeseries_view<T> bindings ────────────
    //
    // Numpy-ctor: take any handle and copy through numpy_to_vector.
    // The timeseries() method returns numpy arrays directly; the
    // 2-D overload for vector<double> goes through the
    // ts_to_numpy_vector_rows helper.
    #define ALPS_PY_EXPORT_MCTIMESERIES_SCALAR(Value, PyName)                             \
        nb::class_<alps::alea::mctimeseries<Value>>(m, PyName)                            \
            .def(nb::init<>())                                                            \
            .def(nb::init<alps::alea::mcdata<Value>>())                                   \
            .def("__init__",                                                              \
                 [](alps::alea::mctimeseries<Value> * self, nb::handle a) {               \
                     new (self) alps::alea::mctimeseries<Value>(                          \
                         numpy_to_vector<Value>(a));                                      \
                 })                                                                       \
            .def("timeseries", [](alps::alea::mctimeseries<Value> const & self) {         \
                    return ts_to_numpy_scalar<Value>(self);                               \
                })                                                                        \
            .def_prop_ro("size", &alps::alea::mctimeseries<Value>::size)                  \
            .def("__repr__", &stream_repr<alps::alea::mctimeseries<Value>>);              \
        nb::class_<alps::alea::mctimeseries_view<Value>>(m, PyName "View")                \
            .def(nb::init<alps::alea::mctimeseries<Value>>())                             \
            .def(nb::init<alps::alea::mctimeseries_view<Value>>())                        \
            .def("timeseries", [](alps::alea::mctimeseries_view<Value> const & self) {    \
                    return ts_to_numpy_scalar<Value>(self);                               \
                })                                                                        \
            .def_prop_ro("size", &alps::alea::mctimeseries_view<Value>::size)             \
            .def("__repr__", &stream_repr<alps::alea::mctimeseries_view<Value>>)
    ALPS_PY_EXPORT_MCTIMESERIES_SCALAR(double, "MCScalarTimeseries");
    #undef ALPS_PY_EXPORT_MCTIMESERIES_SCALAR
    // Vector-valued mctimeseries: ctor from a 2-D numpy array, rows
    // are time samples. timeseries() returns 2-D.
    using VecTs    = alps::alea::mctimeseries<std::vector<double>>;
    using VecTsV   = alps::alea::mctimeseries_view<std::vector<double>>;
    using VecMcD   = alps::alea::mcdata<std::vector<double>>;
    nb::class_<VecTs>(m, "MCVectorTimeseries")
        .def(nb::init<>())
        // Register typed overloads before the catch-all array constructor.
        .def(nb::init<VecMcD>())
        .def("__init__", [](VecTs * self, nb::handle a) {
                auto view = alps::python::as_contiguous<double>(a);
                if (view.ndim() != 2)
                    throw std::invalid_argument(
                        "MCVectorTimeseries ctor: expected 2-D array");
                std::size_t nrows = static_cast<std::size_t>(view.shape(0));
                std::size_t ncols = static_cast<std::size_t>(view.shape(1));
                std::vector<std::vector<double>> rows(nrows);
                double const * data = view.data();
                for (std::size_t i = 0; i < nrows; ++i) {
                    rows[i].resize(ncols);
                    for (std::size_t j = 0; j < ncols; ++j)
                        rows[i][j] = data[i * ncols + j];
                }
                new (self) VecTs(rows);
            })
        .def("timeseries", [](VecTs const & self) { return ts_to_numpy_vector_rows(self); })
        .def_prop_ro("size", &VecTs::size)
        .def("__repr__", &stream_repr<VecTs>);
    nb::class_<VecTsV>(m, "MCVectorTimeseriesView")
        .def(nb::init<VecTs>())
        .def(nb::init<VecTsV>())
        .def("timeseries", [](VecTsV const & self) { return ts_to_numpy_vector_rows(self); })
        .def_prop_ro("size", &VecTsV::size)
        .def("__repr__", &stream_repr<VecTsV>);
    // ─── alps::alea free functions over mcdata / mctimeseries ────────
    #define DEF_ALL(name, fn)                                                             \
        /* scalar-valued */                                                               \
        m.def(name, static_cast<double (*)(alps::alea::mcdata<double> const &)>(&fn));    \
        m.def(name, static_cast<double (*)(alps::alea::mctimeseries<double> const &)>(&fn)); \
        m.def(name, static_cast<double (*)(alps::alea::mctimeseries_view<double> const &)>(&fn));
    // size — works for both scalar and vector value types.
    m.def("size", static_cast<std::size_t (*)(alps::alea::mcdata<double> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries<double> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries_view<double> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mcdata<std::vector<double>> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries<std::vector<double>> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries_view<std::vector<double>> const &)>(&alps::size));
    // mean — scalar overloads return double, vector overloads return numpy.
    DEF_ALL("mean", alps::alea::mean)
    m.def("mean", &mean_vector<alps::alea::mcdata<std::vector<double>>>);
    m.def("mean", &mean_vector<alps::alea::mctimeseries<std::vector<double>>>);
    m.def("mean", &mean_vector<alps::alea::mctimeseries_view<std::vector<double>>>);
    // variance — same pattern.
    DEF_ALL("variance", alps::alea::variance)
    m.def("variance", &variance_vector<alps::alea::mcdata<std::vector<double>>>);
    m.def("variance", &variance_vector<alps::alea::mctimeseries<std::vector<double>>>);
    m.def("variance", &variance_vector<alps::alea::mctimeseries_view<std::vector<double>>>);
    // integrated_autocorrelation_time — scalar only. The C++ signature
    // takes the (slope, intercept) pair by const-ref.
    //
    // NOTE: four overloads in total — the two std::pair forms here
    // (satisfied by any 2-tuple via <nanobind/stl/pair.h>) and the two
    // StdPairDouble forms further down. They are disjoint today because
    // StdPairDouble's implicit conversion to std::pair is invisible to
    // nanobind; keep the pair overloads registered FIRST and do not add
    // an implicitly_convertible between the two, or the dispatch order
    // silently changes.
    m.def("integrated_autocorrelation_time",
          static_cast<double (*)(alps::alea::mctimeseries<double> const &,
                                 std::pair<double, double> const &)>(
              &alps::alea::integrated_autocorrelation_time));
    m.def("integrated_autocorrelation_time",
          static_cast<double (*)(alps::alea::mctimeseries_view<double> const &,
                                 std::pair<double, double> const &)>(
              &alps::alea::integrated_autocorrelation_time));
    // running_mean / reverse_running_mean — scalar only, mctimeseries-valued.
    m.def("running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mcdata<double> const &)>(
              &alps::alea::running_mean));
    m.def("running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries<double> const &)>(
              &alps::alea::running_mean));
    m.def("running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries_view<double> const &)>(
              &alps::alea::running_mean));
    m.def("reverse_running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mcdata<double> const &)>(
              &alps::alea::reverse_running_mean));
    m.def("reverse_running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries<double> const &)>(
              &alps::alea::reverse_running_mean));
    m.def("reverse_running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries_view<double> const &)>(
              &alps::alea::reverse_running_mean));
    #undef DEF_ALL
    // ─── mcanalyze free functions consumed by pyalps.alea ────────────
    //
    // autocorrelation / cut_head / cut_tail / error in alea.py dispatch
    // to these through alea_detail; the instantiation set matches the
    // Boost.Python module.
    #define DEF_TS_SCALAR(name, fn)                                                      \
        m.def(name, &fn<alps::alea::mcdata<double>>);                                    \
        m.def(name, &fn<alps::alea::mctimeseries<double>>);                              \
        m.def(name, &fn<alps::alea::mctimeseries_view<double>>);
    #define DEF_TS_VECTOR(name, fn)                                                      \
        m.def(name, &fn<alps::alea::mcdata<std::vector<double>>>);                       \
        m.def(name, &fn<alps::alea::mctimeseries<std::vector<double>>>);                 \
        m.def(name, &fn<alps::alea::mctimeseries_view<std::vector<double>>>);
    // autocorrelation — by distance or by decay limit.
    DEF_TS_SCALAR("autocorrelation_distance", alps::alea::autocorrelation_distance)
    DEF_TS_VECTOR("autocorrelation_distance", alps::alea::autocorrelation_distance)
    DEF_TS_SCALAR("autocorrelation_limit", alps::alea::autocorrelation_limit)
    DEF_TS_VECTOR("autocorrelation_limit", alps::alea::autocorrelation_limit)
    // head/tail cuts — views by distance for scalar and vector series;
    // by decay limit for scalar series only.
    DEF_TS_SCALAR("cut_head_distance", alps::alea::cut_head_distance)
    DEF_TS_VECTOR("cut_head_distance", alps::alea::cut_head_distance)
    DEF_TS_SCALAR("cut_tail_distance", alps::alea::cut_tail_distance)
    DEF_TS_VECTOR("cut_tail_distance", alps::alea::cut_tail_distance)
    DEF_TS_SCALAR("cut_head_limit", alps::alea::cut_head_limit)
    DEF_TS_SCALAR("cut_tail_limit", alps::alea::cut_tail_limit)
    // error estimates — scalar overloads return float, vector overloads numpy.
    DEF_TS_SCALAR("uncorrelated_error", alps::alea::uncorrelated_error)
    DEF_TS_SCALAR("binning_error", alps::alea::binning_error)
    m.def("uncorrelated_error", &uncorrelated_error_vector<alps::alea::mcdata<std::vector<double>>>);
    m.def("uncorrelated_error", &uncorrelated_error_vector<alps::alea::mctimeseries<std::vector<double>>>);
    m.def("uncorrelated_error", &uncorrelated_error_vector<alps::alea::mctimeseries_view<std::vector<double>>>);
    m.def("binning_error", &binning_error_vector<alps::alea::mcdata<std::vector<double>>>);
    m.def("binning_error", &binning_error_vector<alps::alea::mctimeseries<std::vector<double>>>);
    m.def("binning_error", &binning_error_vector<alps::alea::mctimeseries_view<std::vector<double>>>);
    #undef DEF_TS_SCALAR
    #undef DEF_TS_VECTOR
    // exponential_autocorrelation_time fits — return StdPairDouble so the
    // fit.first / fit.second attribute API documented in pyalps.alea
    // is preserved.
    m.def("exponential_autocorrelation_time_distance",
          [](alps::alea::mctimeseries<double> const & ts, int from, int to) {
              std::pair<double, double> fit =
                  alps::alea::exponential_autocorrelation_time_distance(ts, from, to);
              return StdPairDouble(fit.first, fit.second);
          });
    m.def("exponential_autocorrelation_time_distance",
          [](alps::alea::mctimeseries_view<double> const & ts, int from, int to) {
              std::pair<double, double> fit =
                  alps::alea::exponential_autocorrelation_time_distance(ts, from, to);
              return StdPairDouble(fit.first, fit.second);
          });
    m.def("exponential_autocorrelation_time_limit",
          [](alps::alea::mctimeseries<double> const & ts, double max, double min) {
              std::pair<double, double> fit =
                  alps::alea::exponential_autocorrelation_time_limit(ts, max, min);
              return StdPairDouble(fit.first, fit.second);
          });
    m.def("exponential_autocorrelation_time_limit",
          [](alps::alea::mctimeseries_view<double> const & ts, double max, double min) {
              std::pair<double, double> fit =
                  alps::alea::exponential_autocorrelation_time_limit(ts, max, min);
              return StdPairDouble(fit.first, fit.second);
          });
    // integrated_autocorrelation_time also accepts the StdPairDouble
    // returned by the fit helpers, in addition to a plain 2-tuple.
    m.def("integrated_autocorrelation_time",
          [](alps::alea::mctimeseries<double> const & ts, StdPairDouble const & fit) {
              return alps::alea::integrated_autocorrelation_time(
                  ts, std::pair<double, double>(fit.first, fit.second));
          });
    m.def("integrated_autocorrelation_time",
          [](alps::alea::mctimeseries_view<double> const & ts, StdPairDouble const & fit) {
              return alps::alea::integrated_autocorrelation_time(
                  ts, std::pair<double, double>(fit.first, fit.second));
          });
}
