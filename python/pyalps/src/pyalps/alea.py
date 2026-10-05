# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Olivier Parcollet
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

"""Native ALEA accumulators and results, plus remaining legacy analysis readers.

Use BatchAccumulator for bounded weighted histories, AutocorrelationAccumulator
for binning diagnostics, and NumPy arrays for full chronological sample history.
"""

from .cxx.pymcdata_c import *
from .cxx.pyalea_c import BatchAccumulator, BatchResult, ComplexBatchAccumulator, ComplexBatchResult
from .cxx.pyalea_c import (
    MeanAccumulator, MeanResult, ComplexMeanAccumulator, ComplexMeanResult,
    VarianceAccumulator, VarianceResult, ComplexVarianceAccumulator, ComplexVarianceResult,
    CovarianceAccumulator, CovarianceResult, ComplexCovarianceAccumulator, ComplexCovarianceResult,
    AutocorrelationAccumulator, AutocorrelationResult,
    ComplexAutocorrelationAccumulator, ComplexAutocorrelationResult,
    EllipticVarianceAccumulator, EllipticVarianceResult,
    EllipticCovarianceAccumulator, EllipticCovarianceResult,
    MeanTest, merge, ratio_real_imag,
)
import operator
import numpy
import pyalps.dataset


def _series(values, minimum=0):
    """Real observations, with time on axis zero; never infer samples from bins."""
    values = numpy.asarray(values)
    if values.dtype.kind not in "iuf" or values.ndim not in (1, 2):
        raise ValueError("expected a real sample array of shape (time,) or (time, components)")
    if len(values) < minimum or (values.ndim == 2 and values.shape[1] == 0):
        raise ValueError(f"expected at least {minimum} samples and nonempty components")
    if not numpy.isfinite(values).all():
        raise ValueError("samples must be finite")
    return values.astype(float, copy=False)


def size(timeseries):
    """Number of chronological observations (not number of vector components)."""
    return len(_series(timeseries))


def mean(timeseries):
    return _series(timeseries, 1).mean(axis=0)


def variance(timeseries):
    """Sample variance with N-1 normalization."""
    return _series(timeseries, 2).var(axis=0, ddof=1)


def make_dataset(timeseries):
    """Plot a sample array using one-based sample indices."""
    d = pyalps.dataset.DataSet()
    d.y = _series(timeseries)
    d.x = numpy.arange(1, len(d.y)+1)
    d.props['line'] = 'scatter'
    return d


def _selection(distance, limit):
    if (distance is None) == (limit is None):
        raise ValueError("specify exactly one of distance or limit")
    if limit is not None and (not numpy.isfinite(limit) or not 0 <= limit <= 1):
        raise ValueError("limit must be between zero and one")


def autocorrelation(timeseries, _distance=None, _limit=None):
    """Return positive-lag correlations, excluding lag zero.

    Each lag k uses (N-k) times the sample variance as its denominator.
    With _limit, include the first lag below that fraction of lag one's value
    in any component. Constant components have undefined correlation.
    """
    _selection(_distance, _limit)
    values = _series(timeseries, 2)
    n = len(values)
    distance = n-1 if _distance is None else operator.index(_distance)
    if distance < 0:
        distance += n
    if not 0 <= distance < n:
        raise ValueError("distance must be between zero and N-1")
    centered = values - values.mean(axis=0)
    var = centered.var(axis=0, ddof=1)
    if numpy.any(var == 0):
        raise ValueError("constant samples have undefined autocorrelation")
    # Zero padding prevents circular wraparound while avoiding quadratic work.
    spectrum = numpy.fft.rfft(centered, n=2*n, axis=0)
    sums = numpy.fft.irfft(spectrum * spectrum.conj(), n=2*n, axis=0)[1:distance+1]
    counts = numpy.arange(n-1, n-distance-1, -1).reshape((-1,) + (1,)*(values.ndim-1))
    correlation = sums / (counts * var)
    if _limit is not None:
        crossed = (correlation < _limit * correlation[0]).reshape(len(correlation), -1).any(axis=1)
        indices = numpy.flatnonzero(crossed)
        if len(indices):
            correlation = correlation[:indices[0]+1]
    return correlation


def _crossing(values, limit):
    if values.ndim != 1 or not len(values):
        raise ValueError("threshold selection requires a nonempty scalar series")
    indices = numpy.flatnonzero(values <= limit * values[0])
    return int(indices[0]) if len(indices) else len(values)


def _cut(timeseries, distance, limit, head):
    _selection(distance, limit)
    values = _series(timeseries)
    if distance is not None:
        count = operator.index(distance)
        if count < 0:
            count += len(values)
        if not 0 <= count <= len(values):
            raise ValueError("cut distance must be between zero and N")
        return values[count:] if head else values[:len(values)-count]
    stop = min(_crossing(values, limit)+1, len(values))
    return values[stop:] if head else values[:stop]


def cut_head(timeseries, _distance=None, _limit=None):
    """Remove initial samples, including the first threshold crossing if given."""
    return _cut(timeseries, _distance, _limit, True)


def cut_tail(timeseries, _distance=None, _limit=None):
    """Remove final samples, retaining the first threshold crossing if given."""
    return _cut(timeseries, _distance, _limit, False)


def exponential_autocorrelation_time(autocorrelation, _from=None, _to=None, _max=None, _min=None):
    """Fit A*exp(b*lag), returning (A, b), by least squares in log space.

    Bounds _from/_to are inclusive one-based lags; negative bounds add N.
    Alternatively select from
    the first _max crossing up to (excluding) the first _min crossing, where
    thresholds are fractions of the first correlation. At least two positive
    values are required. The amplitude always refers to lag zero.
    """
    values = _series(autocorrelation, 2)
    if values.ndim != 1:
        raise ValueError("exponential fitting requires a scalar series")
    if _from is not None and _to is not None and _max is None and _min is None:
        first, stop = operator.index(_from), operator.index(_to)
        first = first + len(values) if first < 0 else first
        stop = stop + len(values) if stop < 0 else stop
        start = first - 1
    elif _max is not None and _min is not None and _from is None and _to is None:
        if not 0 <= _min < _max <= 1:
            raise ValueError("fit limits require 0 <= min < max <= 1")
        start, stop = _crossing(values, _max), _crossing(values, _min)
    else:
        raise ValueError("specify either from/to or max/min")
    if not 0 <= start < stop <= len(values) or stop-start < 2:
        raise ValueError("fit range must contain at least two samples")
    if numpy.any(values[start:stop] <= 0):
        raise ValueError("exponential fit requires positive correlations")
    slope, intercept = numpy.polyfit(numpy.arange(start+1, stop+1), numpy.log(values[start:stop]), 1)
    return float(numpy.exp(intercept)), float(slope)


def integrated_autocorrelation_time(autocorrelation, fit=None):
    """Sum positive-lag correlations (without 1/2), optionally fitting the tail.

    The tail is the continuum integral of A*exp(b*t) from N+1/2 to infinity.
    The corresponding variance inflation factor is 1 + 2*tau.
    """
    values = _series(autocorrelation, 1)
    total = values.sum(axis=0)
    if fit is not None:
        amplitude, exponent = fit
        if values.ndim != 1 or not numpy.isfinite([amplitude, exponent]).all() or amplitude <= 0 or exponent >= 0:
            raise ValueError("tail fit requires a scalar series, positive amplitude and negative exponent")
        total -= amplitude / exponent * numpy.exp(exponent * (len(values)+0.5))
    return total


def running_mean(timeseries):
    values = _series(timeseries)
    counts = numpy.arange(1, len(values)+1).reshape((-1,) + (1,)*(values.ndim-1))
    return values.cumsum(axis=0) / counts


def reverse_running_mean(timeseries):
    return running_mean(_series(timeseries)[::-1])[::-1]


binning = "binning"
uncorrelated = "uncorrelated"


def error(timeseries, selector=uncorrelated):
    """Error of the mean; binning uses native ALEA's automatic level selection.

    Use AutocorrelationAccumulator directly to inspect convergence and levels.
    """
    values = _series(timeseries, 2)
    if selector == uncorrelated:
        return numpy.sqrt(values.var(axis=0, ddof=1) / len(values))
    if selector != binning:
        raise ValueError("error selector must be 'uncorrelated' or 'binning'")
    accumulator = AutocorrelationAccumulator(1 if values.ndim == 1 else values.shape[1])
    for sample in values:
        accumulator << sample
    result = accumulator.result().error
    return result[0] if values.ndim == 1 else result


class ReportedEstimate:
    """Published statistics with no inferred weights, covariance or history.

    Optional statistics are absent attributes when not reported. This record is
    not an accumulator and cannot be treated as a set of independent samples.
    """
    _format = 'alps.reported-estimate.v1'
    _fields = dict(mean='mean/value', error='mean/error', variance='variance/value',
                   tau='tau/value', converged_errors='mean/error_convergence')

    def __init__(self, *, count, mean, **statistics):
        if isinstance(count, (bool, numpy.bool_)) or not isinstance(count, (int, numpy.integer)) or not 0 <= count <= numpy.iinfo('u8').max:
            raise ValueError("reported count must be a uint64 integer")
        if statistics.keys() - (self._fields.keys() - {'mean'}):
            raise ValueError("unknown reported statistic")
        self.count = int(count)
        for name, value in dict(mean=mean, **statistics).items():
            value = numpy.asarray(value)
            kinds = 'iufc' if name == 'mean' else 'iu' if name == 'converged_errors' else 'iuf'
            if value.ndim != 1 or not value.size or value.dtype.kind not in kinds:
                raise ValueError(f"{name}: expected a nonempty numeric component vector")
            dtype = 'c16' if value.dtype.kind == 'c' else 'i8' if name == 'converged_errors' else 'f8'
            setattr(self, name, value.astype(dtype, copy=True))
        if any(getattr(self, name).shape != self.mean.shape for name in statistics):
            raise ValueError("reported statistics have different component counts")

    @property
    def size(self):
        return self.mean.size

    @classmethod
    def read(cls, archive, path):
        base = path.rstrip('/')+'/' if path else ''
        if archive[base+'@format'] != cls._format or archive.is_attribute(base+'@kind'):
            raise ValueError("expected a reported-estimate record")
        return cls(count=archive[base+'count'], **{
            name: archive[base+field] for name, field in cls._fields.items()
            if archive.is_data(base+field)})

    def save(self, archive, path=''):
        # Validate mutable arrays before changing a previously saved record.
        value = type(self)(**vars(self))
        base = path.rstrip('/')+'/' if path else ''
        if archive.is_attribute(base+'@kind'):
            raise ValueError("cannot overwrite a native estimator with a reported estimate")
        if archive.is_attribute(base+'@format') and archive[base+'@format'] != self._format:
            raise ValueError("cannot overwrite a different statistical format")
        archive.create_group(path or '.')
        archive[base+'@format'] = self._format
        archive[base+'count'] = numpy.uint64(value.count)
        for name, field in self._fields.items():
            if hasattr(value, name):
                archive[base+field] = getattr(value, name)
            elif archive.is_data(base+field):
                archive.delete_data(base+field)


def read_result(archive, path):
  """Read a native ALEA result or an explicitly marked reported estimate.

  ``archive`` is a pyalps.hdf5 archive (or a borrowed native archive). Released
  legacy encodings must first pass through the offline archive converter.
  Accumulator checkpoints are deliberately not interpreted as analysis results.
  """
  if archive.is_attribute(path + '/@format') and (
      not archive.is_attribute(path + '/@kind') or
      archive[path + '/@format'] == ReportedEstimate._format):
    return ReportedEstimate.read(archive, path)
  if not archive.is_attribute(path + '/@kind'):
    raise ValueError(f"{path}: missing native ALEA result metadata; use alps-hdf5-convert first")
  kind = archive[path + '/@kind']
  types = {1: (MeanResult, ComplexMeanResult),
           2: (VarianceResult, ComplexVarianceResult),
           3: (CovarianceResult, ComplexCovarianceResult),
           4: (AutocorrelationResult, ComplexAutocorrelationResult),
           5: (BatchResult, ComplexBatchResult)}
  if not isinstance(kind, (int, numpy.integer)) or isinstance(kind, (bool, numpy.bool_)) or kind not in types:
    raise ValueError(f"{path}: expected a native ALEA result, got kind {kind}")
  complex_values = archive.is_complex(path + '/mean/value')
  result_type = types[kind][complex_values]
  if complex_values and kind in (2, 3):
    moment_path = path + ('/var' if kind == 2 else '/cov')
    if len(archive.extent(moment_path)) == (3 if kind == 2 else 4):
      result_type = EllipticVarianceResult if kind == 2 else EllipticCovarianceResult
  return result_type.read(archive, path)
