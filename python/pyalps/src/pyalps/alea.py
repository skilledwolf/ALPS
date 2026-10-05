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
from .cxx.pyalea_c import MCScalarTimeseries, MCScalarTimeseriesView, MCVectorTimeseries, MCVectorTimeseriesView, ValueWithError, StdPairDouble, size, mean, variance, integrated_autocorrelation_time, running_mean, reverse_running_mean
from . import alea_detail as detail
import numpy
import pyalps.dataset


# TODO: __repr__ function for StdPairDouble
#def std_pair_double.__repr__():
#  return str(self.first) + ", " + str(self.second)

def make_dataset(MCTimeseries):
  "Makes a pyalps.dataset from a MCTimeseries object making it possible to plot one using pyalps.plot"
  d = pyalps.dataset.DataSet()
  d.y = MCTimeseries.timeseries()
  d.x = numpy.arange(1,len(d.y)+1,1)
  d.props['line'] = 'scatter'
  return d

def autocorrelation(timeseries, _distance = None, _limit = None):
  "Calculates the autocorrelation of a given timeseries. Can be invoked with _distance or _limit.\n\
timeseries: Any MCTimeseries or MCData object\n\
_distance: Calculates the autocorrelation until a specific length\n\
_limit: Calculates the autocorrelation until it has reached _limit of its initial value\n\
returns: MCTimeseries object with the autocorrelation"
  if _distance is not None:
    return detail.autocorrelation_distance(timeseries, _distance)
  if _limit is not None:
    return detail.autocorrelation_limit(timeseries, _limit)
  print("Usage: autocorrelation(timeseries, [_distance = XXX | _limit = XXX] )")

def cut_head(timeseries, _distance = None, _limit = None):
  "Creates a MCTimeseriesView object. Can be invoked with _distance or _limit.\n\
timeseries: Any MCTimeseries object\n\
_distance: Cuts the first _distance entries\n\
_limit: Cuts the front until the timeseries reaches _limit of its initial value\n\
returns: MCTimeseriesView object with the smaller timeseries\n\
  Note: does not copy the data, only creates a reference."
  if _distance is not None:
    return detail.cut_head_distance(timeseries, _distance)
  if _limit is not None:
    return detail.cut_head_limit(timeseries, _limit)
  print("Usage: cut_head(timeseries, [_distance = XXX | _limit = XXX] )")

def cut_tail(timeseries, _distance = None, _limit = None):
  "Creates a MCTimeseriesView object. Can be invoked with _distance or _limit.\n\
timeseries: Any MCTimeseries object\n\
_distance: Cuts the last _distance entries\n\
_limit: Cuts the tail until the timeseries only decays from its initial value to _limit of its initial value\n\
returns: MCTimeseriesView object with the smaller timeseries\n\
  Note: does not copy the data, only creates a reference."
  if _distance is not None:
    return detail.cut_tail_distance(timeseries, _distance)
  if _limit is not None:
    return detail.cut_tail_limit(timeseries, _limit)
  print("Usage: cut_head(timeseries, [_distance = XXX | _limit = XXX] )")

def exponential_autocorrelation_time(autocorrelation, _from = None, _to = None, _max = None, _min = None):
  "Fits a timeseries exponentially. Can be invoked with _from and _to or with _max and _min.\n\
autocorrelation: A MCTimeseries object\n\
_from & _to: fits the autocorrelation between _from and _to\n\
_max & _min: fits the autocorrelation between the values where it is at _max and where it is at _min from its initial value\n\
returns: StdPairDouble object FIT with the parameters of the fit\n\
  Note: The equation is FIT.first * exp(FIT.second * t)"
  if (_from is not None and _to is not None):
    return detail.exponential_autocorrelation_time_distance(autocorrelation, _from, _to)
  if (_max is not None and _min is not None):
    return detail.exponential_autocorrelation_time_limit(autocorrelation, _max, _min)
  print("Usage: exponential_autocorrelation_time(autocorrelation, [_from = XXX, _to = XXX | _max = XXX, _min = XXX] )")

binning = "binning"
uncorrelated = "uncorrelated"

def error(timeseries, selector = uncorrelated):
  "Calculates an error estimate of a timeseries. Can be invoked with 'uncorrelated' or 'binning'\n\
timeseries: Any MCTimeseries or MCData object\n\
uncorrelated (default): The error estimate when uncorrelated data is assumed\n\
binning: The error estimate is calculated using a binning analysis.\n\
returns: A float or numpyarray (depending on the dimension of the timeseries) with the error(s)"
  if selector == "binning":
    return detail.binning_error(timeseries)
  if selector == "uncorrelated":
    return detail.uncorrelated_error(timeseries)


def read_result(archive, path):
  """Read a typed native ALEA result without discarding batches or covariance.

  ``archive`` is a pyalps.hdf5 archive (or a borrowed native archive). Released
  legacy encodings must first pass through the offline archive converter.
  Accumulator checkpoints are deliberately not interpreted as analysis results.
  """
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
