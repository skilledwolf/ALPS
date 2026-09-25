# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Bela Bauer <bauerb@phys.ethz.ch>
#                            Brigitte Surer <surerb@phys.ethz.ch> 
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import urllib, copy, os, traceback
import numpy as np

import pyalps.hdf5 as h5
import pyalps.alea as pa

from .dataset import ResultFile
from .dataset import DataSet
#from floatwitherror import FloatWithError as fwe

import pyalps.pytools as pt # the C++ conversion functions

# or the C++ class as alternative
from pyalps.alea import MCScalarData as fwe
from pyalps.alea import MCVectorData as vwe

def log(m):
    print(m)
    
def parse_label(label):
    if '--' in label:
      vals = label.rsplit('--')
      ret = ()
      for val in vals:
          ret = ret + (eval(val),)
      return ret
    else:
      return eval(str(label))
 
 
def parse_labels(labels):
    if type(labels)==int: 
      return np.array([labels])
    larr=[]
    allsame = True
    first = None
    for x in labels:
      v = parse_label(x)
      larr.append(v)
      if '--' in x:      
        if first is None:
          first = v[0]
        else:
          if first != v[0] or len(v) != 2:
            allsame = False
      else:
        allsame = False
    if allsame:
      larr = [x[1] for x in larr]
    return np.array(larr)
       
class Hdf5Missing(Exception):
    def __init__(self,what):
        self.what = what
    
    def __str__(self):
        return 'Failed to find ' + self.what

class Hdf5Loader:
    """Load HDF5 simulation results into hierarchical datasets."""

    def GetFileNames(self, flist):
        files = []
        for filename in flist:
            if filename.endswith('.xml'):
                filename = filename[:-3] + 'h5'
            elif not filename.endswith('.h5'):
                filename += '.h5'
            if os.path.exists(filename):
                files.append(filename)
            else:
                log("FILE " + filename + "DOES NOT EXIST!")
        return files

    def _read_files(self, flist, proppath, verbose, read_file, propagate=()):
        """Own archive lifetime and the common per-file error policy."""
        results = []
        for filename in self.GetFileNames(flist):
            try:
                with h5.archive(filename, 'r') as self.h5f:
                    self.h5fname = filename
                    if verbose:
                        log("Loading from file " + filename)
                    result = read_file(self.ReadParameters(proppath))
                results.append(result)
            except Exception as error:
                if isinstance(error, propagate):
                    raise
                log(error)
                log(traceback.format_exc())
        return results

    def ReadParameters(self, proppath):
        params = {'filename': self.h5fname}
        for name in self.h5f.list_children(proppath):
            try:
                params[name] = self.h5f[proppath + '/' + name]
                try:
                    params[name] = float(params[name])
                except Exception:
                    params[name] = list(map(float, params[name]))
            except ValueError:
                pass
        return params

    def GetObservableList(self, respath):
        return self.h5f.list_children(respath) if self.h5f.is_group(respath) else []

    def _observables(self, respath, measurements):
        available = self.GetObservableList(respath)
        if measurements is None:
            return available
        names = map(pt.hdf5_name_encode, measurements)
        return [name for name in names if name in available]

    @staticmethod
    def _dataset(path, observable, values, params, x=None):
        return DataSet(
            x=np.arange(len(values)) if x is None else x,
            y=values,
            props={'hdf5_path': path, 'observable': observable, **params},
        )

    def _axis(self, path, values, size=None):
        if 'labels' in self.h5f.list_children(path):
            return parse_labels(self.h5f[path + '/labels'])
        return np.arange(len(values) if size is None else size)

    def _quantumnumbers(self, path, verbose):
        try:
            return self.ReadParameters(path)
        except Exception:
            if verbose:
                log("no quantumnumbers stored ")
            return {}

    def GetProperties(self, flist, proppath='/parameters', respath='/simulation/results', verbose=False):
        def read(params):
            result = ResultFile(self.h5fname)
            result.props = params
            try:
                result.props['ObservableList'] = [
                    pt.hdf5_name_decode(name) for name in self.GetObservableList(respath)
                ]
            except Exception:
                pass
            return result
        return self._read_files(flist, proppath, verbose, read)

    def ReadSpectrumFromFile(self, flist, proppath='/parameters', respath='/spectrum', verbose=False):
        def read(params):
            datasets = []
            children = self.h5f.list_children(respath)
            paths = [(respath, 'quantumnumbers')] if 'energies' in children else []
            if 'sectors' in children:
                paths.extend((respath + '/sectors/' + sector,
                              respath + '/sectors/' + sector + '/quantumnumbers')
                             for sector in self.h5f.list_children(respath + '/sectors'))
            for path, quantum_path in paths:
                try:
                    values = self.h5f[path + '/energies']
                    dataset = self._dataset(path, 'spectrum', values, params, range(len(values)))
                    dataset.props.update(self._quantumnumbers(quantum_path, verbose))
                    datasets.append(dataset)
                except AttributeError:
                    if path != respath:
                        log("Could not create DataSet")
            return datasets
        return self._read_files(flist, proppath, verbose, read)

    def _diagonal_observables(self, respath, params, names, index, verbose):
        datasets = []
        for name in names:
            path = respath + '/' + name
            if 'mean' not in self.h5f.list_children(path):
                continue
            try:
                if verbose:
                    log("Loading " + name)
                values = np.array([])
                if index is None:
                    values = self.h5f[path + '/mean/value']
                else:
                    try:
                        values = self.h5f[path + '/mean/value'][index]
                    except Exception:
                        pass
                datasets.append(self._dataset(
                    path, pt.hdf5_name_decode(name), values, params,
                    self._axis(path, values),
                ))
            except AttributeError:
                log("Could not create DataSet")
        return datasets

    def GetIterations(self, current_path, params=None, measurements=None, index=None, verbose=False):
        iterations = []
        for iteration in self.h5f.list_children(current_path + '/iteration'):
            path = current_path + '/iteration/' + iteration
            props = self.ReadParameters(path + '/parameters') if 'parameters' in self.h5f.list_children(path) else {}
            props['iteration'] = iteration
            respath = path + '/results'
            datasets = self._diagonal_observables(
                respath, {}, self._observables(respath, measurements), index, verbose,
            )
            for dataset in datasets:
                # Iteration metadata overrides global parameters; measurement
                # identity overrides both, as in the original iteration reader.
                dataset.props = {**(params or {}), **props, **dataset.props}
            iterations.append(datasets)
        return iterations

    def ReadDiagDataFromFile(self, flist, proppath='/parameters', respath='/spectrum', measurements=None, index=None, loadIterations=False, verbose=False):
        def read(params):
            datasets = []
            children = self.h5f.list_children(respath)
            if 'results' in children:
                path = respath + '/results'
                if loadIterations:
                    if 'iteration' in self.h5f.list_children(path):
                        datasets.append(self.GetIterations(path, params, measurements, index, verbose))
                else:
                    datasets.extend(self._diagonal_observables(
                        path, params, self._observables(path, measurements), index, verbose,
                    ))
            if loadIterations and 'iteration' in children:
                datasets.append(self.GetIterations(respath, params, measurements, index, verbose))
            if 'sectors' in children:
                names = self._observables(respath + '/sectors/0/results', measurements)
                for sector in self.h5f.list_children(respath + '/sectors'):
                    path = respath + '/sectors/' + sector
                    sector_sets = self._diagonal_observables(path + '/results', params, names, index, verbose)
                    for dataset in sector_sets:
                        dataset.props.update(self._quantumnumbers(path + '/quantumnumbers', verbose))
                    datasets.append(sector_sets)
            return datasets
        return self._read_files(flist, proppath, verbose, read, propagate=(RuntimeError,))

    def ReadBinningAnalysis(self, flist, measurements=None, proppath='/parameters', respath=None, verbose=False):
        respath = '/simulation/results' if respath is None else respath
        def read(params):
            datasets = []
            for name in self._observables(respath, measurements):
                try:
                    path = respath + '/' + name
                    if 'timeseries' not in self.h5f.list_children(path):
                        continue
                    series = path + '/timeseries'
                    required = {'logbinning', 'logbinning2', 'logbinning_counts'}
                    if not required.issubset(self.h5f.list_children(series)):
                        continue
                    if verbose:
                        log("Loading " + name)
                    bins = self.h5f[series + '/logbinning'][:-7]
                    bins2 = self.h5f[series + '/logbinning2'][:-7]
                    counts = self.h5f[series + '/logbinning_counts'][:-7]
                    scale = 1
                    for i in range(len(counts)):
                        mean = bins[i] / (counts[i] * scale)
                        mean2 = bins2[i] / counts[i]
                        bins2[i] = np.sqrt((mean2 - mean * mean) / counts[i])
                        scale *= 2
                    datasets.append(self._dataset(
                        respath + name, 'binning analysis of ' + pt.hdf5_name_decode(name), bins2, params,
                    ))
                    if verbose:
                        log("  loaded binning analysis for " + name)
                except AttributeError:
                    log("Could not create DataSet")
            return datasets
        return self._read_files(flist, proppath, verbose, read)

    def _measurement(self, path, params):
        """Decode one MC observable while retaining its scalar/vector form."""
        if 'histogram' in self.h5f.list_children(path):
            values = self.h5f[path + '/histogram']
            xmin = self.h5f[path + '/@min']
            step = self.h5f[path + '/@stepsize']
            return values, np.arange(xmin, xmin + step * len(values), step)
        scalar = self.h5f.is_scalar(path + '/mean/value')
        if 'error' in self.h5f.list_children(path + '/mean'):
            if scalar:
                observable = pa.MCScalarData()
                observable.load(self.h5fname, path)
                values, size = np.array([observable]), 1
                if observable.count == 0:
                    values = None
            elif not self.h5f.is_group(path + '/timeseries'):
                values = np.array(self.h5f[path + '/mean/value'])
                if 'L' in params:
                    length = int(params['L'])
                    if length == values.size:
                        params['origin'] = [(length - 1.) / 2.]
                    if length ** 2 == values.size:
                        values = values.reshape([length, length])
                        params['origin'] = [(length - 1.) / 2.] * 2
                    elif length ** 3 == values.size:
                        values = values.reshape([length, length, length])
                        params['origin'] = [(length - 1.) / 2.] * 3
                size = values.size
            else:
                values = pa.MCVectorData()
                values.load(self.h5fname, path)
                size = len(values.mean)
                if values.count == 0:
                    values = None
        else:
            values = self.h5f[path + '/mean/value']
            if scalar:
                values = np.array([values])
            size = len(values)
        return values, self._axis(path, values, size)

    def ReadMeasurementFromFile(self, flist, proppath='/parameters', respath='/simulation/results', measurements=None, verbose=False):
        def read(params):
            datasets = []
            for name in self._observables(respath, measurements):
                if verbose:
                    log("Loading " + name)
                path = respath + '/' + name
                values, axis = self._measurement(path, params)
                try:
                    if values is not None:
                        datasets.append(self._dataset(path, pt.hdf5_name_decode(name), values, params, axis))
                except AttributeError:
                    log("Could not create DataSet")
            return datasets
        return self._read_files(flist, proppath, verbose, read)

    def ReadDMFTIterations(self, flist, observable='G_tau', measurements='0', proppath='/parameters', respath='/simulation/iteration', verbose=False):
        def read(params):
            names = self._observables(respath + '/1/results/' + observable + '/', measurements)
            iterations = []
            for iteration in self.h5f.list_children(respath):
                datasets = []
                for name in names:
                    try:
                        if verbose:
                            log("Loading " + name)
                        path = respath + '/' + iteration + '/results/' + observable + '/' + name
                        values = np.array([])
                        if 'mean' in self.h5f.list_children(path):
                            values = self.h5f[path + '/mean/value']
                            if self.h5f.is_scalar(path + '/mean/value'):
                                values = np.array([values])
                        datasets.append(self._dataset(
                            path, pt.hdf5_name_decode(name), values, {'iteration': iteration, **params},
                        ))
                    except AttributeError:
                        log("Could not create DataSet")
                iterations.append(datasets)
            return iterations
        return self._read_files(flist, proppath, verbose, read)


def loadBinningAnalysis(files,what=None,verbose=False,respath='/simulation/results'):
    """ loads MC binning analysis from ALPS HDF5 result files
    
        this function loads results of a MC binning analysis from ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            what (list): optional argument that is either a string or list of strings, specifying the names of the observables for which the binning analysis should be loaded
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            a list of list of DataSet objects: loaded binning analysis.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the inner list are each for a different observable.
            The x-values of the DataSet objects are the logarithmic binning level and the y-values the error estimates at that binning level.
    """
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadBinningAnalysis(files,measurements=what,verbose=verbose)

def loadMeasurements(files,what=None,verbose=False,respath='/simulation/results'):
    """ loads ALPS measurements from ALPS HDF5 result files
    
        this function loads results of ALPS simulations ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            what (list): optional argument that is either a string or list of strings, specifying the names of the observables which should be loaded
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            a list of list of DataSet objects: loaded measurements.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the inner list are each for a different observable.
            The y-values of the DataSet objects are the measurements and the x-values optionally the labels (indices) of array-valued measurements
    """
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadMeasurementFromFile(files,measurements=what,verbose=verbose,respath=respath)
    
    
def loadEigenstateMeasurements(files,what=None, verbose=False):
    """ loads ALPS eigenstate measurements from ALPS HDF5 result files
    
        this function loads results of ALPS diagonalization or DMRG simulations from an HDF5 file
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            what (list): an optional argument that is either a string or list of strings, specifying the names of the observables which should be loaded
            verbose (bool): an optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            list of list of (lists of) DataSet objects: loaded measurements.
            The elements of the outer list each correspond to the file names specified as input
            The elements of the next level are different quantum number sectors, if any exists
            The elements of the inner-most list are each for a different observable
            The y-values of the DataSet objects is an array of the measurements in all eigenstates calculated in this sector, and the x-values optionally the labels (indices) of array-valued measurements
    """
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadDiagDataFromFile(files,measurements=what,verbose=verbose)
    
def loadIterationMeasurements(files,what=None,verbose=False):
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadDiagDataFromFile(files,measurements=what,loadIterations=True,verbose=verbose)

def loadSpectra(files,verbose=False):
    """ loads ALPS spectra from ALPS HDF5 result files
    
        This function loads the spectra calculated in ALPS diagonalization or DMRG simulations from an HDF5 file.
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded.
        
        Returns:
            list of (lists of) DataSet objects: Loaded spectra.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the next level are different quantum number sectors, if any exists.
            The y-values of the DataSet objects are the energies in that quantum number sector.
    """
    ll = Hdf5Loader()
    return ll.ReadSpectrumFromFile(files,verbose=verbose)

def loadDMFTIterations(files,observable='G_tau',measurements='0',verbose=False):
    """ loads ALPS measurements from ALPS HDF5 result files
    
        this function loads results of ALPS simulations ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS HDF5 result files.
            observable (str): optional argument specifying the name of the observables which should be loaded
            measurements (list): optional argument that is either a string or list of strings, specifying the names of the measurements which should be loaded
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            list of list of list of DataSet objects: loaded iteration measurements.
            The elements of the outer list each correspond to the file names specified as input.
            The elements of the next level are different iterations.
            The elements of the inner list contains a DataSet for each measurement.
            The y-values of the DataSet objects are the measurements and the x-values optionally the labels (indices) of array-valued measurements
    """
    ll = Hdf5Loader()
    if isinstance(measurements,str):
      measurements = [measurements]
    return ll.ReadDMFTIterations(files,observable=observable,measurements=measurements,verbose=verbose)

def loadProperties(files,proppath='/parameters',respath='/simulation/results',verbose=False):
    """ loads properties (parameters) of simulations from ALPS HDF5 result files
    
        this function loads the properties (parameters) of ALPS simulations ALPS HDF5 result files
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
        
        Returns:
            list of dicts: properties contained in each file.
    """
    ll = Hdf5Loader()
    res = ll.GetProperties(files,proppath,respath,verbose=verbose)
    results = []
    for x in res:
      results.append(x.props)
    return results

def loadObservableList(files,proppath='/parameters',respath='/simulation/results',verbose=False):   
    """ loads lists of existing measurements from ALPS HDF5 result files
    
        The function returns a list of lists, containing the names of measurements that are stored in the result files 
        
        Parameters:
            files (list): ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
            verbose (bool): optional argument that if set to True causes more output to be printed as the data is loaded
    """
    ll = Hdf5Loader()
    res = ll.GetProperties(files,proppath,respath,verbose=verbose)
    results = []
    for x in res:
      results.append(x.props['ObservableList'])
    return results

def loadTimeEvolution( flist,globalproppath='/parameters',resroot='/timesteps/',localpropsuffix='/parameters', measurements=None):
    ll=Hdf5Loader()
    data=[]
    #loop over files
    for f in flist:
        try:
            with h5.archive(f, 'r') as h5file:
                L=h5file.list_children(resroot)
            #Create an iterator of length the number of subgroups
            stepper=[i+1 for i in range(len(L))]
            #Read in global props
            globalprops=ll.GetProperties([f],globalproppath)
            for d in stepper:
                #Get the measurements from the numbered subgroups
                locdata=ll.ReadMeasurementFromFile([f],proppath=resroot+str(d)+localpropsuffix, \
                respath=resroot+str(d)+'/results', measurements=measurements)
                #Append the global props to the local props
                for i in range(len(locdata[0])):
                    locdata[0][i].props.update(globalprops[0].props)
                #Extend the total dataset with this data
                data.extend(locdata)
        except Exception as e:
            log(e)
            log( traceback.format_exc())
    return data
    

   
