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

import ast
import urllib, copy, os, traceback
import numpy as np

import pyalps.hdf5 as h5
import pyalps.alea as pa

from .dataset import ResultFile
from .dataset import DataSet

import pyalps.pytools as pt # the C++ conversion functions


def log(m):
    print(m)
    
def parse_label(label):
    text = str(label)
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        # Released lattice labels encode a site pair as "source--target".
        return tuple(ast.literal_eval(part) for part in text.split('--'))


def parse_labels(labels):
    if isinstance(labels, (int, np.integer)):
        return np.array([labels])
    parsed = [parse_label(label) for label in labels]
    # Compress legacy site pairs with a common origin to their target sites.
    # A quoted string containing '--' is an ordinary label, not a site pair.
    if parsed and all(isinstance(value, tuple) and len(value) == 2 and
                      '--' in str(label) and value[0] == parsed[0][0]
                      for label, value in zip(labels, parsed)):
        parsed = [value[1] for value in parsed]
    return np.asarray(parsed)


class Hdf5Missing(Exception):
    def __init__(self,what):
        self.what = what
    
    def __str__(self):
        return 'Failed to find ' + self.what

class Hdf5Loader:
    """The Hdf5Loader class loads simulation parameters and observables from hdf5-files and returns them as hierarchical datasets"""
    def _read_replicas(self, filename, reader, suffix, measurements, verbose):
        replicas = sorted(self.h5f.list_children('/simulation/replicas'), key=int)
        self.h5f.close()
        result = []
        for replica in replicas:
            base = '/simulation/replicas/' + replica
            loaded = reader([filename], proppath=base + '/parameters', respath=base + suffix,
                            measurements=measurements, verbose=verbose)
            for group in loaded:
                for dataset in group:
                    dataset.props['replica'] = int(replica)
            result.extend(loaded)
        return result

    def GetFileNames(self, flist):
        files = []
        for f in flist:
          if f[-4:]=='.xml':
            f = f[:-3]+'h5'
          else:
            if f[-3:]!='.h5':
              f += '.h5'
          if os.path.exists(f):
            files.append(f)
          else:
            log( "FILE "+ f+ "DOES NOT EXIST!")
        return files
        
    def ReadParameters(self,proppath):
        dict = {'filename' : self.h5fname}
        if self.h5f.is_data(proppath + '/format') and str(self.h5f[proppath + '/format']).startswith('alps.params.'):
            from .ngs import params
            dict.update(params(self.h5f, proppath))
            return dict
        # Unmigrated applications still write legacy Parameters result groups.
        LOP=self.h5f.list_children(proppath)
        for m in LOP:
                try:
                    dict[m] = self.h5f[proppath+'/'+m]
                    try:
                        dict[m] = float(dict[m])
                    except:
                        dict[m] = list(map(float,dict[m]))
                except ValueError:
                    pass
        return dict 
        
    def GetProperties(self,flist,proppath='/parameters',respath='/simulation/results',verbose=False):
        fs = self.GetFileNames(flist)
        resultfiles = []
        for f in fs:
            try:
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log( "Loading from file" + f)
                rfile = ResultFile(f)
                rfile.props = self.ReadParameters(proppath)
                try:
                    obs = self.GetObservableList(respath)
                    rfile.props["ObservableList"] = [pt.hdf5_name_decode(x) for x in obs]
                except: pass
                resultfiles.append(rfile)
            except Exception as e:
                log(e)
                log(traceback.format_exc())
        return resultfiles
        
    def GetObservableList(self,respath):
        if self.h5f.is_group(respath):
            olist = self.h5f.list_children(respath)
        else:
            olist = []
        return olist

# Pre: file is a hdf5 file descriptor
# Post: returns DataSet with all parameters set
    
    def read_one_spectrum(self,path):
        pass
        
    def ReadSpectrumFromFile(self,flist,proppath='/parameters',respath='/spectrum',verbose=False):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            try:
                fileset=[]
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log("Loading from file " + f)
                params = self.ReadParameters(proppath)
                if 'energies' in self.h5f.list_children(respath):
                        try:
                            d = DataSet()
                            d.props['hdf5_path'] = respath 
                            d.props['observable'] = 'spectrum'
                            d.y = self.h5f[respath+'/energies']
                            d.x = range(len(d.y))
                            d.props.update(params)
                            try:
                                d.props.update(self.ReadParameters('quantumnumbers'))
                            except:
                                if verbose: log("no quantumnumbers stored ")
                                pass
                            fileset.append(d)
                        except AttributeError:
                            pass
                if 'sectors' in self.h5f.list_children(respath):
                    for secnum in self.h5f.list_children(respath+'/sectors'):
                        try:
                            d = DataSet()
                            secpath = respath+'/sectors/'+secnum
                            d.props['hdf5_path'] = secpath 
                            d.props['observable'] = 'spectrum'
                            d.y = self.h5f[secpath+'/energies']
                            d.x = range(len(d.y))
                            d.props.update(params)
                            try:
                                d.props.update(self.ReadParameters(secpath+'/quantumnumbers'))
                            except:
                                if verbose: log("no quantumnumbers stored ")
                                pass
                            fileset.append(d)
                        except AttributeError:
                            log( "Could not create DataSet")
                            pass
                sets.append(fileset)
            except Exception as e:
                log(e)
                log(traceback.format_exc())
        return sets
        
    def GetIterations(self, current_path, params={}, measurements=None, index=None, verbose=False):
        iterationset=[]
        #iteration_grp = self.h5f.require_group(respath+'/iteration')
        for it in self.h5f.list_children(current_path+'/iteration'):
            obsset=[]
            iteration_props = {}
            if 'parameters' in self.h5f.list_children(current_path+'/iteration/'+it):
                iteration_props = self.ReadParameters(current_path+'/iteration/'+it+'/parameters')
            iteration_props['iteration'] = it
            
            respath = current_path+'/iteration/'+it+'/results'
            list_ = self.GetObservableList(respath)
            if measurements is None:
                obslist = list_
            else:
                obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
            for m in obslist:
                if m in self.h5f.list_children(respath):
                    if "mean" in self.h5f.list_children(respath+'/'+m):
                        try:
                            d = DataSet()
                            itresultspath = respath+'/'+m
                            if verbose: log("Loading "+ m)
                            measurements_props = {}
                            measurements_props['hdf5_path'] = itresultspath 
                            measurements_props['observable'] = pt.hdf5_name_decode(m)
                            if index is None:
                                d.y = self.h5f[itresultspath+'/mean/value']
                                d.x = np.arange(0,len(d.y))
                            else:
                                try:
                                    d.y = self.h5f[itresultspath+'/mean/value'][index]
                                except:
                                    pass
                            if "labels" in self.h5f.list_children(itresultspath):
                                d.x = parse_labels(self.h5f[itresultspath+'/labels'])
                            else:
                                d.x = np.arange(0,len(d.y))
                            d.props.update(params)
                            d.props.update(iteration_props)
                            d.props.update(measurements_props)
                        except AttributeError:
                            log( "Could not create DataSet")
                    obsset.append(d)
            iterationset.append(obsset)
        return iterationset
        
    def ReadDiagDataFromFile(self,flist,proppath='/parameters',respath='/spectrum', measurements=None, index=None, loadIterations=False,verbose=False):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            try:
                fileset=[]
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log("Loading from file"+ f)
                params = self.ReadParameters(proppath)
                if 'results' in self.h5f.list_children(respath):
                    list_ = self.GetObservableList(respath+'/results')
                    if measurements is None:
                        obslist = list_
                    else:
                        obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                    if loadIterations==True:
                        if "iteration" in self.h5f.list_children(respath+'/results'):
                            fileset.append(self.GetIterations(respath+'/results', params, measurements, index, verbose))
                    else:        
                        for m in obslist:
                            if "mean" in self.h5f.list_children(respath+'/results/'+m):
                                try:
                                    if verbose: log("Loading" + m)
                                    d = DataSet()
                                    secresultspath = respath+'/results/'+m
                                    d.props['hdf5_path'] = secresultspath 
                                    d.props['observable'] = pt.hdf5_name_decode(m)
                                    if index is None:
                                        d.y = self.h5f[secresultspath+'/mean/value']
                                        d.x = np.arange(0,len(d.y))
                                    else:
                                        try:
                                            d.y = self.h5f[secresultspath+'/mean/value'][index]
                                        except:
                                            pass
                                    if "labels" in self.h5f.list_children(secresultspath):
                                        d.x = parse_labels(self.h5f[secresultspath+'/labels'])
                                    else:
                                        d.x = np.arange(0,len(d.y))
                                    d.props.update(params)
                                    
                                    fileset.append(d)
                                except AttributeError:
                                    log("Could not create DataSet")
                if loadIterations==True:
                    if "iteration" in self.h5f.list_children(respath):
                        fileset.append(self.GetIterations(respath, params, measurements, index, verbose))
                if 'sectors' in self.h5f.list_children(respath):
                    list_ = self.GetObservableList(respath+'/sectors/0/results')
                    if measurements is None:
                        obslist = list_
                    else:
                        obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                    for secnum in self.h5f.list_children(respath+'/sectors'):
                        sector_sets=[]
                        for m in obslist:
                            if "mean" in self.h5f.list_children(respath+'/sectors/'+secnum+'/results/'+m):
                                try:
                                    if verbose: log("Loading" + m)
                                    d = DataSet()
                                    secpath = respath+'/sectors/'+secnum
                                    secresultspath = respath+'/sectors/'+secnum+'/results/'+m
                                    d.props['hdf5_path'] = secresultspath 
                                    d.props['observable'] = pt.hdf5_name_decode(m)
                                    if index is None:
                                        d.y = self.h5f[secresultspath+'/mean/value']
                                        d.x = np.arange(0,len(d.y))
                                    else:
                                        try:
                                            d.y = self.h5f[secresultspath+'/mean/value'][index]
                                        except:
                                            pass
                                    if "labels" in self.h5f.list_children(secresultspath):
                                        d.x = parse_labels(self.h5f[secresultspath+'/labels'])
                                    else:
                                        d.x = np.arange(0,len(d.y))
                                    d.props.update(params)
                                    try:
                                        d.props.update(self.ReadParameters(secpath+'/quantumnumbers'))
                                    except:
                                        if verbose: log("no quantumnumbers stored ")
                                        pass
                                    sector_sets.append(d)

                                except AttributeError:
                                    log( "Could not create DataSet")
                                    pass
                        fileset.append(sector_sets)
                sets.append(fileset)
            except RuntimeError:
                raise
            except Exception as e:
                log(e)
                log(traceback.format_exc())
        return sets
        
    # Pre: file is a hdf5 file descriptor
    # Post: returns DataSet with the evaluated binning analysis set
    def ReadBinningAnalysis(self,flist,measurements=None,proppath='/parameters',respath=None,verbose=False):
        sets = []
        for filename in self.GetFileNames(flist):
            with h5.archive(filename, 'r') as archive:
                self.h5f, self.h5fname = archive, filename
                if verbose: log('Loading from file ' + filename)
                base = respath or '/simulation/results'
                if base == '/simulation/results' and archive.is_group('/simulation/replicas'):
                    sets.extend(self._read_replicas(filename, self.ReadBinningAnalysis,
                        '/realizations/0/clones/0/autocorrelation', measurements, verbose))
                    continue
                diagnostics = '/simulation/realizations/0/clones/0/autocorrelation'
                if base == '/simulation/results' and archive.is_group(diagnostics):
                    base = diagnostics
                names, params = self.GetObservableList(base), self.ReadParameters(proppath)
                selected = names if measurements is None else [pt.hdf5_name_encode(name)
                    for name in measurements if pt.hdf5_name_encode(name) in names]
                fileset = []
                for name in selected:
                    path = base+'/'+name
                    if archive.is_attribute(path+'/@kind') and archive[path+'/@kind'] == 4:
                        result = pa.read_result(archive, path)
                        levels = np.arange(result.levels)
                        errors = np.array([result.level(i).error for i in levels])
                        data = DataSet(levels, errors[:, 0] if result.size == 1 else errors)
                        data.native_result = result
                        data.props.update(params)
                        data.props.update(hdf5_path=path, observable='binning analysis of '+pt.hdf5_name_decode(name))
                        fileset.append(data)
                    elif archive.is_data(path+'/timeseries/logbinning'):
                        raise ValueError(f"{path}: legacy binning diagnostics require alps-hdf5-convert --alea-autocorr")
                sets.append(fileset)
        return sets

    # Pre: file is a hdf5 file descriptor
    # Post: returns DataSet with all parameters set
    def ReadMeasurementFromFile(self,flist,proppath='/parameters',respath='/simulation/results',measurements=None,verbose=False):
        sets = []
        for filename in self.GetFileNames(flist):
            with h5.archive(filename, 'r') as archive:
                self.h5f, self.h5fname = archive, filename
                if verbose: log("Loading from file " + filename)
                if respath == '/simulation/results' and archive.is_group('/simulation/replicas'):
                    sets.extend(self._read_replicas(filename, self.ReadMeasurementFromFile,
                        '/results', measurements, verbose))
                    continue
                names = self.GetObservableList(respath)
                params = self.ReadParameters(proppath)
                selected = names if measurements is None else [pt.hdf5_name_encode(name)
                    for name in measurements if pt.hdf5_name_encode(name) in names]
                fileset = []
                for name in selected:
                    path = respath+'/'+name
                    if verbose: log("Loading " + name)
                    if archive.is_data(path+'/histogram'):
                        values = archive[path+'/histogram']
                        x = archive[path+'/@min'] + archive[path+'/@stepsize']*np.arange(len(values))
                        data = DataSet(x, values)
                    elif archive.is_attribute(path+'/@kind') or archive.is_attribute(path+'/@format'):
                        result = pa.read_result(archive, path)
                        if not result.count:
                            continue
                        data = DataSet.from_result(result)
                    else:
                        children = archive.list_children(path)
                        if (any(field in children for field in ('count', 'timeseries', 'jacknife')) or
                                archive.is_data(path+'/mean/error')):
                            raise ValueError(f"{path}: legacy statistics require alps-hdf5-convert --alea-results "
                                             f"{path.rsplit('/', 1)[0] or '/'} (or --alea-batches/--alea-summary)")
                        # Deterministic eigenstate measurements carry plain values.
                        values = np.atleast_1d(archive[path+'/mean/value'])
                        data = DataSet(np.arange(len(values)), values)
                    if archive.is_data(path+'/labels'):
                        data.x = parse_labels(archive[path+'/labels'])
                    data.props.update(params)
                    data.props.update(hdf5_path=path, observable=pt.hdf5_name_decode(name))
                    fileset.append(data)
                sets.append(fileset)
        return sets

    # Pre: file is a hdf5 file descriptor
    # Post: returns DataSet with all parameters set
    def ReadDMFTIterations(self,flist,observable='G_tau',measurements='0',proppath='/parameters',respath='/simulation/iteration',verbose=False):
        fs = self.GetFileNames(flist)
        fileset = []
        for f in fs:
            try:
                self.h5f = h5.archive(f, 'r')
                self.h5fname = f
                if verbose: log("Loading from file "+ f)
                list_ = self.GetObservableList(respath+'/1/results/'+observable+'/')
                #grp = self.h5f.require_group(respath)
                params = self.ReadParameters(proppath)
                obslist = [pt.hdf5_name_encode(obs) for obs in measurements if pt.hdf5_name_encode(obs) in list_]
                iterationset=[]
                for it in  self.h5f.list_children(respath):
                    obsset=[]
                    for m in obslist:
                        try:
                            if verbose: log( "Loading "+ m)
                            d = DataSet()
                            size=0
                            path=it+'/results/'+observable+'/'+m
                            if "mean" in self.h5f.list_children(respath+'/'+path):
                                if self.h5f.is_scalar(respath+'/'+path+'/mean/value'):
                                    size=1
                                    obs=self.h5f[respath+'/'+path+'/mean/value']
                                    d.y = np.array([obs]) 
                                else:
                                    obs=self.h5f[respath+'/'+path+'/mean/value']
                                    size=len(obs)
                                    d.y = obs
                                d.x = np.arange(0,size)
                            d.props['hdf5_path'] = respath +"/"+ path
                            d.props['observable'] = pt.hdf5_name_decode(m)
                            d.props['iteration'] = it
                            d.props.update(params)
                        except AttributeError:
                            log( "Could not create DataSet")
                            pass
                        obsset.append(d)
                    iterationset.append(obsset)
                fileset.append(iterationset)
            except Exception as e:
                log( e)
                log( traceback.format_exc())
        return fileset
        
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
    return ll.ReadBinningAnalysis(files,measurements=what,verbose=verbose,respath=respath)

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
            #open the file and open the results root group
            h5file = h5.archive(f, 'r')
            #enumerate the subgroups
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
    

   
