# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2012 by Sebastian Keller
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import pyalps.load
import shutil
import tempfile
import socket
import subprocess
import sys
import os
import os.path
import numpy as np
from xml.etree import ElementTree
from glob import glob
from time import gmtime, strftime

#*******************************************************************************
# Misc
#*******************************************************************************

def read_testprop_xml(testinput):
    """ Returns a dictionary with the
    testinfo contents of the test xml file, 
    a list of reference files and a list of observables"""

    root = ElementTree.parse(testinput).getroot()
    props = {}
    props['SCRIPT'] = root.find('SCRIPT').attrib['file']
    infos = root.find('TESTINFOS')
    for tinfo in list(infos):
        if tinfo.tag == 'TESTINFO':
            props[tinfo.attrib['name']] = tinfo.text

    files = root.find('REFERENCE')    
    refFileList = [ rf.attrib['file'] for rf in list(files) ]
    if not refFileList or None in refFileList:
        print("Could not create list of Reference Files")
        sys.exit(1)

    obsentries = root.find('OBSERVABLES')
    obslist = []
    if obsentries is not None:
        obslist = [ obs.text for obs in list(obsentries) ]

    try:
        inputs = [ inpf.attrib['file'] for inpf in list(root.find('SCRIPTINPUTFILES')) ]

    except AttributeError:
        inputs = None

    return props, inputs, refFileList, obslist

def tasknr( fname ):
    """ Extract the task number from alps filenames """
    nstart = fname.rfind('task')+4
    nend   = fname[nstart:].find('.')+nstart
    try: return int( fname[nstart:nend] )
    except ValueError: return ''

def getPrefix( fname ):
    fstring = os.path.basename( fname )
    ind = fstring.rfind( '.task' )
    return fstring[:ind]

def recursive_mkdir( path ):
    """ Create directory path including parents if necessary """
    dirs = []; path = os.path.abspath( path )
    while not os.path.exists( path ):
        dirs.append( path )
        path = path[:path.rindex('/')]
    for d in reversed( dirs ):
        os.mkdir( d )

#*******************************************************************************
# Output
#*******************************************************************************

def writeTest2stdout(comparelist, fobj):
    """ Write test results to stdout in XML Format

    compare_list[task][obs] is a nested list of dicts, which hold the
    test results """
    
    sys.stdout = fobj
    for comptask in comparelist:

        obs_list = []
        for obs in comptask:
            if obs['passed']:
                continue
            else:
                obs_list.append( obs )

        if not obs_list:
            print("<TESTTASK file='%s' status='PASS'>" % os.path.basename(comptask[0]['filename']))
            print("</TESTTASK>")
        else:
            print("<TESTTASK file='%s' status='FAIL'>" % comptask[0]['filename'])
            print("\t<OBSERVABLES>")
            for obs in obs_list:
                print("\t<OBSERVABLE name='%s' diff='%s' tol='%s'/>" % \
                    (obs['observable'], obs['difference'], obs['tolerance']))

            print("\t</OBSERVABLES>")
            print("</TESTTASK>")

    sys.stdout = sys.__stdout__


#*******************************************************************************
# Comparison
#*******************************************************************************
def obsdict(tol, diff, props):
    return {'difference':diff, 'tolerance':tol, \
            'observable': props['observable'], \
            'filename': props['filename'], \
            'passed'  : bool(tol - diff >= 0) }
            #TODO: Include sector information for epsilon-type measurements

def _mc_difference(test, ref, vector):
    tolerance = np.sqrt(test.error**2 + ref.error**2)
    difference = np.abs(test.mean - ref.mean)
    return tolerance, difference


def _epsilon_difference(test, ref, vector):
    return max(10e-12, (ref if vector else np.abs(ref)) * 10e-12), np.abs(test - ref)


def _compare_values(test, ref, tol_factor, difference, *, scalar_index=False):
    if pyalps.size(test) == 1:
        if scalar_index:
            test, ref = test[0], ref[0]
        tol, diff = difference(test, ref, False)
    else:
        pairs = [difference(ty, ry, True) for ty, ry in zip(test, ref)]
        differences = [pair[1] for pair in pairs]
        diff = max(differences)
        tol = pairs[differences.index(diff)][0]
        # Preserve the historical MC vector tolerance (factor squared).
        if difference is _mc_difference:
            tol *= tol_factor
    return tol * tol_factor, diff


def _compare_mc(test, ref, factor):
    return _compare_values(test.y[0], ref.y[0], factor, _mc_difference)


def _compare_epsilon(test, ref, factor):
    return _compare_values(test.y[0], ref.y[0], factor, _epsilon_difference)


def _compare_mixed(test, ref, factor):
    try:
        return _compare_values(test.y, ref.y, factor, _mc_difference, scalar_index=True)
    except AttributeError:
        return _compare_values(test.y, ref.y, factor, _epsilon_difference, scalar_index=True)


def _select_observables(observables, whatlist, filename, source):
    if not whatlist:
        return observables
    names = [obs.props['observable'] for obs in observables]
    missing = [name for name in whatlist if name not in names]
    if missing:
        print("The following observables specified for comparison\nhave not been found in %s results:" % source)
        print("File:", filename)
        print(missing)
        sys.exit(1)
    return [obs for obs in observables if obs.props['observable'] in whatlist]


def _compare_data(testdata, refdata, factor, whatlist, compare, *, sectors=False):
    if not sectors and len(testdata) != len(refdata):
        raise Exception("Comparison Error: test and reference data differ in number of tasks")
    comparisons = []
    for testtask, reftask in zip(testdata, refdata):
        if sectors:
            try:
                testfile = testtask[0][0].props['filename']
                reffile = reftask[0][0].props['filename']
            except AttributeError:
                testtask, reftask = [testtask], [reftask]
                testfile = testtask[0][0].props['filename']
                reffile = reftask[0][0].props['filename']
        else:
            testfile = testtask[0].props['filename']
            reffile = reftask[0].props['filename']
        if len(testtask) != len(reftask):
            unit = "sectors" if sectors else "observables"
            raise Exception("Comparison Error: test and reference data have different number of " + unit)
        if not sectors:
            testtask, reftask = [testtask], [reftask]
        results = []
        for testsector, refsector in zip(testtask, reftask):
            testsector = _select_observables(testsector, whatlist, testfile, "test")
            refsector = _select_observables(refsector, whatlist, reffile, "reference")
            for testobs, refobs in zip(testsector, refsector):
                tol, diff = compare(testobs, refobs, factor)
                results.append(obsdict(tol, diff, testobs.props))
        comparisons.append(results)
    return all(obs['passed'] for task in comparisons for obs in task), comparisons


def compareMC(testfiles, reffiles, tol_factor='auto', whatlist=None):
    """Compare Monte Carlo results; return success and per-task comparisons."""
    return _compare_data(pyalps.loadMeasurements(testfiles), pyalps.loadMeasurements(reffiles),
                         2.0 if tol_factor == 'auto' else tol_factor, whatlist, _compare_mc)


def compareMixed(testfiles, reffiles, tol_factor='auto', whatlist=None):
    """Compare mixed Monte Carlo and deterministic results (QWL, DMRG)."""
    testdata = pyalps.loadMeasurements(testfiles)
    refdata = pyalps.loadMeasurements(reffiles)
    if len(testdata) != len(refdata):
        raise Exception("Comparison Error: test and reference data differ in number of tasks")
    try:
        testeig = pyalps.loadEigenstateMeasurements(testfiles)
        refeig = pyalps.loadEigenstateMeasurements(reffiles)
        for test, ref, teig, reig in zip(testdata, refdata, testeig, refeig):
            test += teig
            ref += reig
    except RuntimeError:
        pass
    return _compare_data(testdata, refdata, 2.0 if tol_factor == 'auto' else tol_factor,
                         whatlist, _compare_mixed)


def compareEpsilon(testfiles, reffiles, tol_factor='auto', whatlist=None):
    """Compare diagonalization results; return success and per-task comparisons."""
    testdata = pyalps.loadEigenstateMeasurements(testfiles)
    refdata = pyalps.loadEigenstateMeasurements(reffiles)
    if not testdata or not refdata:
        if not testdata:
            print("loadEigenstateMeasurements of file %s returned an empty list" % testfiles)
        if not refdata:
            print("loadEigenstateMeasurements of file %s returned an empty list" % reffiles)
        return
    return _compare_data(testdata, refdata, 1.0 if tol_factor == 'auto' else tol_factor,
                         whatlist, _compare_epsilon, sectors=True)

#*******************************************************************************

def detectDataType( fname ):

    fname = pyalps.make_list( fname )

    # Monte Carlo results
    try:
        data = pyalps.loadMeasurements( fname )
        if len(data[0]) == 0:
            raise RuntimeError

        for task in data:
            for obs in task:
                tmp = obs.y[0].error

    except (RuntimeError, AttributeError, IndexError):
        pass
    else:
        return compareMC

    # mixed type (QWL)
    try:
        data = pyalps.loadMeasurements( fname )
        if len(data[0]) == 0:
            raise RuntimeError

    except (RuntimeError, AttributeError, IndexError):
        pass
    else:
        return compareMixed

    # Epsilon-precise results
    try:
        data = pyalps.loadEigenstateMeasurements( fname )
    except RuntimeError:
        pass    
    else:
        return compareEpsilon

    raise Exception("Measurement data type couldn't be detected")

def checkProperties( testfile, reffile ):
    
    ll = pyalps.load.Hdf5Loader()
    try:
        tprop = ll.GetProperties( pyalps.make_list(testfile) )[0].props
        rprop = ll.GetProperties( pyalps.make_list(reffile)  )[0].props

    except IndexError:
        print("""WARNING: Simulation properties of files %s and %s 
        couldn't be loaded. Check for matching properties
        has been skipped""" % ( testfile, reffile ))
        return True

    if 'SEED' in tprop and 'SEED' in rprop:
        del tprop['SEED']
        del rprop['SEED']

    if 'filename' in tprop and 'filename' in rprop:
        del tprop['filename']
        del rprop['filename']

    if tprop.keys() != rprop.keys():
        return False
    return all(np.array_equal(tprop[key], rprop[key]) for key in tprop)

def compareTest( testinputfile, outputs, tmpdir, tstart, compMethod='auto' ):

    # read input
    testinputfile = os.path.expandvars( testinputfile )
    props, inputs, refFileList, whatlist = read_testprop_xml(testinputfile)
    try:
        tol = float(props['TOLERANCE'])
    except ValueError:
        tol = 'auto'

    testname = props['TESTNAME']
    #if 'yes' in props['WRITE_RESULTS_TO_FILE'].lower():
        # Redirect stdout to a file if file output is desired
        #sys.stdout = open( testname+'.testout.xml','w' ) 
    # Redirect stdout to a file to include all information in the .testout.xml file
    outfile = os.path.join(tmpdir, testname+'.testout.xml')
    fobj = open(outfile, 'w') 

    fobj.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    fobj.write('<?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>\n')
    fobj.write('<TEST xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.comp-phys.org/2002/10/ALPS.xsd">\n')
    fobj.write('<INPUT file="%s"/>\n' % testinputfile)

    # call comparison functions
    test_success = True
    if compMethod == 'auto':
        missmatch = [(tf,rf) for (tf,rf) in zip(outputs,refFileList) if not checkProperties(tf,rf)]
        if missmatch:
            for ft in missmatch:
                print("Simulation parameter mismatch in files %s and %s" % (ft[0],ft[1]))
                test_success = False
        else:
            for tfile, rfile in zip(outputs, refFileList):
                # call the function object returned by detectDataType()
                Success, compare_list = detectDataType(tfile)([tfile], [rfile], tol, whatlist = whatlist)
                writeTest2stdout(compare_list, fobj)
                if not Success: test_success = False
    else:
        missmatch = [(tf,rf) for (tf,rf) in zip(outputs,refFileList) if not checkProperties(tf,rf)]
        if missmatch:
            for ft in missmatch:
                print("Simulation parameter mismatch in files %s and %s" % (ft[0],ft[1]))
                test_success = False
        else:
            for tfile, rfile in zip(outputs, refFileList):
                # call the function specified by compMethod (input by user)
                string2fobj = { 'compareMC' : compareMC, 'compareMixed' : compareMixed, \
                                'compareEpsilon' : compareEpsilon }
                Success, compare_list = string2fobj[compMethod](outputs, refFileList, tol, whatlist = whatlist)
                writeTest2stdout(compare_list, fobj)
                if not Success: test_success = False

    fobj.write('<EXECUTED>\n')
    fobj.write('    <FROM>%s</FROM>\n' % tstart)
    fobj.write('    <TO>%s</TO>\n' % strftime("%Y-%b-%d %H:%M:%S", gmtime()))
    fobj.write('    <MACHINE><NAME>%s</NAME></MACHINE>\n' % socket.gethostname())
    fobj.write('</EXECUTED>\n')


    # if test unsuccessful, archive all test data
    if not test_success and 'yes' in props['SAVE_OUT_IF_FAIL'].lower():
        archive_name = props['TESTNAME'] + '.failed'
        shutil.make_archive(archive_name, 'bztar', tmpdir)
        fobj.write('<OARCHIVE file="%s"/>\n' % (archive_name+'.tar.bz2'))

    fobj.write('</TEST>\n')

    # reset stdout, otherwise interactive usage doesn't work anymore
    fobj.close()
    #sys.stdout.close()
    #sys.stdout = sys.__stdout__

    # write test results to a file or to stdout
    if 'yes' in props['WRITE_RESULTS_TO_FILE'].lower():
        shutil.copy(outfile, '.')
    else:
        f = open(outfile)
        print(f.read())
        f.close()
    return test_success


def runTest(testinputfile, outputs='auto', compMethod='auto', pyexec='auto'):
    """ Run the test according to testinputfile and compare
        its output against reference files

    inputs are:
    -----------
    script: a script to be run

    testinputfile: file containing reference info

    outputs: - default is 'auto' -> assume file names are identical with reference files
           list of script outputs
         - glob pattern
         - explicit list of files ( order has to be the same as in the .testin.xml file )

    compMethod: - default is 'auto' -> try to detect comparison Method
              to be used (e.g. Monte Carlo or epsilon precise data)
    """
    tstart = strftime("%Y-%b-%d %H:%M:%S", gmtime())
    testinputfile = os.path.expandvars(testinputfile)
    props, inputs, refFileList, whatlist = read_testprop_xml(testinputfile)
    script = props['SCRIPT']

    tmpdir = tempfile.mkdtemp()

    # execute given script in tmpdir
    pardir = os.getcwd()
    shutil.copy(script, tmpdir)
    if inputs is not None:
        for f in inputs: shutil.copy(f, tmpdir)
    
    # get python executable
    if pyexec == 'auto':
        pyexec = sys.executable
    
    os.chdir( tmpdir )
    cmdline = [pyexec, os.path.basename(script)]
    pyalps.executeCommand(cmdline)
    if inputs is not None:
        for f in inputs: os.remove(f)

    os.chdir(pardir)
    
    # copy stylesheet
    pyalps.tools.copyStylesheet(pardir)

    # Guess outputs from reference files
    if outputs == 'auto':
        outputs = [ os.path.join( tmpdir, os.path.basename(x) ) for x in refFileList ]

    elif type(outputs) == str:
        print("Using glob '%s' to find outputs" % outputs)
        outputs = pyalps.getResultFiles( pattern=os.path.basename(outputs), \
            dirname=os.path.dirname(outputs) )

    if not outputs:
        print("\nList of output files of %s is empty\n" % script)
        test_success = False

    else:
        missing = [ x for x in outputs if not os.path.exists(x) ]
        if missing:
            for f in missing:
                print("Output file '%s' does not exist" % f)
            test_success = False
        else:
            # Start test
            test_success = compareTest( testinputfile, outputs, tmpdir, tstart, compMethod=compMethod )
    
    # if something goes wrong above, tmpdir will not be removed
    # for the moment this is useful, later maybe use try/except
    shutil.rmtree(tmpdir)
    return test_success

#*******************************************************************************
# Test creation
#*******************************************************************************

def writeTestInputFile( script, inputs, refparms, refFileList, what, dirname='.', postfix='' ):
    """ Write a *.testin.xml input file """

    refparms['TESTNAME'] = refparms['TESTNAME'].replace('.in.xml','') + postfix
    testinfile = os.path.join( os.path.expandvars(dirname), refparms['TESTNAME']+'.testin.xml')
    f = open( testinfile, 'w' )
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>\n')
    f.write('<TEST xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.comp-phys.org/2002/10/ALPS.xsd">\n') 

    f.write('   <OUTPUT file="%s"/>\n' % (refparms['TESTNAME'] + '.testout.xml'))
    f.write('   <SCRIPT file="%s"/>\n' % script)
    if inputs is not None:
        f.write('   <SCRIPTINPUTFILES>\n')
        for inpf in inputs:
            f.write('       <SCRIPTINPUTFILE file="%s"/>\n' % inpf)

        f.write('   </SCRIPTINPUTFILES>\n')
        
    f.write('   <TESTINFOS>\n')
    for key in refparms:
        f.write('       <TESTINFO name="%s">%s</TESTINFO>\n' % (key, refparms[key]) )
        #f.write('      <TESTINFO name="TOLERANCE">%s</TESTINFO>\n' % refparms['TOLERANCE'])

    f.write('   </TESTINFOS>\n')

    f.write('   <REFERENCE>\n')
    for reffile in refFileList:
        f.write('       <FILE file="%s" />\n' % reffile )

    f.write('   </REFERENCE>\n')

    if what is not None:

        f.write('   <OBSERVABLES>\n')
        for obs in what:
            f.write('       <OBSERVABLE>' + obs + '</OBSERVABLE>\n')

        f.write('   </OBSERVABLES>\n')

    f.write('</TEST>\n')
    f.close()

    return testinfile

def createTest( script, inputs=None, outputs=None, prefix=None, refdir='./ref' ):
    """ Create reference data, .testin.xml file and execute_test.py

    inputs are:
    -----------
    script: computes results to be tested 

    inputs: Optional list of input files if the application(s)
            called in 'script' rely on them and the input files are in the
            same directory as 'script'. If you specified
            relative paths to another directory, it won't work.

    outputs or prefix: outputs of script can either be specified with
               a complete list of output files or as a prefix 

    creates a script called apptest_name_of_script.py, which can be used to execute the test
    """
 
    if outputs is not None and prefix is not None:
        raise Exception("Cannot both define outputs and prefix")
    elif outputs is None and prefix is None:
        raise Exception("Script output has to be specified")
    script = os.path.expandvars(script)
    scriptdir = os.path.dirname(script)

    if not os.path.exists(refdir):  recursive_mkdir(refdir)

    # Copy input files to refdir to allow execution of script there
    if inputs is not None:

        for f in inputs:
            if not os.path.expandvars(os.path.dirname(f)) == scriptdir:
                print("Input files to %s should be in the same directory as %s" % (script, script))
                sys.exit(1)

            shutil.copy(f, refdir)

    # execute given script in refdir ( creates reference data )
    pardir = os.getcwd()
    os.chdir(refdir)
    cmdline = [ sys.executable, os.path.join(pardir, script) ]
    pyalps.executeCommand(cmdline)
    if inputs is not None:
        for f in inputs: os.remove(f)
    os.chdir(pardir)

    if prefix is None:
        reffiles = [ os.path.join(refdir, os.path.basename(f)) for f in outputs ]
    else:
        reffiles = pyalps.getResultFiles(prefix=prefix, dirname=refdir)

    if not reffiles:
        print("Reference files not found. (If you use 'loop' or 'dmrg', try to delete old result files.)")
        sys.exit(1)

    # acquire a list of all observables
    allobs = []
    try:
        eigenstatedata = pyalps.loadEigenstateMeasurements(reffiles)   
    except RuntimeError:
        pass
    else:
        try:
            allobs += [ o.props['observable'] for o in eigenstatedata[0][0] ]

        # DMRG eigenstate data has one level of nesting less
        except TypeError:
            allobs += [ o.props['observable'] for o in eigenstatedata[0] ]

    try:
        mcdata = pyalps.loadMeasurements(reffiles)
    except RuntimeError:
        pass
    else:
        allobs += [ o.props['observable'] for o in mcdata[0] ]

    allobs = list(set(allobs))

    scriptname = os.path.basename(script)
    scriptname = os.path.splitext( scriptname )[0]
    scriptname_prefixed = 'apptest_%s.py' % scriptname

    # Write .xml test-input file
    refparms = {
        "TESTNAME"                : scriptname,
        "TOLERANCE"                : "auto",
        "WRITE_RESULTS_TO_FILE" : "yes",
        "SAVE_OUT_IF_FAIL"      : "yes"
    }

    testinputfile = writeTestInputFile(script, inputs, refparms, reffiles, allobs)
    pyalps.tools.copyStylesheet(pardir)

    # Write .py test-start script
    f = open( scriptname_prefixed, 'w' )
    f.write( '#!/usr/bin/env python\n\n' )
    f.write( 'import sys\n' )
    f.write( 'from pyalps import apptest\n' )

    f.write('# Explicitly specify "compMethod=..." and "outputs=..." if needed\n')
    f.write("ret = apptest.runTest( '%s', outputs='auto', compMethod='auto', pyexec='auto' )\n" % testinputfile)
    f.write('if not ret: sys.exit(1)\n')

    f.close()
    os.chmod(scriptname_prefixed, 0o755)
