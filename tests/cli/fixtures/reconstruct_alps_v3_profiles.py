"""Reconstruct the pinned ALPS v3.0.0 writer layouts without building ALPS.

Run with Python, h5py and NumPy; this intentionally records reconstruction
provenance rather than claiming output from a compiled release.
"""

# SPDX-License-Identifier: MIT

import hashlib,json
from pathlib import Path
import h5py,numpy as np
folder=Path(__file__).resolve().parent
fixture=folder/'alps-v3.0.0-profiles.h5'
ascii_string=h5py.string_dtype('ascii')
def string(group,name,value):
    group.create_dataset(name,data=value.encode('utf-8'),dtype=ascii_string)
def attr(group,name,value):
    if isinstance(value,str): group.attrs.create(name,value.encode('utf-8'),dtype=ascii_string)
    else: group.attrs[name]=value
with h5py.File(fixture,'w') as f:
    p=f.create_group('ngs-parameters')
    p['N']=np.int32(7); p['BETA']=np.float64(2.5); p['ENABLED']=np.int8(1)
    p['VECTOR']=np.array([1,2,3],dtype='i4')
    p['EMPTY_INT']=h5py.Empty('i4'); p['EMPTY_REAL']=h5py.Empty('f8')
    p['EMPTY_TEXT']=h5py.Empty(ascii_string)
    string(p,'TITLE','Δ experiment')
    string(p,'EXPRESSION','sqrt(2) + unknown')
    d=p.create_dataset('DELTA',data=np.array([1.25,-0.0],dtype='f8'))
    d.attrs['__complex__']=np.int8(1)
    p=f.create_group('legacy-parameters')
    p['L']=np.int32(4); p['T']=np.float64(.5)
    string(p,'EXPRESSION','2 * unresolved')
    s=f.create_group('scalar-mcdata')
    s['count']=np.uint64(0); s['mean/value']=0.; s['mean/error']=0.
    attr(s,'cannotrebin',np.int8(0))
    d=s.create_dataset('timeseries/data',data=h5py.Empty('f8'))
    for name,value in [('binsize',np.uint64(0)),('maxbinnum',np.uint64(128)),('binningtype','linear')]: attr(d,name,value)
    d=s.create_dataset('jacknife/data',data=h5py.Empty('f8')); attr(d,'binningtype','linear')
    s=f.create_group('scalar-evaluator')
    s['count']=np.float64(6); s['mean/value']=1.; s['mean/error']=.25
    s['mean/error_convergence']=np.int32(0)
    attr(s,'changed',np.int8(0)); attr(s,'nonlinearoperations',np.int8(1))
    for field,data in [('data',np.array([1.,2.,3.])),('data2',h5py.Empty('f8'))]:
        d=s.create_dataset('timeseries/'+field,data=data)
        for name,value in [('discard',np.uint32(0)),('maxbinnum',np.uint32(128)),('binningtype','linear')]: attr(d,name,value)
    s=f.create_group('scalar-observable')
    s['count']=np.uint64(1); s['mean/value']=2.5
    string(s,'labels','Δ energy')
    s=f.create_group('vector-observable')
    s['count']=np.uint64(1); s['mean/value']=np.array([1.,2.]); s['sum']=np.array([1.,2.]); s['sum2']=np.array([1.,4.])
    s.create_dataset('labels',data=np.array([b'x','λ'.encode()],dtype=object),dtype=ascii_string)
    for field,data in [('logbinning',[[1.,2.]]),('logbinning2',[[1.,4.]]),('logbinning_lastbin',[[1.,2.]]),('logbinning_counts',np.array([1],dtype='u8'))]:
        d=s.create_dataset('timeseries/'+field,data=data); attr(d,'binningtype','logarithmic')
    for field,data in [('partialbin',[1.,2.]),('partialbin2',[1.,4.])]:
        d=s.create_dataset('timeseries/'+field,data=data); attr(d,'count',np.uint32(1))
    for field in ['data','data2']:
        d=s.create_dataset('timeseries/'+field,data=h5py.Empty('i4'))
        for name,value in [('binsize',np.uint32(1)),('minbinsize',np.uint32(1)),('maxbinnum',np.uint32(128)),('binningtype','linear')]: attr(d,name,value)
    s=f.create_group('unrecoverable-empty-vector-observable'); s['count']=np.uint64(0)
    for field in ['logbinning','logbinning2','logbinning_lastbin','data','data2']:
        s['timeseries/'+field]=h5py.Empty('i4')
    s['timeseries/logbinning_counts']=h5py.Empty('u8')
metadata={
  'generator':Path(__file__).name,'fixture':fixture.name,'sha256':hashlib.sha256(fixture.read_bytes()).hexdigest(),
  'provenance':'Reconstructed with h5py from the pinned C++ writer contracts; not emitted by a compiled ALPS release.',
  'release':'ALPS v3.0.0','source_revision':'1950cc6f682d7c4c1deae8b816f283857b1819d1',
  'source_url':'https://github.com/ALPSim/ALPS/tree/v3.0.0',
  'writers':{
    'ngs-parameters':'src/alps/ngs/lib/params.cpp:114-116; src/alps/ngs/lib/paramvalue.cpp:115-119',
    'legacy-parameters':'src/alps/parameter/parameters.C:126-145',
    'scalar-mcdata':'src/alps/alea/mcdata.hpp:462-489',
    'scalar-evaluator':'src/alps/alea/simpleobsdata.h:935-972',
    'scalar-observable':'src/alps/alea/abstractsimpleobservable.h:84; src/alps/alea/abstractsimpleobservable.ipp:30-46 (base observable fields, scalar label_type=std::string)',
    'vector-observable':'src/alps/alea/abstractsimpleobservable.ipp:30-59; src/alps/alea/simplebinning.h:652-669; src/alps/alea/detailedbinning.h:286-327',
    'unrecoverable-empty-vector-observable':'Same binning writers with count=0 and no component-shape exemplar.'},
  'encoding_notes':'Old archive native std::string is variable-length H5T_C_S1 (ASCII declaration, UTF-8 bytes). Native bool is signed byte. Empty primitive vectors use typed NULL; empty vectors of valarrays use INT NULL. Floating values are doubles, int is represented as signed 32-bit and count widths are explicit in the fixture.',
  'coverage_boundary':'Released writer layouts only. Flagless legacy sum-to-mean semantics, arbitrary lost array ranks, and private WIP params.v1 envelopes are outside the automatic profile.'}
(folder/'alps-v3.0.0-profiles.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(fixture,fixture.stat().st_size)
