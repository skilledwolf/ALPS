"""Inputs for the remaining legacy readers, independent of removed producers."""
import h5py
import numpy as np
import pytest


@pytest.fixture
def legacy_alea_file(tmp_path):
    def write(samples, path='/simulation/results/energy', raw=False):
        # Reconstruct the released MCData/base-observable field contracts with
        # ordinary HDF5 primitives. This is not a release-generated fixture.
        samples = np.asarray(samples, dtype=float)
        filename = str(tmp_path/'legacy.h5')
        with h5py.File(filename, 'a') as archive:
            group = archive.create_group(path)
            group['count'] = np.uint64(len(samples))
            group['mean/value'] = samples.mean(axis=0)
            if not raw:
                group.attrs['cannotrebin'] = False
                group['mean/error'] = samples.std(axis=0, ddof=1)/np.sqrt(len(samples))
                group['variance/value'] = samples.var(axis=0, ddof=1)
                group['tau/value'] = np.zeros(samples.shape[1:])
                data = group.create_dataset('timeseries/data', data=samples)
                data.attrs.update(binsize=np.uint64(1), maxbinnum=np.uint64(len(samples)), binningtype='linear')
        return filename
    return write
