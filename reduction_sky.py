#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
from reduction_functions import make_array, make_transmission_map

import matplotlib.pyplot as plt

import healpy

with h5py.File('/data2/talens/Oct2014LPE.hdf5') as f:

    stars = f['Stars']
    ascc = stars['ASCC'].value
    ra = stars['RAJ2000'].value
    dec = stars['DEJ2000'].value
    Nobs = stars['Nobs'].value.astype('int')

    data = f['Data']
    jd = data['JDmid'].value
    lst = data['lpst'].value
    lst_id = data['lpstid'].value.astype('int')
    flags = data['Flag'].value

with h5py.File('/data2/talens/Oct2014LPE_Result.hdf5') as f:
    tflux0 = f['Pipeline/tFlux0'].value
    etflux0 = f['Pipeline/etFlux0'].value
    myflags = f['Pipeline/MyFlags'].value

# Set bad data to NaN.
here, = np.where((flags > 0)|(myflags > 0))
tflux0[here] = np.nan
etflux0[here] = np.nan

# Modify the LST indices.
lst_id = np.mod(lst_id-13500/2, 13500) # bad periodicity hack...
jd_id = np.floor(jd).astype('int')-2456931
time_id = np.ravel_multi_index((jd_id, lst_id), (31, 13500)) # I think this is unique (but only because we don't observe 24hours), but it's not ordered. NEED AN ALTERNATIVE
_, time_id = np.unique(time_id, return_inverse=True)

# Cast the data and the errors to arrays suitable for sysrem.
data = make_array(tflux0, Nobs, time_id) # still neglecting periodicity...
error = make_array(etflux0, Nobs, time_id)

# Create grid on the sky.
binnum = healpy.ang2pix(16, (dec+90)*np.pi/180., ra*np.pi/180.)
offset_binnum = np.amin(binnum)

transmission, normalization, niter, chi2, sflags = make_transmission_map(data, error, binnum-offset_binnum)

binnum = healpy.ang2pix(16, (np.repeat(dec, Nobs)+90)*np.pi/180., np.repeat(ra, Nobs)*np.pi/180.)
_, binnum = np.unique(binnum, return_inverse=True)

scurve = transmission[binnum, time_id]
myflags += 4*sflags[binnum, time_id]

stflux0 = tflux0/scurve
estflux0 = etflux0/scurve

with h5py.File('/data2/talens/Oct2014LPE_Result.hdf5') as f:

    grp = f['Pipeline']
    grp.create_dataset('stFlux0', data=stflux0)
    grp.create_dataset('estFlux0', data=estflux0)
    grp.create_dataset('Sky', data=scurve)
    grp.create_dataset('MyFlags2', data=myflags)
