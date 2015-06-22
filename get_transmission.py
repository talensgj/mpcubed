#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
from reduction_functions import make_array, make_transmission_map

import matplotlib.pyplot as plt

# Read data.
with h5py.File('/data2/talens/20150203LPE.hdf5') as f:

    stars = f['Stars']
    ascc = stars['ASCC'].value
    ra = stars['RAJ2000'].value
    dec = stars['DEJ2000'].value
    Nobs = stars['Nobs'].value.astype('int')

    data = f['Data']
    jd = data['JDmid'].value
    lst = data['lpst'].value
    lst_id = data['lpstid'].value.astype('int')
    flux0 = data['Flux0'].value
    eflux0 = data['eFlux0'].value
    flags = data['Flag'].value

# Set bad data to NaN.
here, = np.where(flags > 0)
flux0[here] = np.nan
eflux0[here] = np.nan

# Create HA indices.
ra_id = np.floor(ra/15.*3600./6.4).astype('int')
ha_id = lst_id - np.repeat(ra_id, Nobs)
ha_id = np.mod(ha_id, 13500)
offset_ha = np.amin(ha_id) # neglects periodicity...

# Cast the data and the errors to arrays suitable for sysrem.
data = make_array(flux0, Nobs, ha_id-offset_ha) # still neglecting periodicity...
error = make_array(eflux0, Nobs, ha_id-offset_ha)

# Create Dec indices.
dec_id = np.digitize(dec, bins=np.linspace(-90, 90, 1441))
offset_dec = np.amin(dec_id)

# Obtain the transmission map.
transmission, normalization, niter, chi2, flags = make_transmission_map(data, error, dec_id-offset_dec)

with h5py.File('/data2/talens/Feb2015LPE_Trans.hdf5') as f:
    
    dset = f.create_dataset('20150203/Transmission', data=transmission)
    dset.attrs['offset_HA'] = offset_ha
    dset.attrs['offset_Dec'] = offset_dec
    dset = f.create_dataset('20150203/Flags', data=flags)
    dset.attrs['offset_HA'] = offset_ha
    dset.attrs['offset_Dec'] = offset_dec


