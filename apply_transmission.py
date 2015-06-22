#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

with h5py.File('/data2/talens/Feb2015LPE.hdf5') as f:

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

# Create HA indices.
ra_id = np.floor(ra/15.*3600./6.4).astype('int')
ha_id = lst_id - np.repeat(ra_id, Nobs)
ha_id = np.mod(ha_id, 13500)

# Create Dec indices.
dec_id = np.digitize(dec, bins=np.linspace(-90, 90, 1441))

# Read transmission map.
with h5py.File('/data2/talens/Feb2015LPE_Trans.hdf5') as f:
    dset = f['20150203/Transmission']
    transmission = dset.value
    offset_ha = dset.attrs['offset_HA']
    offset_dec = dset.attrs['offset_Dec']
    flags_t = f['20150203/Flags'].value
    
ind1 = np.repeat(dec_id-offset_dec, Nobs)
ind2 = ha_id-offset_ha

# Find indices outside of transmission array.
bad_ind = np.where((ind1 < 0) | (ind2 < 0) | (ind1 >= transmission.shape[0]) | (ind2 >= transmission.shape[1]))
ind1[bad_ind] = 0 # Set to existing index
ind2[bad_ind] = 0 # Set to existing index

# Expand the map to match the data.
tcurve = transmission[ind1, ind2] 
fcurve = flags_t[ind1, ind2]

# Set the values of the bad indices to their defaults.
tcurve[bad_ind] = np.nan
fcurve[bad_ind] = 2

tflux0 = flux0/tcurve
etflux0 = eflux0/tcurve

with h5py.File('/data2/talens/Feb2015LPE_Result.hdf5') as f:
    grp = f.create_group('Pipeline')
    grp.create_dataset('tFlux0', data=tflux0)
    grp.create_dataset('etFlux0', data=etflux0)
    grp.create_dataset('Transmission', data=tcurve)
    grp.create_dataset('MyFlags', data=fcurve)
