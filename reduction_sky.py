#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
from reduction_functions import make_array, make_transmission_map

import matplotlib.pyplot as plt

import healpy

# Read data.
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
    flux0 = data['Flux0'].value
    eflux0 = data['eFlux0'].value
    flags = data['Flag'].value

# Create HA indices.
ra_id = np.floor(ra/15.*3600./6.4).astype('int')
ha_id = lst_id - np.repeat(ra_id, Nobs)
ha_id = np.mod(ha_id, 13500)

# Create Dec indices.
dec_id = np.digitize(dec, bins=np.linspace(-90, 90, 361))

# Read transmission map.
with h5py.File('/data2/talens/Oct2014LPE_Trans.hdf5') as f:
    dset = f['20141009/Transmission']
    transmission = dset.value
    offset_ha = dset.attrs['offset_HA']
    offset_dec = dset.attrs['offset_Dec']
    flags_t = f['20141009/Flags'].value
    
# Expand the transmission map to match the 1d data vector.
ind1 = np.repeat(dec_id-offset_dec, Nobs)
ind2 = ha_id-offset_ha

# Need to deal with indices not in transmission!!!!!!

# First indices outside of transmission array.
bad_ind = np.where((ind1 < 0) | (ind2 < 0) | (ind1 >= transmission.shape[0]) | (ind2 >= transmission.shape[1]))
ind1[bad_ind] = 0 # Set to existing index
ind2[bad_ind] = 0 # Set to existing index

tcurve = transmission[ind1, ind2] 
tcurve[bad_ind] = np.nan # Set the bad ones to NaN, so everywhere no transmission is known is now NaN

# This now includes all data for which no transmission value was known
bad_ind = np.where(np.isnan(tcurve))

fcurve = flags_t[ind1, ind2]
fcurve[np.isnan(fcurve)] = 0 
fcurve[bad_ind] += 2

plt.hist(fcurve, bins=np.linspace(-0.5,3.5,5), log=True)
plt.show()

exit()

# Remove the transmission effect.
tflux0 = flux0/tcurve
etflux0 = eflux0/tcurve

# Set bad data to NaN.
here, = np.where(flags > 0)
flux0[here] = np.nan
eflux0[here] = np.nan

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

transmission, normalization, niter, chi2 = make_transmission_map(data, error, binnum-offset_binnum)

binnum = healpy.ang2pix(16, (np.repeat(dec, Nobs)+90)*np.pi/180., np.repeat(ra, Nobs)*np.pi/180.)
_, binnum = np.unique(binnum, return_inverse=True)

scurve = transmission[binnum, time_id]

stflux0 = tflux0/scurve
estflux0 = etflux0/scurve

with h5py.File('/data2/talens/Oct2014LPE_Result.hdf5') as f:

    grp = f.create_group('Pipeline')
    grp.create_dataset('tFlux0', data=tflux0)
    grp.create_dataset('etFlux0', data=etflux0)
    grp.create_dataset('tcurve', data=tcurve)
    grp.create_dataset('stFlux0', data=stflux0)
    grp.create_dataset('estFlux0', data=estflux0)
    grp.create_dataset('scurve', data=scurve)

exit()
binnum = np.unique(binnum)
for i in range(transmission.shape[1]):
    result = np.full(healpy.nside2npix(16), fill_value=np.nan)
    result[binnum] = transmission[:,i]
    healpy.mollview(result, min=0.5, max=1.5)
    plt.savefig('/data2/talens/cloud_%.4i.png'%i)
    plt.close()

exit()

    
    
    
    
