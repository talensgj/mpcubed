#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from coordinate_grids import HealpixGrid
from index_functions import index_statistics

import sysrem_dev

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 
from matplotlib.ticker import FuncFormatter, MultipleLocator, NullLocator
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def periodic(x, pos):
        'The two args are the value and tick position'
        x = x%24
        x1 = np.floor(x)
        x2 = (x-np.floor(x))*60
        
        return '%02.fh%02.fm' % (x1, x2)

formatter = FuncFormatter(periodic)
hg = HealpixGrid(8)

with h5py.File('/data2/talens/Jul2015/fLC_20150714LPC.hdf5', 'r') as f, h5py.File('/data2/talens/Jul2015/red_20150714LPC.hdf5', 'r') as g:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    nobs = f['header_table/nobs'].value.astype('int')
    
    skyidx = hg.find_gridpoint(ra, dec)
    print skyidx[ascc=='807144']

    here = (skyidx == 266)
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    nobs = nobs[here]
    
    lstidx = np.array([])
    flux = np.array([])
    eflux = np.array([])
    sflux = np.array([])
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]]
        rc = g['data/'+ascc[i]]
        lstidx = np.append(lstidx, lc['lstidx'])
        flux = np.append(flux, rc['ipcflux0'])
        eflux = np.append(eflux, rc['eipcflux0'])
        sflux = np.append(sflux, rc['sipcflux0'])
    
eflux = sflux

lstidx = lstidx.astype('int')
staridx = np.repeat(np.arange(len(ascc)), nobs)

here = (flux>0)&(eflux>0)
staridx = staridx[here]
lstidx = lstidx[here]
flux = flux[here]
eflux = eflux[here]

data = np.full((len(ascc), 13500), fill_value=np.nan)
data[staridx, lstidx] = flux
data = data/np.nanmedian(data, axis=1, keepdims=True)

sort = np.argsort(ra)
data = data[sort]

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_locator(NullLocator())
plt.title(r'14-07-2015 LPC')
plt.imshow(data, aspect='auto', cmap=viridis, extent=(0,24,0,1), vmin=0, vmax=1.2)
cb = plt.colorbar()
plt.xlim(16.5,23)
cb.set_label('Normalized Flux')
plt.xlabel('Local Sidereal Time')
plt.tight_layout()
#plt.show()
plt.savefig('Pipeline_part2_0_data.png')

F, S, niter, chisq, npoints, npars = sysrem_dev.sysrem(staridx, lstidx, flux, eflux, verbose=True)

x = np.arange(len(S))
x = x/13500.*24.

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
plt.title(r'$\rm{niter}=%i$, $\chi^{2}=%.2f$'%(niter, chisq))
plt.plot(x, S, '.')
plt.xlim(16.5,23)
plt.ylim(0,1.2)
plt.xlabel('Local Sidereal Time')
plt.ylabel('Sky')
#plt.show()
plt.savefig('Pipeline_part2_1_sky.png')

data = np.full((len(ascc), 13500), fill_value=np.nan)
data[staridx, lstidx] = flux/(F[staridx]*S[lstidx])
data = data/np.nanmedian(data, axis=1, keepdims=True)

sort = np.argsort(ra)
data = data[sort]

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_locator(NullLocator())
plt.title(r'14-07-2015 LPC')
plt.imshow(data, aspect='auto', cmap=viridis, extent=(0,24,0,1), vmin=0, vmax=1.2)
cb = plt.colorbar()
plt.xlim(16.5,23)
cb.set_label('Reduced Flux')
plt.xlabel('Local Sidereal Time')
plt.tight_layout()
#plt.show()
plt.savefig('Pipeline_part2_3_result.png')
