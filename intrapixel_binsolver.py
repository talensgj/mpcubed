#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt

from coordinate_grids import PolarGrid
from index_functions import index_statistics

fLC = '/data2/talens/Jul2015/fLC_20150716LPC.hdf5'
red = '/data2/talens/Jul2015/red_20150716LPC.hdf5'

from scipy.optimize import minimize
from scipy.ndimage.filters import gaussian_filter1d

from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def binsolver(ind1, ind2, y, flux, eflux):
    
    cs = np.cos(2*np.pi*y)
    sn = np.sin(2*np.pi*y)
    w = 1/eflux**2
    
    a = .1*np.ones(np.amax(ind2)+1)
    b = np.zeros(np.amax(ind2)+1)
    
    for i in range(600):
        print i
        f = np.bincount(ind1, flux*(a[ind2]*cs-b[ind2]*sn+1)*w)/np.bincount(ind1, (a[ind2]*cs-b[ind2]*sn+1)**2*w)
        a = np.bincount(ind2, (flux+f[ind1]*b[ind2]*sn-f[ind1])*f[ind1]*cs*w)/np.bincount(ind2, (f[ind1]*cs)**2*w)
        b = np.bincount(ind2, (flux-f[ind1]*a[ind2]*cs-f[ind1])*f[ind1]*sn*w)/np.bincount(ind2, -(f[ind1]*sn)**2*w)
    
        if i == 0:
            solo = f[ind1]*(a[ind2]*cs-b[ind2]*sn+1)
            aold = np.copy(a) 
            bold = np.copy(b)
            
        else:
            sol = f[ind1]*(a[ind2]*cs-b[ind2]*sn+1)
            crit = np.nanmax((sol-solo)/solo)
            
            crita = np.nanmax((a-aold)/aold)
            critb = np.nanmax((b-bold)/bold)
            
            print crita, critb
            if (crit<1e-3)&(crita<1e-3)&(critb<1e-3):
                break
                
            solo = np.copy(sol)
            aold = np.copy(a)
            bold = np.copy(b)
    
    return f, a, b
    
pg = PolarGrid(270, 720)

with h5py.File(fLC, 'r') as f:
    
    ascc = f['table_header/ascc'].value
    ra = f['table_header/ra'].value
    dec = f['table_header/dec'].value
    nobs = f['table_header/nobs'].value.astype('int')

decidx, decuni = pg.find_decidx(dec, compact=True)

ascc = ascc[decuni==82]
ra = ra[decuni==82]
dec = dec[decuni==82]
nobs = nobs[decuni==82]

ndata = np.sum(nobs)
select = np.append(0, np.cumsum(nobs))

flux = np.zeros(ndata)
eflux = np.zeros(ndata)
y = np.zeros(ndata)
lst = np.zeros(ndata)

with h5py.File(fLC, 'r') as f:
    
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]].value
        
        flux[select[i]:select[i+1]] = lc['flux0']
        eflux[select[i]:select[i+1]] = lc['eflux0']
        y[select[i]:select[i+1]] = lc['y']
        lst[select[i]:select[i+1]] = lc['lst']
        
ha = np.mod(lst*15.-np.repeat(ra, nobs), 360.)
haidx, hauni = pg.find_raidx(ha, compact=True)
staridx = np.repeat(np.arange(len(ascc)), nobs)

idx = np.ravel_multi_index((hauni, staridx), (len(haidx), len(ascc)))

f, a, b = binsolver(idx, hauni, y, flux, eflux)

plt.subplot(211)
plt.plot(np.mod(haidx-135, 270), np.abs(np.arctan(b/a)), '.')
plt.subplot(212)
plt.plot(np.mod(haidx-135, 270), np.sqrt(a**2+b**2), '.')
plt.show()

solution = f[idx]*(a[hauni]*np.cos(2*np.pi*y)-b[hauni]*np.sin(2*np.pi*y)+1)
solution1 = f[idx]
solution2 = (a[hauni]*np.cos(2*np.pi*y)-b[hauni]*np.sin(2*np.pi*y)+1)

for i in range(len(ascc)):
    here = staridx == i
    
    ax = plt.subplot(311)
    plt.plot(flux[here], '.')
    plt.plot(solution[here], '.')
    
    ax2 = plt.subplot(312, sharex=ax)
    plt.plot(flux[here]/solution1[here], '.')
    plt.plot(solution2[here], '.')
    
    plt.subplot(313, sharex=ax, sharey=ax2)
    plt.plot(flux[here]/solution[here], '.')
    
    plt.show()
    plt.close()
  
    
