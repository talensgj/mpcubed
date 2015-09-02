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

def binsolver(ind1, ind2, ind3, y, flux, eflux):
    
    cs = np.cos(2*np.pi*y)
    sn = np.sin(2*np.pi*y)
    w = 1/eflux**2
    
    F = np.bincount(ind1, flux*w)/np.bincount(ind1, w)
    T = np.bincount(ind2, flux*F[ind1]*w)/np.bincount(ind2, F[ind1]**2*w)
    
    a = np.zeros(np.amax(ind3)+1)
    b = .1*np.ones(np.amax(ind3)+1)
    
    for i in range(50):
        print i
        F = np.bincount(ind1, flux*T[ind2]*(a[ind3]*sn+b[ind3]*cs+1)*w)/np.bincount(ind1, (T**2)[ind2]*(a[ind3]*sn+b[ind3]*cs+1)**2*w)
        T = np.bincount(ind2, flux*F[ind1]*(a[ind3]*sn+b[ind3]*cs+1)*w)/np.bincount(ind2, (F**2)[ind1]*(a[ind3]*sn+b[ind3]*cs+1)**2*w)
        
        a = np.bincount(ind3, (flux-F[ind1]*T[ind2]*(b[ind3]*cs+1))*F[ind1]*T[ind2]*sn*w)/np.bincount(ind3, (F[ind1]*T[ind2]*sn)**2*w)
        b = np.bincount(ind3, (flux-F[ind1]*T[ind2]*(a[ind3]*sn+1))*F[ind1]*T[ind2]*cs*w)/np.bincount(ind3, (F[ind1]*T[ind2]*cs)**2*w)
    
        if i == 0:
            solo = F[staridx]*T[hauni2]*(a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1)
            aold = np.copy(a) 
            bold = np.copy(b)
            
        else:
            sol = F[staridx]*T[hauni2]*(a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1)
            crit = np.nanmax((sol-solo)/solo)
            
            crita = np.nanmax((a-aold)/aold)
            critb = np.nanmax((b-bold)/bold)
            
            print crita, critb
            if (crit<1e-3)&(crita<1e-3)&(critb<1e-3):
                break
                
            solo = np.copy(sol)
            aold = np.copy(a)
            bold = np.copy(b)
            
    return F, T, a, b
    

with h5py.File(fLC, 'r') as f:
    
    ascc = f['table_header/ascc'].value
    ra = f['table_header/ra'].value
    dec = f['table_header/dec'].value
    nobs = f['table_header/nobs'].value.astype('int')

pg = PolarGrid(270, 720)
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
haidx1, hauni1 = pg.find_raidx(ha, compact=True)
staridx = np.repeat(np.arange(len(ascc)), nobs)

pg = PolarGrid(13500, 720)
haidx2, hauni2 = pg.find_raidx(ha, compact=True)

F, T, a, b = binsolver(staridx, hauni2, hauni1, y, flux, eflux)

sol = F[staridx]*T[hauni2]*(a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1)

np.savez('joint_solution.npz', camtrans=T, solution=sol)
exit()




plt.subplot(311)
plt.plot(T, '.')
plt.subplot(312)
plt.plot(T, '.')
plt.subplot(313)
plt.plot(T, '.')
plt.show()

plt.plot(np.mod(haidx2-13500/2, 13500), T, '.')
plt.show()

plt.subplot(211)
plt.plot(np.mod(haidx1-270/2, 270), np.abs(np.arctan(a/b)), '.')
plt.subplot(212)
plt.plot(np.mod(haidx1-270/2, 270), np.sqrt(a**2+b**2), '.')
plt.show()

solution1 = F[staridx]*T[hauni2]
solution2 = a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1

for i in range(len(ascc)):
    here = staridx == i
    
    ax = plt.subplot(311)
    plt.plot(flux[here], '.')
    plt.plot(solution1[here]*solution2[here], '.')
    
    ax2 = plt.subplot(312, sharex=ax)
    plt.plot(flux[here]/solution1[here], '.')
    plt.plot(solution2[here], '.')
    
    plt.subplot(313, sharex=ax, sharey=ax2)
    plt.plot(flux[here]/(solution1[here]*solution2[here]), '.')
    
    plt.show()
    plt.close()
  
    
