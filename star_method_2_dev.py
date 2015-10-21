#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from core import sysrem
from core import coarse_decor
from core.coordinate_grids import HealpixGrid

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

P = 2.21857312
Tp = 2454037.612

filelist = glob.glob('/data2/talens/3mEast/fLC_201506??LPE.hdf5')
filelist = np.sort(filelist)

mascc = np.array([])
mnobs = np.array([])

mjdmid = np.array([])
mlstidx = np.array([])
mflux0 = np.array([])
meflux0 = np.array([])
mx = np.array([])
my = np.array([])
mflags = np.array([])

for filename in filelist:
    print filename
    with h5py.File(filename, 'r') as f:
        
        ascc = f['header_table/ascc'].value
        ra = f['header_table/ra'].value
        dec = f['header_table/dec'].value
        nobs = f['header_table/nobs'].value
        
        hg = HealpixGrid(8)
        skyidx = hg.find_gridpoint(ra, dec)
        here = (skyidx == 266)
        
        ascc = ascc[here]
        nobs = nobs[here]
        mascc = np.append(mascc, ascc)
        mnobs = np.append(mnobs, nobs)
        
        select = np.append(0, np.cumsum(nobs)).astype('int')
        jdmid = np.zeros(np.sum(nobs))
        lstidx = np.zeros(np.sum(nobs))
        flux0 = np.zeros(np.sum(nobs))
        eflux0 = np.zeros(np.sum(nobs))
        x = np.zeros(np.sum(nobs))
        y = np.zeros(np.sum(nobs))
        flags = np.zeros(np.sum(nobs))
        
        for i in range(len(ascc)):
        
            lc = f['data/'+ascc[i]].value
           
            jdmid[select[i]:select[i+1]] = lc['jdmid']
            lstidx[select[i]:select[i+1]] = lc['lstidx']
            flux0[select[i]:select[i+1]] = lc['flux0']
            eflux0[select[i]:select[i+1]] = lc['eflux0']
            x[select[i]:select[i+1]] =  lc['x']
            y[select[i]:select[i+1]] =  lc['y']
            flags[select[i]:select[i+1]] = lc['flag']
 
    mjdmid = np.append(mjdmid, jdmid)
    mlstidx = np.append(mlstidx, lstidx)
    mflux0 = np.append(mflux0, flux0)
    meflux0 = np.append(meflux0, eflux0)
    mx = np.append(mx, x)
    my = np.append(my, y)
    mflags = np.append(mflags, flags)
 
ascc = mascc
nobs = mnobs.astype('int')

ascc = np.repeat(ascc, nobs)

jdmid = mjdmid
lstidx = mlstidx
flux0 = mflux0
eflux0 = meflux0
x = mx
y = my
flags = mflags
            
lstidx = lstidx.astype('int')
dayidx = np.floor(jdmid)
dayidx = dayidx - np.amin(dayidx)
dayidx = dayidx.astype('int')
ascc, staridx = np.unique(ascc, return_inverse = True)

print ascc

here = (flux0 > 0) & (eflux0 > 0) & (flags < 1)
jdmid = jdmid[here]
dayidx = dayidx[here]
lstidx = lstidx[here]
staridx = staridx[here]
flux0 = flux0[here]
eflux0 = eflux0[here]
x = x[here]
y = y[here]

mag0 = -2.5*np.log10(flux0)
emag0 = 2.5/np.log(10.)*eflux0/flux0

#midx = np.ravel_multi_index((dayidx, staridx), (31, len(ascc)))
midx = staridx
tidx = np.ravel_multi_index((staridx, lstidx), (len(ascc), 13500))
sidx = np.ravel_multi_index((dayidx, lstidx), (31, 13500))

weights = 1/(emag0)**2
t = np.zeros(np.amax(tidx)+1)
s = np.zeros(np.amax(sidx)+1)

for i in range(15):
    m = np.bincount(midx, weights*(mag0 - t[tidx] - s[sidx]))/np.bincount(midx, weights)
    t = np.bincount(tidx, weights*(mag0 - m[midx] - s[sidx]))/np.bincount(tidx, weights)
    s = np.bincount(sidx, weights*(mag0 - m[midx] - t[tidx]))/np.bincount(sidx, weights)

    sol = m[midx] + t[tidx] + s[sidx]
    res = mag0 - sol
    sigma1 = coarse_decor.find_sigma(midx, res, emag0)
    sigma2 = coarse_decor.find_sigma(dayidx, res, emag0)
    weights = 1/(mag0**2 + (sigma1**2)[midx] + (sigma2**2)[dayidx])

    if (i > 0):
        print np.nanmax(np.abs(m-m_old)), np.nanmax(np.abs(t-t_old)), np.nanmax(np.abs(s-s_old))

    m_old = np.copy(m)
    t_old = np.copy(t)
    s_old = np.copy(s)

plt.plot(t, '.')
plt.show()    

plt.plot(s, '.')
plt.show()

sol1 = m[midx]
sol2 = t[tidx]
sol3 = s[sidx]

for i in range(len(ascc)):
    here = (staridx == i)

    binidx = np.ravel_multi_index((dayidx[here], lstidx[here]//50), (31, 270))
    count = np.bincount(binidx)
    bin_mag0 = np.bincount(binidx, mag0[here] - sol1[here] - sol2[here] - sol3[here])/count
    bin_jdmid = np.bincount(binidx, jdmid[here])/count

    here = (count > 49)
    bin_mag0 = bin_mag0[here]
    bin_jdmid = bin_jdmid[here]
    
    plt.title(ascc[i])
    plt.plot(bin_mag0, '.')
    plt.ylim(-.1, .1)
    plt.show()
    plt.close()





