#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

from coordinate_grids import PolarGrid
from index_functions import index_statistics

from sysrem import sysrem


def bin_intrapixel(ind1, ind2, y, flux, eflux, maxiter=500, eps=1e-3):
    
    cs = np.cos(2*np.pi*y)
    sn = np.sin(2*np.pi*y)
    w = 1/eflux**2
    
    a = np.zeros(np.amax(ind2)+1)
    b = np.zeros(np.amax(ind2)+1)
    
    niter = 0
    while (niter < maxiter):
        
        print niter
        
        f = np.bincount(ind1, flux*(a[ind2]*sn+b[ind2]*cs+1)*w)/np.bincount(ind1, (a[ind2]*sn+b[ind2]*cs+1)**2*w)
        a = np.bincount(ind2, (flux-f[ind1]*(b[ind2]*cs+1))*f[ind1]*sn*w)/np.bincount(ind2, (f[ind1]*sn)**2*w)
        b = np.bincount(ind2, (flux-f[ind1]*(a[ind2]*sn+1))*f[ind1]*cs*w)/np.bincount(ind2, (f[ind1]*cs)**2*w)
        
        if niter == 0:
            a_old = np.copy(a) 
            b_old = np.copy(b)
            
        else:
            
            crita = np.nanmax((a-a_old)/a_old)
            critb = np.nanmax((b-b_old)/b_old)
            
            print crita, critb
            
            if (crita < eps) & (critb < eps):
                break
                
            a_old = np.copy(a)
            b_old = np.copy(b)
    
        niter += 1
    
    return f, a, b
    
    
def trans_intrapixel(ind1, ind2, ind3, y, flux, eflux, maxiter=500, eps=1e-3):
    
    cs = np.cos(2*np.pi*y)
    sn = np.sin(2*np.pi*y)
    w = 1/eflux**2
    
    T = np.ones(np.amax(ind2)+1)
    
    a = np.zeros(np.amax(ind3)+1)
    b = np.zeros(np.amax(ind3)+1)
    
    niter = 0
    while (niter < maxiter):
        
        print niter
        
        F = np.bincount(ind1, flux*T[ind2]*(a[ind3]*sn+b[ind3]*cs+1)*w)/np.bincount(ind1, (T**2)[ind2]*(a[ind3]*sn+b[ind3]*cs+1)**2*w)
        T = np.bincount(ind2, flux*F[ind1]*(a[ind3]*sn+b[ind3]*cs+1)*w)/np.bincount(ind2, (F**2)[ind1]*(a[ind3]*sn+b[ind3]*cs+1)**2*w)
        
        a = np.bincount(ind3, (flux-F[ind1]*T[ind2]*(b[ind3]*cs+1))*F[ind1]*T[ind2]*sn*w)/np.bincount(ind3, (F[ind1]*T[ind2]*sn)**2*w)
        b = np.bincount(ind3, (flux-F[ind1]*T[ind2]*(a[ind3]*sn+1))*F[ind1]*T[ind2]*cs*w)/np.bincount(ind3, (F[ind1]*T[ind2]*cs)**2*w)
    
        if niter == 0:
            T_old = np.copy(T)
            a_old = np.copy(a) 
            b_old = np.copy(b)
            
        else:
            
            critT = np.nanmax((T-T_old)/T_old)
            crita = np.nanmax((a-a_old)/a_old)
            critb = np.nanmax((b-b_old)/b_old)
            
            print critT, crita, critb
            
            if (critT < eps) & (crita < eps) & (critb < eps):
                break
            
            if np.isnan(critT):
                break
            
            T_old = np.copy(T)
            a_old = np.copy(a)
            b_old = np.copy(b)
    
        niter += 1
            
    return F, T, a, b
    
fLC = '/data2/talens/Jul2015/fLC_20150716LPC.hdf5'
red = '/data2/talens/Jul2015/red_20150716LPC.hdf5'

with h5py.File(fLC, 'r') as f:
    
    ascc1 = f['table_header/ascc'].value
    ra1 = f['table_header/ra'].value
    dec1 = f['table_header/dec'].value
    nobs1 = f['table_header/nobs'].value.astype('int')

pg1 = PolarGrid(270, 720)
decidx, decuni = pg1.find_decidx(dec1, compact=True)

array1 = np.full((272, 722), fill_value=np.nan)
array2 = np.full((272, 722), fill_value=np.nan)
for ind in range(len(decidx)):

    ascc = ascc1[decuni==ind]
    ra = ra1[decuni==ind]
    dec = dec1[decuni==ind]
    nobs = nobs1[decuni==ind]

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
    haidx1, hauni1 = pg1.find_raidx(ha, compact=True)
    staridx = np.repeat(np.arange(len(ascc)), nobs)
    idx = np.ravel_multi_index((staridx, hauni1), (len(ascc), len(haidx1)))

    pg = PolarGrid(13500, 720)
    haidx2, hauni2 = pg.find_raidx(ha, compact=True)

    #f, a, b = bin_intrapixel(idx, hauni1, y, flux, eflux)
    here = (flux>0)&(eflux>0)
    F, T, a1, b1 = trans_intrapixel(staridx[here], hauni2[here], hauni1[here], y[here], flux[here], eflux[here])
    
    array1[haidx1, [ind]*len(a1)] = np.sqrt(a1**2+b1**2)
    array2[haidx1, [ind]*len(a1)] = np.arctan(a1/b1)%np.pi


plt.imshow(np.roll(array1[1:-1,1:-1], 135, axis=0).T, aspect='auto', vmin=0, vmax=.1)
plt.colorbar()
plt.show()
    
plt.imshow(np.roll(array2[1:-1,1:-1], 135, axis=0).T, aspect='auto', vmin=0, vmax=np.pi)
plt.colorbar()
plt.show()
    
exit()

plt.subplot(211)
plt.plot(np.mod(haidx1-135, 270), np.arctan(a/b)%np.pi, '.')
plt.plot(np.mod(haidx1-135, 270), np.arctan(a1/b1)%np.pi, '.')
plt.subplot(212)
plt.plot(np.mod(haidx1-135, 270), np.sqrt(a**2+b**2), '.')
plt.plot(np.mod(haidx1-135, 270), np.sqrt(a1**2+b1**2), '.')
plt.show()
solution = f[idx]*(a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1)
intrapix = (a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1)
solution1 = F[staridx]*T[hauni2]*(a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1)
camtrans1 = F[staridx]*T[hauni2]
intrapix1 = (a[hauni1]*np.sin(2*np.pi*y)+b[hauni1]*np.cos(2*np.pi*y)+1)

for i in range(len(ascc)):
    
    here = (staridx == i)
    
    ax = plt.subplot(311)
    plt.plot(flux[here], '.', c='k')
    
    plt.plot(solution[here], '.', c='r')
    plt.plot(solution1[here], '.', c='g')
    
    ax1 = plt.subplot(312, sharex=ax)
    plt.plot(flux[here]/solution[here], '.', c='r')
    
    plt.subplot(313, sharex=ax, sharey=ax1)
    plt.plot(flux[here]/camtrans1[here], '.', c='k')
    plt.plot(flux[here]/solution1[here], '.', c='g')
    
    plt.show()
    plt.close()
