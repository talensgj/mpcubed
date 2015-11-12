#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

P = 2.21857312
Tp = 2454037.612

def running_median(ind, x, size=50):
    
    result = np.zeros(np.amax(ind)+1)
    for i in np.unique(ind):
        here = (ind > i - 25) & (ind < i+25)
        result[i] = np.median(x[here])
        
    return result
        
def running_mean(ind, x, size=50):
    
    result = np.zeros(np.amax(ind)+1)
    for i in np.unique(ind):
        here = (ind > i - 25) & (ind < i+25)
        result[i] = np.mean(x[here])
        
    return result

def solve_linear_exact(x, y, yerr):
    
    arr = np.zeros((2,2,len(y)))
    aim = np.zeros((2,len(y)))
    
    arr[0,0] = x**2
    arr[1,1] = 1
    arr[1,0] = x
    arr[0,1] = x

    aim[0] = x
    aim[1] = 1
    
    weights = 1/yerr**2
    arr = np.sum(weights*arr, axis=2)
    aim = np.sum(weights*y*aim, axis=1)
    
    sol = np.linalg.solve(arr, aim)
    
    return sol

def solve_quadratic_exact(x, y, yerr):

    arr = np.zeros((3,3,len(y)))
    aim = np.zeros((3,len(y)))
    
    arr[0,0] = x**4
    arr[1,1] = x**2
    arr[2,2] = 1
    arr[1,0] = x**3
    arr[0,1] = x**3
    arr[2,0] = x**2
    arr[0,2] = x**2
    arr[2,1] = x
    arr[1,2] = x
    
    aim[0] = x**2
    aim[1] = x
    aim[2] = 1
    
    weights = 1/yerr**2
    arr = np.sum(weights*arr, axis=2)
    aim = np.sum(weights*y*aim, axis=1)
    
    sol = np.linalg.solve(arr, aim)

    return sol

# Initialize reader and coordinate grids.
f = fLCfile('/data2/talens/3mEast/LBtests/June1.hdf5')
pg = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)
hg = HealpixGrid(8)

# Read transmap.
with h5py.File('/data2/talens/3mEast/LBtests/camip_June1_iter5.hdf5', 'r') as g:
    mz = g['data/magnitudes/m'].value
    idx1 = g['data/camtrans/idx'].value
    z = g['data/camtrans/z'].value
    idx2 = g['data/intrapix/idx'].value
    a = g['data/intrapix/a'].value
    b = g['data/intrapix/b'].value
    c = g['data/intrapix/c'].value
    d = g['data/intrapix/d'].value

z = pg.put_values_on_grid(z, idx1, np.nan)
a = pg2.put_values_on_grid(a, idx2, np.nan)
b = pg2.put_values_on_grid(b, idx2, np.nan)
c = pg2.put_values_on_grid(c, idx2, np.nan)
d = pg2.put_values_on_grid(d, idx2, np.nan)

z = np.ravel(z)
a = np.ravel(a)
b = np.ravel(b)
c = np.ravel(c)
d = np.ravel(d)

# Read skymap.
with h5py.File('/data2/talens/3mEast/LBtests/skyip_June1_iter5.hdf5', 'r') as g:
    sigma1 = g['data/magnitudes/sigma'].value
    
    idx = g['data/skytrans/idx'].value
    lstseq = g['data/skytrans/lstseq'].value
    s = g['data/skytrans/s'].value
    sigma2 = g['data/skytrans/sigma'].value
    
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp

tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = sigma2
sigma2 = tmp

# Read header data.
Mascc, Mra, Mdec, Mnobs, Mvmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
Mnobs = Mnobs.astype('int')

# Build indices for the stars.
Mstaridx = np.arange(len(Mascc))
Mskyidx = hg.find_gridpoint(Mra, Mdec)

for ind in range(len(Mascc)):
        
    if Mascc[ind] not in ['807144']: continue
        
    here = (Mstaridx == ind)
    ascc = Mascc[here]
    ra = Mra[here]
    dec = Mdec[here]
    nobs = Mnobs[here]
    staridx = Mstaridx[here]
    skyidx = Mskyidx[here]

    jdmid, lstidx, lst, flux, eflux, sky, flag, x, y = f.read_data(['jdmid', 'lstidx', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], ascc, nobs)
    lstidx = lstidx.astype('int')

    # Build indices.    
    staridx = np.repeat(staridx, nobs)
    
    dayidx = np.floor(jdmid).astype('int')
    dayidx = dayidx - 2457175 #June1
    #dayidx = dayidx - 2457190
    
    skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
        
    ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
    camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))
    intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec, nobs))
        
    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    jdmid = jdmid[here]
    lst = lst[here]
    flux = flux[here]
    eflux = eflux[here]
    x = x[here]
    y = y[here]
    
    lstidx = lstidx[here]
    dayidx = dayidx[here]
    staridx = staridx[here]
    skytransidx = skytransidx[here]
    camtransidx = camtransidx[here]
    intrapixidx = intrapixidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    mag = mag - mz[ind]

    sol = z[camtransidx] + a[intrapixidx]*np.sin(2*np.pi*x) + b[intrapixidx]*np.cos(2*np.pi*x) + c[intrapixidx]*np.sin(2*np.pi*y) + d[intrapixidx]*np.cos(2*np.pi*y)
    mag = mag - sol

    sol = s[skyidx, skytransidx]
    mag = mag - sol

    a1, a2, niter, chisq, npoints, npars = systematics_dev.trans(dayidx, lstidx, mag, emag)
    option1a = a1[dayidx] + a2[lstidx]
    
    a1, a2, niter, chisq, npoints, npars = systematics_dev.trans(dayidx, lstidx, mag, np.sqrt(emag**2 + sigma1[ind]**2 + sigma2[skyidx, skytransidx]**2))
    option1b = a1[dayidx] + a2[lstidx]
   
    plt.figure(figsize=(16,8))
   
    ax = plt.subplot(211)
    plt.plot(option1a, '.')
    plt.plot(option1b, '.')
    
    plt.subplot(212, sharex=ax, sharey=ax)
    plt.plot(option1a - option1b, '.')
   
    plt.show()
    
    a1, a2, niter, chisq, npoints, npars = systematics_dev.trans(dayidx, lstidx//50, mag, emag)
    option2a = a1[dayidx] + a2[lstidx//50]
   
    a1, a2, niter, chisq, npoints, npars = systematics_dev.trans(dayidx, lstidx//50, mag, np.sqrt(emag**2 + sigma1[ind]**2 + sigma2[skyidx, skytransidx]**2))
    option2b = a1[dayidx] + a2[lstidx//50]
   
    plt.figure(figsize=(16,8))
   
    ax = plt.subplot(211)
    plt.plot(option2a, '.')
    plt.plot(option2b, '.')
    
    plt.subplot(212, sharex=ax, sharey=ax)
    plt.plot(option2a - option2b, '.')
   
    plt.show()
   
    # Fit a polynomial.
    lst = lst - (np.amax(lst)+np.amin(lst))/2

    sol = solve_quadratic_exact(lst, mag, emag)
    
    a = sol[0]
    b = sol[1]
    c = sol[2]*np.ones(15)
    weights = 1/emag**2
    for niter in range(50):
        c = np.bincount(dayidx, weights*(mag-a*lst**2-b*lst))/np.bincount(dayidx, weights)
        b = np.sum(weights*(mag-a*lst**2-c[dayidx])*lst)/np.sum(weights*lst**2)
        a = np.sum(weights*(mag-b*lst-c[dayidx])*lst**2)/np.sum(weights*lst**4)

    option3a = a*lst**2 + b*lst + c[dayidx]
    
    a = sol[0]
    b = sol[1]
    c = sol[2]*np.ones(15)
    weights = 1/(emag**2 + sigma1[ind]**2 + sigma2[skyidx, skytransidx]**2)
    for niter in range(50):
        c = np.bincount(dayidx, weights*(mag-a*lst**2-b*lst))/np.bincount(dayidx, weights)
        b = np.sum(weights*(mag-a*lst**2-c[dayidx])*lst)/np.sum(weights*lst**2)
        a = np.sum(weights*(mag-b*lst-c[dayidx])*lst**2)/np.sum(weights*lst**4)

    option3b = a*lst**2 + b*lst + c[dayidx]
    
    plt.figure(figsize=(16,8))
    
    ax = plt.subplot(211)
    plt.plot(option3a, '.')
    plt.plot(option3b, '.')
    
    plt.subplot(212, sharex=ax, sharey=ax)
    plt.plot(option3a - option3b, '.')
   
    plt.show()
    
    result = np.zeros(13500)
    for i in range(50):
        norm = index_statistics(dayidx, mag-result[lstidx], statistic='median', keeplength=True)
        #result = running_median(lstidx, mag-norm)
        result = running_mean(lstidx, mag-norm)
        
    option4 = norm + result[lstidx]
    
    plt.figure(figsize=(16,8))
    
    plt.subplot(211)
    plt.plot(option1a, '.')
    plt.plot(option2a, '.')
    plt.plot(option3a, '.')
    plt.plot(option4, '.')
    
    plt.subplot(212)
    plt.plot(option1a-option3a, '.')
    plt.plot(option2a-option3a, '.')
    plt.plot(option4-option3a, '.')
    
    plt.show()
    
    mag = mag - option4
    
    binidx = np.ravel_multi_index((dayidx, lstidx//50), (15, 270))
    
    count = index_statistics(binidx, None, statistic='count')
    bin_mag = index_statistics(binidx, mag, statistic='mean')
    bin_emag = index_statistics(binidx, mag, statistic='std')/np.sqrt(count)
    bin_jdmid = index_statistics(binidx, jdmid, statistic='mean')
    
    here = count == 50
    bin_mag = bin_mag[here]
    bin_emag = bin_emag[here]
    bin_jdmid = bin_jdmid[here]
    
    phase = (bin_jdmid - Tp)/P
    phase = np.mod(phase+.5, 1)-.5
    
    plt.figure(figsize=(16,8))
    
    plt.subplot(211)
    plt.plot(bin_mag, '.')
    plt.ylim(.1, -.1)
    plt.xlabel('Time')
    plt.ylabel('Magnitudes')
    
    plt.subplot(212)
    plt.plot(phase, bin_mag, '.')
    plt.xlim(-.5, .5)
    plt.ylim(.1, -.1)
    plt.xlabel('Phase')
    plt.ylabel('Magnitudes')
    
    plt.tight_layout()
    plt.show()
