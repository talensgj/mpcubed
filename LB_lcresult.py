#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev
from core import sysrem

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

pg1 = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)

with h5py.File('/data2/talens/3mEast/LBtests/camip_June1_iter5.hdf5', 'r') as f:
    
    idx1 = f['data/camtrans/idx'].value
    z = f['data/camtrans/z'].value
    
    idx2 = f['data/intrapix/idx'].value
    a = f['data/intrapix/a'].value
    b = f['data/intrapix/b'].value
    c = f['data/intrapix/c'].value
    d = f['data/intrapix/d'].value
    
z = pg1.put_values_on_grid(z, idx1, np.nan).ravel()
a = pg2.put_values_on_grid(a, idx2, np.nan).ravel()
b = pg2.put_values_on_grid(b, idx2, np.nan).ravel()
c = pg2.put_values_on_grid(c, idx2, np.nan).ravel()
d = pg2.put_values_on_grid(d, idx2, np.nan).ravel()
    
hg = HealpixGrid(8)
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_June1_iter5.hdf5', 'r') as f:
    
    m = f['data/magnitudes/m'].value
    sigma1 = f['data/magnitudes/sigma'].value
    
    idx = f['data/skytrans/idx'].value
    lstseq = f['data/skytrans/lstseq'].value
    s = f['data/skytrans/s'].value
    sigma2 = f['data/skytrans/sigma'].value
    
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp
    
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = sigma2
sigma2 = tmp    
    
f = fLCfile('/data2/talens/3mEast/LBtests/June1.hdf5')
    
ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
nobs = nobs.astype('int')

skyidx = hg.find_gridpoint(ra, dec)

for i in range(len(ascc)):

    if ascc[i] not in ['265281']: continue

    jdmid, lstidx, lst, flux, eflux, sky, flag, x, y = f.read_data(['jdmid', 'lstidx', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], [ascc[i]], [nobs[i]])
    lstidx = lstidx.astype('int')

    # Build indices.    
    staridx = np.repeat(i, nobs[i])
    
    dayidx = np.floor(jdmid).astype('int')
    dayidx = dayidx - 2457175
    print np.unique(dayidx)
    skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
    binidx = np.ravel_multi_index((dayidx, lstidx//50), (15, 270))
        
    ha = np.mod(lst*15 - ra[i], 360.)
    camtransidx = pg1.find_gridpoint(ha, np.repeat(dec[i], nobs[i]))
    intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec[i], nobs[i]))
    
    # Flag bad data.
    #here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    #jdmid = jdmid[here]
    #flux = flux[here]
    #eflux = eflux[here]
    #x = x[here]
    #y = y[here]

    #staridx = staridx[here]
    #skytransidx = skytransidx[here]
    #camtransidx = camtransidx[here]
    #intrapixidx = intrapixidx[here]

    #dayidx = dayidx[here]
    #lstidx = lstidx[here]
    #binidx = binidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    sol1 = m[i] + z[camtransidx] + a[intrapixidx]*np.sin(2*np.pi*x) + b[intrapixidx]*np.cos(2*np.pi*x) + c[intrapixidx]*np.sin(2*np.pi*y) + d[intrapixidx]*np.cos(2*np.pi*y)
    sol2 = s[skyidx[i], skytransidx]

    a1, a2, niter, chisq, npoints, npars = systematics_dev.trans(dayidx, lstidx, mag-sol1-sol2, np.sqrt(emag+sigma1[i]**2+sigma2[skyidx[i], skytransidx]**2))
    trend = a2[lstidx] + a1[dayidx]

    plt.figure(figsize=(16,8))

    plt.subplot(311)
    plt.plot(mag, '.')
    plt.plot(sol1, '.')
    plt.ylim(np.mean(mag)+1, np.mean(mag)-1)
    plt.ylabel('Magnitude')
    
    plt.subplot(312)
    plt.plot(mag - sol1, '.')
    plt.plot(sol2, '.')
    plt.ylim(1, -1)
    plt.ylabel('Magnitude')
    
    plt.subplot(313)
    plt.plot(mag - sol1 - sol2, '.')
    plt.plot(trend, '.')
    plt.ylim(1, -1)
    plt.xlabel('Time')
    plt.ylabel('Magnitude')
    
    plt.tight_layout()
    plt.show()
    
    count = index_statistics(binidx, None, statistic='count')
    bin_jdmid = index_statistics(binidx, jdmid, statistic='mean')
    bin_mag = index_statistics(binidx, (mag-sol1-sol2-trend), statistic='mean')
    bin_emag = index_statistics(binidx, (mag-sol1-sol2-trend), statistic='std')/np.sqrt(count)
    
    bin_jdmid = bin_jdmid[count==50]
    bin_mag = bin_mag[count==50]
    bin_emag = bin_emag[count==50]
    
    P = 2.21857312
    Tp = 2454037.612
    phase = (bin_jdmid - Tp)/P
    phase = np.mod(phase+.5, 1)-.5
    
    args, = np.where(np.abs(phase)<1e-2)
    
    plt.figure(figsize=(16,8))
    
    plt.subplot(211)
    plt.title('HD189733', size='xx-large')
    plt.errorbar(bin_jdmid, bin_mag, yerr=bin_emag, fmt='o')
    
    for i in args:
        plt.axvline(bin_jdmid[i])
    
    plt.xlim(np.amin(bin_jdmid)-.01*np.ptp(bin_jdmid), np.amax(bin_jdmid)+.01*np.ptp(bin_jdmid))
    plt.ylim(.05,-.05)
    plt.ylabel('Magnitude')
    
    plt.subplot(212)
    plt.errorbar(phase, bin_mag, yerr=bin_emag, fmt='o')
    plt.xlim(-.5,.5)
    plt.ylim(.05,-.05)
    plt.xlabel('Time [JD]')
    plt.ylabel('Magnitude')
    
    plt.tight_layout()
    plt.show()
    
    
    
    
    
    
    
