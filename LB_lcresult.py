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

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter3.hdf5', 'r') as f:
    z = f['data/z'].value
    a = f['data/a'].value
    b = f['data/b'].value
    c = f['data/c'].value
    d = f['data/d'].value
    
pg1 = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter3.hdf5', 'r') as f:
    m = f['data/m'].value
    s = f['data/s'].value

hg = HealpixGrid(8)
    
with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    
    for i in range(0, len(ascc), 50):
        
        lc = f['data/'+ascc[i]].value
    
        ha = np.mod(lc['lst']*15 - ra[i], 360)
        camtransidx = pg1.find_gridpoint(ha, np.repeat(dec[i], len(lc)))
        intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec[i], len(lc)))
        skyidx = hg.find_gridpoint(ra[i], dec[i])
        dayidx = np.floor(lc['jdmid']).astype('int')
        dayidx = dayidx - 2457175
        lstidx = lc['lstidx'].astype('int')
        skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))

        mag = 25 - 2.5*np.log10(lc['flux0'])
        emag = 2.5/np.log(10)*lc['flux0']/lc['eflux0']
        zsol = z[camtransidx]
        ssol = s[skyidx, skytransidx]
        ipx = a[intrapixidx]*np.sin(2*np.pi*lc['x']) + b[intrapixidx]*np.cos(2*np.pi*lc['x']) + c[intrapixidx]*np.sin(2*np.pi*lc['y']) + d[intrapixidx]*np.cos(2*np.pi*lc['y']) 
        
        here = np.isfinite(mag)
        mag = mag[here]
        emag = emag[here]
        zsol = zsol[here]
        ssol = ssol[here]
        dayidx = dayidx[here]
        lstidx = lstidx[here]
        jdmid = lc['jdmid'][here]
        lst = lc['lst'][here]
        ipx = ipx[here]
        
        #a = 0
        #for i in range(1):
            #here = dayidx != 1
            #a1, a2, niter, chisq, npoints, npars = systematics_dev.trans(dayidx[here], lstidx[here], (mag - z - s - ipx)[here], emag[here])
            #print a1
            #plt.plot(a2, '.')
            #plt.show()
            
            #a = a + a2[lstidx]
        
        plt.figure(figsize=(16,8))
        ax = plt.subplot(311)
        plt.title(ascc[i])
        plt.plot(mag, '.')
        plt.plot(zsol + m[i] + ipx, '.')
        plt.subplot(312, sharex=ax, sharey=ax)
        plt.plot(mag - zsol - ipx, '.')
        plt.plot(ssol + m[i], '.')
        plt.subplot(313, sharex=ax, sharey=ax)
        plt.plot(mag - zsol - ssol - ipx, '.')
        #plt.plot(a + m, '.')
        plt.ylim(m[i]-1, m[i]+1)
        plt.tight_layout()
        plt.show()
        
        #binidx = np.ravel_multi_index((dayidx, lstidx//50), (15, 270))
        
        #count = np.bincount(binidx)
        #bin_mag = np.bincount(binidx, mag - z - s - a)/count
        #bin_jdmid = np.bincount(binidx, jdmid)/count
        #bin_lst = np.bincount(binidx, lst)/count
        #bin_mag = bin_mag[count==50]
        #bin_jdmid = bin_jdmid[count==50]
        #bin_lst = bin_lst[count==50]
        
        #phase = (bin_jdmid - 2454037.612)/2.21857312
        #phase = np.mod(phase+.5, 1)-.5
        
        #plt.figure(figsize=(16,8))
        #ax = plt.subplot(311)
        #plt.scatter(bin_lst, bin_mag)
        #plt.xlabel('LST')
        #plt.ylabel('magnitude')
        #plt.subplot(312, sharey=ax)
        #plt.scatter(bin_jdmid, bin_mag)
        #plt.xlabel('JD')
        #plt.ylabel('magnitude')
        #plt.subplot(313, sharey=ax)
        #plt.scatter(phase, bin_mag)
        #plt.xlabel('Phase')
        #plt.ylabel('magnitude')
        #plt.tight_layout()
        #plt.show()
    
    
    
    
