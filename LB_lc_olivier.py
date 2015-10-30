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

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter1.hdf5', 'r') as f:
    z = f['data/z'].value
    
pg = PolarGrid(13500, 720)
    
with h5py.File('/data2/talens/3mEast/LBtests/sky_15day_iter1.hdf5', 'r') as f:
    s = f['data/s'].value

hg = HealpixGrid(8)
    
with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    
    lijst = [1342743, 437467, 550807, 810093, 725586, 1094380, 224263, 650758, 690250, 1145083, 607190, 542390, 832916, 744576, 647497, 457661, 965558, 386207, 690485, 433284, 869594, 829540]
    for sID in lijst:
        sID = str(sID)
        i, = np.where(ascc == sID)
    
        lc = f['data/' + sID].value
        
        ha = np.mod(lc['lst']*15 - ra[i], 360)
        camtransidx = pg.find_gridpoint(ha, np.repeat(dec[i], len(lc)))
        
        skyidx = hg.find_gridpoint(ra[i], dec[i])
        dayidx = np.floor(lc['jdmid']).astype('int')
        dayidx = dayidx - 2457175
        lstidx = lc['lstidx'].astype('int')
        skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
        
        mag = 25 - 2.5*np.log10(lc['flux0'])
        emag = 2.5/np.log(10)*lc['eflux0']/lc['flux0']
        zsol = z[camtransidx]
        ssol = s[skyidx, skytransidx]
        m = np.nanmean(mag[np.isfinite(mag)])
        
        here = (lc['flux0'] > 0) & (lc['eflux0'] > 0) & (lc['sky'] > 0) & (lc['flag'] < 1)
        flag = np.where(here, 0, 1)

        #mag = mag[here]
        #emag = emag[here]
        
        #zsol = zsol[here]
        #ssol = ssol[here]
        
        #dayidx = dayidx[here]
        #lstidx = lstidx[here]
        
        #jdmid = lc['jdmid'][here]
        #lst = lc['lst'][here]
        
        #a1, a2, niter, chisq, npoints, npars = systematics_dev.trans(dayidx, lstidx, mag - zsol - ssol, emag)
       
        #plt.plot(a2, '.')
        #plt.show()
        
        #plt.figure(figsize=(16,8))
        #ax = plt.subplot(311)
        #plt.plot(mag, '.')
        #plt.plot(zsol + m, '.')
        #plt.subplot(312, sharex=ax, sharey=ax)
        #plt.plot(mag - zsol, '.')
        #plt.plot(ssol + m, '.')
        #plt.subplot(313, sharex=ax, sharey=ax)
        #plt.plot(mag - zsol - ssol, '.')
        #plt.tight_layout()
        #plt.show()
        
        with h5py.File('lightcurves.hdf5') as g:
            
            g.create_dataset('data/'+sID+'/jdmid', data=lc['jdmid'])
            g.create_dataset('data/'+sID+'/mag', data=mag-zsol-ssol)
            g.create_dataset('data/'+sID+'/emag', data=emag)
            g.create_dataset('data/'+sID+'/sky', data=lc['sky'])
            g.create_dataset('data/'+sID+'/flag', data=flag)
            
            
        
