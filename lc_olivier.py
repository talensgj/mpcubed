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

# Read camera transmission.
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

# Read skytransmission
hg = HealpixGrid(8)
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_June1_iter5.hdf5', 'r') as f:
    
    idx = f['data/skytrans/idx'].value
    lstseq = f['data/skytrans/lstseq'].value
    s = f['data/skytrans/s'].value
    
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp
    
with h5py.File('/data2/talens/3mEast/LBtests/June1.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    nobs = f['header_table/nobs'].value
    vmag = f['header_table/vmag'].value
    
    lijst = np.genfromtxt('/home/talens/star_list.dat')
    lijst = lijst.astype('int')
    lijst = lijst.astype('str')
    
    #lijst = [1342743, 437467, 550807, 810093, 725586, 1094380, 224263, 650758, 690250, 1145083, 607190, 542390, 832916, 744576, 647497, 457661, 965558, 386207, 690485, 433284, 869594, 829540]
    for sID in lijst:
        sID = str(sID)
        i, = np.where(ascc == sID)
    
        lc = f['data/' + sID].value
        
        # Build indices.
        ha = np.mod(lc['lst']*15 - ra[i], 360)
        
        camtransidx = pg1.find_gridpoint(ha, np.repeat(dec[i], len(lc)))
        intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec[i], len(lc)))
        
        skyidx = hg.find_gridpoint(ra[i], dec[i])
        
        dayidx = np.floor(lc['jdmid']).astype('int')
        dayidx = dayidx - 2457175
        
        lstidx = lc['lstidx'].astype('int')
        
        skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
        
        # Compute corrected magnitude.
        mag = 25 - 2.5*np.log10(lc['flux0'])
        emag = 2.5/np.log(10)*lc['eflux0']/lc['flux0']
        
        cam = z[camtransidx]
        ipx = a[intrapixidx]*np.sin(2*np.pi*lc['x']) + b[intrapixidx]*np.cos(2*np.pi*lc['x']) + c[intrapixidx]*np.sin(2*np.pi*lc['y']) + d[intrapixidx]*np.cos(2*np.pi*lc['y'])
        sky = s[skyidx, skytransidx]
        
        mag = mag - cam - ipx - sky
        
        weights = 1/emag**2
        t = np.bincount(lstidx, weights*mag)/np.bincount(lstidx, weights)
        
        # Create simple flags.
        here = (lc['flux0'] > 0) & (lc['eflux0'] > 0) & (lc['sky'] > 0) & (lc['flag'] < 1)
        flag = np.where(here, 0, 1)
        
        #ax = plt.subplot(211)
        #plt.plot(mag, '.')
        #plt.subplot(212, sharex=ax)
        #plt.plot(t[lstidx], '.')
        #plt.show()
        
        mag = mag - t[lstidx]
        
        with h5py.File('detrended_lightcurves.hdf5') as g:
            
            g.create_dataset('data/'+sID+'/jdmid', data=lc['jdmid'])
            g.create_dataset('data/'+sID+'/mag', data=mag)
            g.create_dataset('data/'+sID+'/emag', data=emag)
            g.create_dataset('data/'+sID+'/sky', data=lc['sky'])
            g.create_dataset('data/'+sID+'/flag', data=flag)
    
    here = np.in1d(ascc, lijst)
        
    with h5py.File('detrended_lightcurves.hdf5') as g:
        
        g.create_dataset('header_table/ascc', data=ascc[here])
        g.create_dataset('header_table/ra', data=ra[here])
        g.create_dataset('header_table/dec', data=dec[here])
        g.create_dataset('header_table/vmag', data=vmag[here])
        g.create_dataset('header_table/nobs', data=nobs[here])
        
