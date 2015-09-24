#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from coordinate_grids import PolarGrid
from coarse_decor import coarse_positions

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/Jul2015/coarsecam.hdf5', 'r') as f:
    
    decidx = f['data'].keys()
    camtrans = np.full((722, 13502), fill_value=np.nan)
    pointcount_cam = np.full((722, 13502), fill_value=np.nan)
    a = np.full((722, 272), fill_value=np.nan)
    b = np.full((722, 272), fill_value=np.nan)
    for ind in decidx:
        
        cam = f['data/'+ind]
        haidx = cam['haidx_cam'].value
        camtrans[ind, haidx] = cam['camtrans']
        pointcount_cam[ind, haidx] = cam['pointcount_cam']
        haidx = cam['haidx_ipx'].value
        a[ind, haidx] = cam['a']
        b[ind, haidx] = cam['b']

plt.imshow(camtrans, aspect='auto')
plt.show()

camtrans = camtrans.T
pointcount_cam = pointcount_cam.T
a = a.T
b = b.T

filename = '/data2/talens/Jul2015/fLC_20150710LPC.hdf5'

# Polar grid instance.
pgcam = PolarGrid(13500, 720)
pgipx = PolarGrid(270, 720)

with h5py.File(filename, 'r') as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    
    decidx = pgcam.find_decidx(dec)
    
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]]
        
        ha = np.mod(lc['lst']*15.-ra[i], 360)
        haidx_cam = pgcam.find_raidx(ha) 
        haidx_ipx = pgipx.find_raidx(ha)
        
        camtrans0 = camtrans[haidx_cam, decidx[i]]
        cflux0 = -2.5*np.log10(lc['flux0'])-camtrans0
        
        intrapix0 = a[haidx_ipx, decidx[i]]*np.sin(2*np.pi*lc['y'])+b[haidx_ipx, decidx[i]]*np.cos(2*np.pi*lc['y'])
        ipcflux0 = cflux0-intrapix0

        flags = np.where(np.isnan(camtrans0), 1, 0)
        flags = flags + np.where(pointcount_cam[haidx_cam, decidx[i]]<=5, 2, 0)
        flags = flags + np.where(np.isnan(intrapix0), 4, 0)
        
        with h5py.File('/data2/talens/Jul2015/coarsered.hdf5') as g:
            grp = g.create_group('data/'+ascc[i])
            grp.create_dataset('camtrans0', data=camtrans0)
            grp.create_dataset('intrapix0', data=intrapix0)
            grp.create_dataset('cflux0', data=cflux0)
            grp.create_dataset('ipcflux0', data=ipcflux0)
            grp.create_dataset('flags', data=flags)

