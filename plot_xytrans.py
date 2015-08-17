#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from viridis import viridis

from coordinate_grids import PolarGrid, CartesianGrid, HealpixGrid

with h5py.File('/data2/talens/Jul2015/Trans0716LPS_pg2700x720.hdf5') as f:
    bins = f['binnum'].value
    trans = f['trans'].value

pg = PolarGrid(2700, 720)
trans1 = pg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
trans1 = np.ravel(trans1)

with h5py.File('/data2/talens/Jul2015/Trans0716LPS_cg400x300m50.hdf5') as f:
    bins = f['binnum'].value
    trans = f['trans'].value

cg = CartesianGrid(400, 300, margin=50)
trans2 = cg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
trans2 = np.ravel(trans2)

with h5py.File('/data2/talens/Jul2015/Trans0716LPS_hg512.hdf5') as f:
    bins = f['binnum'].value
    trans = f['trans'].value

hg = HealpixGrid(512)
trans3 = hg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
trans3 = np.ravel(trans3)

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPS.hdf5') as f:
    
    stars = f['table_header/ascc'].value
    vmag = f['table_header/vmag'].value
    dec = f['table_header/dec'].value
    ra = f['table_header/ra'].value
    nobs = f['table_header/nobs'].value
    
    for i in range(len(stars)):

        lc = f['data/'+stars[i]]

        ha = np.mod(lc['lst']*15.-np.repeat(ra[i],nobs[i]), 360.)
        dec_temp = np.repeat(dec[i], nobs[i])

        binnum1 = pg.find_gridpoint(ha, dec_temp)
        binnum2 = cg.find_gridpoint(lc['x'], lc['y'])
        binnum3 = hg.find_gridpoint(ha, dec_temp)
    
        ref = np.floor(lc['jdmid'])
    
        ax = plt.subplot(311)
        plt.plot(lc['jdmid']-ref, lc['flux0'], '.')
        
        plt.xlabel('JD')
        plt.ylabel('Flux')
        
        plt.subplot(312, sharex=ax)
        plt.plot(lc['jdmid']-ref, trans1[binnum1], '.', label='polar')
        plt.plot(lc['jdmid']-ref, trans2[binnum2], '.', label='cartesian')
        plt.plot(lc['jdmid']-ref, trans3[binnum3], '.', label='healpix')
        plt.ylim(0,1)
        plt.legend()
        
        plt.subplot(313, sharex=ax)
        plt.plot(lc['jdmid']-ref, lc['flux0']/trans1[binnum1], '.')
        plt.plot(lc['jdmid']-ref, lc['flux0']/trans2[binnum2], '.')
        plt.plot(lc['jdmid']-ref, lc['flux0']/trans3[binnum3], '.')
        
        plt.tight_layout()
        plt.show()

pg = PolarGrid(13500, 720)
trans = pg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)

plt.subplot(111)
plt.imshow(trans[1:-1,1:-1].T, interpolation='None', origin='lower', aspect='auto', extent=(0,360,-90,90), vmin=0, vmax=1.5, cmap=viridis)
plt.colorbar().set_label('Transmission')
plt.xlabel('HA [deg]')
plt.ylabel('Dec [deg]')

plt.show()
