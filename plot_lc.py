#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt

from coordinate_grids import PolarGrid

import os
import glob

path = '/data2/talens/Jul2015/'
date = '20150716'
ascc = '807144'

cameras = ['LPN', 'LPE', 'LPS', 'LPW', 'LPC']
colors = ['b', 'r', 'g', 'y', 'k']

plt.figure(figsize=(16,8))
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)

for camera, color in zip(cameras, colors):

    if os.path.isfile(path+'fLC_'+date+camera+'.hdf5'):
        with h5py.File(path+'fLC_'+date+camera+'.hdf5', 'r') as f:
            try: 
                lc = f['data/'+ascc]
            except: 
                print 'Star %s not found on camera %s'%(ascc, camera)
            else:
                g = h5py.File(path+'Trans0716'+camera+'_pg2700x720.hdf5', 'r')
                bins = g['Data/binnum'].value
                trans = g['Data/trans'].value
                g.close()
                
                pg = PolarGrid(2700, 720)
                trans = pg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
                trans = np.ravel(trans)
                
                vmag = f['header/'+ascc]['vmag']
                ra = f['header/'+ascc]['ra']
                dec = f['header/'+ascc]['dec']
                Nobs = f['header/'+ascc]['nobs'].astype('int')
                ref_jd = np.floor(lc['jdmid'])
                
                ha = np.mod(lc['lst']*15.-np.repeat(ra,Nobs), 360.)
                dec = np.repeat(dec, Nobs)
                
                binnum = pg.find_gridpoint(ha, dec)
                
                ax1.plot(lc['jdmid']-ref_jd, lc['flux0'], '.', label=camera, c=color)
                ax2.plot(lc['jdmid']-ref_jd, trans[binnum], '.', c=color)
                ax3.plot(lc['jdmid']-ref_jd, lc['flux0']/trans[binnum], '.', c=color)
    else:
        print 'No data for camera %s.'%camera
        
plt.sca(ax1)
plt.title('ASCC %s, V = %.2f, RA = %.2f, Dec = %.2f'%(ascc, vmag, ra, dec[0]), size='x-large')
plt.legend(loc=2)
plt.xlim(0.25, 0.75)
plt.xticks(size='large')
plt.yticks(size='large')
plt.xlabel('Time [JD-%i]'%ref_jd[0], size='x-large')
plt.ylabel('Flux', size='x-large')

plt.sca(ax2)
plt.xlim(0.25, 0.75)
plt.xticks(size='large')
plt.yticks(size='large')
plt.xlabel('Time [JD-%i]'%ref_jd[0], size='x-large')
plt.ylabel('Trans', size='x-large')

plt.sca(ax3)
plt.xlim(0.25, 0.75)
plt.xticks(size='large')
plt.yticks(size='large')
plt.xlabel('Time [JD-%i]'%ref_jd[0], size='x-large')
plt.ylabel('cFlux', size='x-large')

plt.tight_layout()
plt.show()
plt.close()

