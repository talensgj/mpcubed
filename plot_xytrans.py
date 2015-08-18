#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from viridis import viridis

from coordinate_grids import PolarGrid, CartesianGrid, HealpixGrid

with h5py.File('/data2/talens/Jul2015/Trans0716LPC_pg2700x720.hdf5') as f:
    bins = f['Data/binnum'].value
    trans = f['Data/trans'].value
    count = f['Data/count'].value

pg = PolarGrid(2700, 720)

trans = pg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
trans = trans[1:-1,1:-1]

wrap = True
if wrap:
    trans = np.roll(trans, 1350, axis=0)

ax = plt.subplot(111)
plt.imshow(trans.T, interpolation='None', origin='lower', aspect='auto', extent=(0,360,-90,90), vmin=0, vmax=1, cmap=viridis)
plt.colorbar().set_label('Transmission')
if wrap:
    xticks = ax.get_xticks()
    plt.xticks(xticks, ['%i'%tick for tick in (xticks-180)%360])
plt.xlabel('HA [deg]')
plt.ylabel('Dec [deg]')

plt.show()

count = pg.put_values_on_grid(count, ind=bins, fill_value=np.nan)
count = count[1:-1,1:-1]

wrap = True
if wrap:
    count = np.roll(count, 1350, axis=0)

plt.subplot(111)
plt.imshow(count.T, interpolation='None', origin='lower', aspect='auto', extent=(0,360,-90,90), vmin=np.nanmin(count), vmax=np.nanmax(count), cmap=viridis)
plt.colorbar().set_label('# points')
if wrap:
    plt.xticks([0,60,120,180,240,300,360], ['180', '240', '300', '0', '60', '120', '180'])
else:
    plt.xticks([0,60,120,180,240,300,360])
plt.yticks([-90, -60, -30, 0, 30, 60, 90])
plt.xlabel('HA [deg]')
plt.ylabel('Dec [deg]')

plt.show()

#jd_ref = np.floor(lc['jdmid'])
        
#ax = plt.subplot(311)
#plt.plot(lc['jdmid']-jd_ref, lc['flux0'], '.')
#plt.xlabel('JD')
#plt.ylabel('Flux')

#plt.subplot(312, sharex=ax)
#plt.plot(lc['jdmid']-jd_ref, trans[binnum], '.')
#plt.ylim(0,1)
#plt.xlabel('JD')
#plt.ylabel('Trans')

#plt.subplot(313, sharex=ax)
#plt.plot(lc['jdmid']-jd_ref, lc['flux0']/trans[binnum], '.')
#plt.xlim(.25,.75)
#plt.xlabel('JD')
#plt.ylabel('Trans')

#plt.tight_layout()
#plt.show()
#plt.close()
