#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 
import matplotlib.gridspec as gridspec

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def delete_allnan(array):
    
    ind1, = np.where(np.all(np.isnan(array), axis=0))
    ind2, = np.where(np.all(np.isnan(array), axis=1))
    
    array = np.delete(array, ind1, axis=1)
    array = np.delete(array, ind2, axis=0)
    
    return array

with h5py.File('/data2/talens/3mEast/LBtests/15day.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value

#with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter5_weights.hdf5', 'r') as f:
    #mz = f['data/m'].value
    #z = f['data/z'].value
    #a = f['data/a'].value
    #b = f['data/b'].value
    #c = f['data/c'].value
    #d = f['data/d'].value
    
with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter5.hdf5', 'r') as f:
    mz = f['data/magnitudes/m'].value
    
    idx1 = f['data/camtrans/idx'].value
    z = f['data/camtrans/z'].value
    
    idx2 = f['data/intrapix/idx'].value
    a = f['data/intrapix/a'].value
    b = f['data/intrapix/b'].value
    c = f['data/intrapix/c'].value
    d = f['data/intrapix/d'].value
    
pg = PolarGrid(13500, 720)
z = pg.put_values_on_grid(z, idx1, np.nan)
z = z - np.nanmean(z)

pg = PolarGrid(270, 720)
a = pg.put_values_on_grid(a, idx2, np.nan)
b = pg.put_values_on_grid(b, idx2, np.nan)
c = pg.put_values_on_grid(c, idx2, np.nan)
d = pg.put_values_on_grid(d, idx2, np.nan)

#ha = pg.bins1
#dec = pg.bins2

#from altaz2hadec import altaz2hadec
#from hadec2altaz import hadec2altaz
#from mpl_toolkits.basemap import Basemap, shiftgrid, cm

#az0 = 87.96
#alt0 = 49.07

#ha0, dec0 = altaz2hadec(alt0, az0, lat = 28.76)
#th0 = 270.84
#x0 = 2004.5
#y0 = 1299.8

#ha = ha-180
#ha0 = ha0-180

#m = Basemap(width=4008*2, height=2672*2, rsphere=24/9e-3, projection='gnom', lat_0=dec0, lon_0=ha0)
#topodat = m.transform_scalar(z[1:-1,1:-1].T,ha,dec,167*3,167*2)
#im = m.imshow(topodat)
#plt.show()

#az0 = 87.96
#alt0 = 49.07
#th0 = 270.84
#x0 = 2004.5
#y0 = 1299.8

#ha, dec = np.meshgrid(ha, dec)

#alt, az = hadec2altaz(ha, dec, lat = 28.76)
#print np.amin(az), np.amax(az)
#plt.contour(alt, [-80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80])
#plt.contour(az, [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345])
#plt.show()

#az0 = az0 - 180
#az = az - 180
#print z.shape, alt.shape, az.shape
#m = Basemap(width=4008, height=2672, rsphere=24/9e-3, projection='gnom', lat_0=alt0, lon_0=az0)
##topodat = m.transform_scalar(z[1:-1,1:-1].T,ha,dec,167*3,167*2)
#im = m.pcolormesh(az, alt, z[1:-1,1:-1].T, latlon=True)
#plt.show()


#exit()


#pg = PolarGrid(13500, 720)  
#zm = np.ma.masked_invalid(z)

#plt.figure(figsize=(13,12))

#plt.suptitle('Transmission', size='xx-large')

#gs = gridspec.GridSpec(2, 2, height_ratios = [.5,10], width_ratios=[10,.5])

#ax1 = plt.subplot(gs[1,0])
#im = plt.pcolormesh(pg.bins1/15., pg.bins2, zm[1:-1,1:-1].T, cmap=viridis, vmin=-.5, vmax=1.5)
#plt.xlim(270./15, 347./15)
#plt.ylim(-20, 65)

#plt.xlabel('Hour Angle')
#plt.ylabel('Declination')

#ax2 = plt.subplot(gs[1,1])
#plt.colorbar(im, cax=ax2).set_label('z')

#plt.tight_layout()
#plt.show()

# Magnitude calibration.
plt.figure(figsize=(13,12))
    
plt.suptitle('Magnitudes', size='xx-large')

gs = gridspec.GridSpec(3, 1, height_ratios = [.5,10,10])
    
plt.subplot(gs[1])
plt.plot(vmag, mz, '.', alpha=.2)
plt.xlabel('V magnitude')
plt.ylabel('m')

plt.subplot(gs[2])
plt.plot(dec, mz-vmag, '.', alpha=.2)
plt.xlabel('Declination')
plt.ylabel('m - V')

plt.tight_layout()
plt.show()

# Transmission profile.
#z = z.reshape((13502, 722))
z = delete_allnan(z)

plt.figure(figsize=(13,12))

plt.suptitle('Transmission', size='xx-large')

gs = gridspec.GridSpec(2, 2, height_ratios = [.5,10], width_ratios=[10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
im = plt.imshow(z.T, aspect='auto', cmap=viridis, vmin=-.5, vmax=1.5)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

ax2 = plt.subplot(gs[1,1])
plt.colorbar(im, cax=ax2).set_label('z')

plt.tight_layout()
plt.show()

# Intrapixel variations.
Ax = np.sqrt(a**2 + b**2)
Ay = np.sqrt(c**2 + d**2)
phix = np.arctan2(b, a)
phiy = np.arctan2(d, c)

#Ax = Ax.reshape((272, 722))
#Ay = Ay.reshape((272, 722))
#phix = phix.reshape((272, 722))
#phiy = phiy.reshape((272, 722))
Ax = delete_allnan(Ax)
Ay = delete_allnan(Ay)
phix = delete_allnan(phix)
phiy = delete_allnan(phiy)

plt.figure(figsize=(13,12))

plt.suptitle('Intrapixel Variations', size='xx-large')

gs = gridspec.GridSpec(3, 3, height_ratios = [.5,10,10], width_ratios=[10,10,.5])

ax1 = plt.subplot(gs[1,0], xticks=[], yticks=[])
plt.title('x position', size='x-large')
plt.imshow(Ax.T, aspect='auto', cmap=viridis, vmin=0, vmax=.1)
plt.ylabel('Declination')

plt.subplot(gs[1,1], sharex=ax1, sharey=ax1)
plt.title('y position', size='x-large')
im = plt.imshow(Ay.T, aspect='auto', cmap=viridis, vmin=0, vmax=.1)

ax2 = plt.subplot(gs[1,2])
plt.colorbar(im, cax=ax2).set_label('Amplitude')

plt.subplot(gs[2,0], sharex=ax1, sharey=ax1)
plt.imshow(phix.T, aspect='auto', cmap=viridis, vmin=-np.pi, vmax=np.pi)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

plt.subplot(gs[2,1], sharex=ax1, sharey=ax1)
im = plt.imshow(phiy.T, aspect='auto', cmap=viridis, vmin=-np.pi, vmax=np.pi)
plt.xlabel('Hour Angle')

ax3 = plt.subplot(gs[2,2])
plt.colorbar(im, cax=ax3).set_label('Phase')

plt.tight_layout()
plt.show()
