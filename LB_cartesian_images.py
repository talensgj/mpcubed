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

from altaz2hadec import altaz2hadec
from hadec2altaz import hadec2altaz
from altaz2xy import altaz2xy
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
    
az0 = 87.96
alt0 = 49.07
th0 = 270.84
x0 = 2004.5
y0 = 1299.8
    
with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter4.hdf5', 'r') as f:
    #mz = f['data/magnitudes/m'].value
    
    idx1 = f['data/camtrans/idx'].value
    z = f['data/camtrans/z'].value
    
    idx2 = f['data/intrapix/idx'].value
    a = f['data/intrapix/a'].value
    b = f['data/intrapix/b'].value
    c = f['data/intrapix/c'].value
    d = f['data/intrapix/d'].value
    
pg = PolarGrid(13500, 720)
z = pg.put_values_on_grid(z, idx1, np.nan)
a = z - np.nanmean(z)

#pg = PolarGrid(270, 720)
#a = pg.put_values_on_grid(a, idx2, np.nan)
#b = pg.put_values_on_grid(b, idx2, np.nan)
#c = pg.put_values_on_grid(c, idx2, np.nan)
#d = pg.put_values_on_grid(d, idx2, np.nan)

#a = np.arctan2(b, a)

a = a[1:-1,1:-1].T
ha = pg.bins1
dec = pg.bins2
ha, dec = np.meshgrid(ha, dec)
    
ind1, = np.where(np.all(np.isnan(a), axis=0))
ind2, = np.where(np.all(np.isnan(a), axis=1))
    
a = np.delete(a, ind1, axis=1)
a = np.delete(a, ind2, axis=0)
ha = np.delete(ha, ind1, axis=1)
ha = np.delete(ha, ind2, axis=0)
dec = np.delete(dec, ind1, axis=1)
dec = np.delete(dec, ind2, axis=0)

alt, az = hadec2altaz(ha, dec, lat=28.76)
#print np.mean(alt[~np.isnan(a)]), np.mean(az[~np.isnan(a)])

#m = Basemap(width=4008, height=2672, rsphere=24/9e-3, projection='gnom', lat_0=alt0, lon_0=az0)
#x, y = m(az, alt)
#x = x - 2004
#y = y - 1336

#xi = x*np.cos(th0) + y*np.sin(th0)
#yi = -x*np.sin(th0) + y*np.cos(th0)

#x = -xi + x0
#y = yi + y0

a = np.ma.masked_invalid(a)
#m.pcolormesh(x, y, a, vmin=-.5, vmax=1.5, cmap=viridis)
#m.contour(x, y, ha, np.linspace(0, 345, 360/15))
#m.contour(x, y, dec, np.linspace(-80, 80, 17))
#plt.show()
tmp = alt.shape
x, y, z = altaz2xy(np.ravel(alt), np.ravel(az), alt0, az0, th0, x0, y0, 24/9e-3)
x = x.reshape(tmp)
y = y.reshape(tmp)

plt.pcolormesh(x, y, a)
plt.show()

#xi = x*np.cos(th0) + y*np.sin(th0)
#yi = -x*np.sin(yh0) + y*np.cos(th0)

plt.imshow(a, aspect='auto')
#plt.contour(ha, np.linspace(0, 345, 360/15), extent=(0,360,-90,90), aspect='auto')
#plt.contour(dec, np.linspace(-80, 80, 17), extent=(0,360,-90,90), aspect='auto')

#plt.contour(az, np.linspace(0, 345, 360/15), extent=(0,360,-90,90), aspect='auto')
#plt.contour(alt, np.linspace(-80, 80, 17), extent=(0,360,-90,90), aspect='auto')

plt.contour(x, np.linspace(0,4008,4008/167+1), aspect='auto')
plt.contour(y, np.linspace(0,2672,2672/167+1), aspect='auto')

plt.show()
exit()


az = np.linspace(-180, 180, 100)
alt = np.linspace(-90, 90, 100)
m = Basemap(width=4008, height=2672, rsphere=24/9e-3, projection='gnom', lat_0=alt0, lon_0=az0)
x, y = m(*np.meshgrid(az, alt))

plt.contour(x, np.linspace(0,4008,4008/167+1), colors='k')
plt.contour(y, np.linspace(0,2672,2672/167+1), colors='r')
plt.show()


#m = Basemap(width=4008*2, height=2672*2, rsphere=24/9e-3, projection='gnom', lat_0=dec0, lon_0=ha0)
#topodat = m.transform_scalar(z[1:-1,1:-1].T,ha,dec,167*3,167*2)
#im = m.imshow(topodat)
#plt.show()

#m = Basemap(width=4008, height=2672, rsphere=24/9e-3, projection='gnom', lat_0=alt0, lon_0=az0)
##topodat = m.transform_scalar(z[1:-1,1:-1].T,ha,dec,167*3,167*2)
#im = m.pcolormesh(az, alt, z[1:-1,1:-1].T, latlon=True)
#plt.show()
