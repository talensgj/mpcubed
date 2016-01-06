#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid, CartesianGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile
from usefull_functions_dev import mad

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis
import matplotlib.gridspec as gridspec

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/3mEast/variations/cam_20150611LPE+flags+ipx+weights.hdf5', 'r') as f:
    
    hdr = f['header']
    niter = hdr.attrs['niter']
    chisq = hdr.attrs['chisq']
    npoints = hdr.attrs['npoints']
    npars = hdr.attrs['npars']
    
    grp = f['data']
    staridx = grp['staridx'].value
    m = grp['m'].value
    
    camtransidx = grp['camtransidx'].value
    z = grp['z'].value

    intrapixidx = grp['intrapixidx'].value
    a = grp['a'].value
    b = grp['b'].value
    c = grp['c'].value
    d = grp['d'].value

pg = PolarGrid(13500, 720)
z = pg.put_values_on_grid(z, camtransidx, np.nan)
z = np.ravel(z)

pg2 = PolarGrid(270, 720)
a = pg2.put_values_on_grid(a, intrapixidx, np.nan)
b = pg2.put_values_on_grid(b, intrapixidx, np.nan)
c = pg2.put_values_on_grid(c, intrapixidx, np.nan)
d = pg2.put_values_on_grid(d, intrapixidx, np.nan)
a = np.ravel(a)
b = np.ravel(b)
c = np.ravel(c)
d = np.ravel(d)

#cg = CartesianGrid(250, 167)
#a = cg.put_values_on_grid(a, intrapixidx, np.nan)
#b = cg.put_values_on_grid(b, intrapixidx, np.nan)
#c = cg.put_values_on_grid(c, intrapixidx, np.nan)
#d = cg.put_values_on_grid(d, intrapixidx, np.nan)
#a = np.ravel(a)
#b = np.ravel(b)
#c = np.ravel(c)
#d = np.ravel(d)

# Read data.
f = fLCfile('/data2/talens/3mEast/fLC_20150611LPE.hdf5')
ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
nobs = nobs.astype('int')

array = np.full(len(ascc), fill_value=np.nan)
array[staridx] = m
m = array

stat = np.zeros(len(ascc))
for i in range(len(ascc)):
    lst, flux, x, y = f.read_data(['lst', 'flux0', 'x', 'y'], [ascc[i]], [nobs[i]])

    ha = np.mod(lst*15 - ra[i], 360.)
    idx1 = pg.find_gridpoint(ha, np.repeat(dec[i], nobs[i]))
    idx2 = pg2.find_gridpoint(ha, np.repeat(dec[i], nobs[i]))
    #idx2 = cg.find_gridpoint(x, y)
    
    mag = 25 - 2.5*np.log10(flux)
    cmag = mag - z[idx1] - a[idx2]*np.sin(2*np.pi*x) - b[idx2]*np.cos(2*np.pi*x) - c[idx2]*np.sin(2*np.pi*y) - d[idx2]*np.sin(2*np.pi*y)
    
    #plt.subplot(211)
    #plt.title(ascc[i] + ', %.2f'%m[i])
    #plt.plot(mag, '.')
    #plt.plot(m[i] + z[idx1] + a[idx2]*np.sin(2*np.pi*x) + b[idx2]*np.cos(2*np.pi*x) + c[idx2]*np.sin(2*np.pi*y) + d[idx2]*np.sin(2*np.pi*y), '.')
    #plt.subplot(212)
    #plt.plot(cmag, '.')
    #plt.show()
    
    print np.nanmean(cmag), np.nanstd(cmag), mad(cmag)
    
    stat[i] = mad(cmag)
    
plt.semilogy(vmag, stat, '.', alpha=.2)
plt.ylim(1e-3, 1)
plt.xlim(8.4, 2)
plt.xlabel('V')
plt.ylabel('MAD')
plt.show()
    
    
    
