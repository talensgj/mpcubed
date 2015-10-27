#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

filelist = glob.glob('/data2/talens/3mEast/fLC_201507??LPE.hdf5')
filelist = np.sort(filelist)
print filelist

with h5py.File('/data2/talens/3mEast/LBtests/camip_201506LPE.hdf5') as f:
    
    camtransidx = f['data/camtransidx'].value
    z = f['data/z'].value
    intrapixidx = f['data/intrapixidx'].value
    a = f['data/a'].value
    b = f['data/b'].value
    c = f['data/c'].value
    d = f['data/d'].value

pg = PolarGrid(13500, 720)
z = pg.put_values_on_grid(z, camtransidx, np.nan)
pg2 = PolarGrid(270, 720)
a = pg.put_values_on_grid(a, intrapixidx, np.nan)
b = pg.put_values_on_grid(b, intrapixidx, np.nan)
c = pg.put_values_on_grid(c, intrapixidx, np.nan)
d = pg.put_values_on_grid(d, intrapixidx, np.nan)

z = np.ravel(z)
a = np.ravel(a)
b = np.ravel(b)
c = np.ravel(c)
d = np.ravel(d)

mascc = np.array([])
mra = np.array([])
mdec = np.array([])
mnobs = np.array([])
mjdmid = np.array([])
mlstidx = np.array([])
mflux = np.array([])
meflux = np.array([])
msky = np.array([])
mflag = np.array([])
mx = np.array([])
my = np.array([])
mlst = np.array([])

for filename in filelist:
    f = fLCfile(filename)
    ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])

    hg = HealpixGrid(8)
    skyidx = hg.find_gridpoint(ra, dec)
    
    here = (skyidx == 266)
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    nobs = nobs[here]

    jdmid, lstidx, flux, eflux, sky, flag, x, y, lst = f.read_data(['jdmid', 'lstidx', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y', 'lst'], ascc, nobs)

    mascc = np.append(mascc, ascc)
    mra = np.append(mra, ra)
    mdec = np.append(mdec, dec)
    mnobs = np.append(mnobs, nobs)
    
    mjdmid = np.append(mjdmid, jdmid)
    mlstidx = np.append(mlstidx, lstidx)
    mflux = np.append(mflux, flux)
    meflux = np.append(meflux, eflux)
    msky = np.append(msky, sky)
    mflag = np.append(mflag, flag)
    mx = np.append(mx, x)
    my = np.append(my, y)
    mlst = np.append(mlst, lst)
    
ascc = mascc
ra = mra
dec = mdec
nobs = mnobs
    
jdmid = mjdmid
lstidx = mlstidx
flux = mflux
eflux = meflux
sky = msky
flag = mflag
x = mx
y = my
    
nobs = nobs.astype('int')
lstidx = lstidx.astype('int')

# Build indices
ascc = np.repeat(ascc, nobs)
ascc, staridx = np.unique(ascc, return_inverse=True)

dayidx = np.floor(jdmid)
dayidx = dayidx - np.amin(dayidx)
dayidx = dayidx.astype('int')

skytransidx = np.ravel_multi_index((dayidx, lstidx), (31, 13500))

ha = np.mod(lst*15 - np.repeat(ra, nobs), 360)
idx1 = pg.find_gridpoint(ha, np.repeat(dec, nobs))
idx2 = pg2.find_gridpoint(ha, np.repeat(dec, nobs))

# Flag bad data.
here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
flux = flux[here]
eflux = eflux[here]

staridx = staridx[here]
skytransidx = skytransidx[here]
idx1 = idx1[here]
idx2 = idx2[here]
x = x[here]
y = y[here]

# Convert flux to magnitudes
mag = 25 - 2.5*np.log10(flux)
emag = 2.5/np.log(10)*eflux/flux

# Get unique indices.
staridx, staruni = np.unique(staridx, return_inverse=True)
skytransidx, skytransuni = np.unique(skytransidx, return_inverse=True)

# Calculate a model fit to the data.
m, z, niter, chisq, npoints, npars = systematics_dev.trans(staruni, skytransuni, mag, emag, verbose=True, use_weights=False)

plt.plot(z, '.')
plt.show()

for i in range(len(staridx)):
    here = (staruni == i)
    
    ax = plt.subplot(211)
    plt.title(ascc[staridx[i]])
    plt.plot(mag[here], '.')
    plt.plot(z[skytransuni[here]]+m[i], '.')
    plt.subplot(212, sharex=ax, sharey=ax)
    plt.plot(mag[here] - z[skytransuni[here]], '.')
    plt.show()
    plt.close()
    
    
    
