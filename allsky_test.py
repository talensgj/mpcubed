#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid, CartesianGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

# Read data.
f = fLCfile('/data2/talens/3mEast/fLC_20150611LPE.hdf5')
ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
lst, x, y, flux, eflux, sky, flag = f.read_data(['lst', 'x', 'y', 'flux0', 'eflux0', 'sky', 'flag'], ascc, nobs)

nobs = nobs.astype('int')

# Build indices.
staridx = np.repeat(np.arange(len(ascc)), nobs)

ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
pg = PolarGrid(13500, 720)
camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))

pg2 = PolarGrid(270, 720)
intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec, nobs))

#cg = CartesianGrid(250, 167, margin=0)
#intrapixidx = cg.find_gridpoint(x, y)

# Flag bad data.
here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
x = x[here]
y = y[here]
flux = flux[here]
eflux = eflux[here]

staridx = staridx[here]
camtransidx = camtransidx[here]
intrapixidx = intrapixidx[here]

# Convert flux to magnitudes.
mag = 25 - 2.5*np.log10(flux)
emag = 2.5/np.log(10)*np.abs(eflux/flux)

# Get unique indices.
staridx, staruni = np.unique(staridx, return_inverse=True)
camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)

# Calculate a model fit to the data.
#m, z, niter, chisq, npoints, npars = systematics_dev.trans(staruni, camtransuni, mag, emag, verbose=True, use_weights=False)
#m, z, sigma, niter, chisq, npoints, npars = systematics_dev.trans(staruni, camtransuni, mag, emag, verbose=True, use_weights=True)
#m, z, a, b, c, d, niter, chisq, npoints, npars = systematics_dev.trans_ipx(staruni, camtransuni, intrapixuni, mag, emag, x, y, verbose=True, use_weights=False)
m, z, sigma, a, b, c, d, niter, chisq, npoints, npars = systematics_dev.trans_ipx(staruni, camtransuni, intrapixuni, mag, emag, x, y, verbose=True, use_weights=True)

with h5py.File('/data2/talens/3mEast/variations/cam_20150611LPE+flags+ipx+weights.hdf5') as f:
    
    hdr = f.create_group('header')
    hdr.attrs['camgrid'] = 'polar'
    hdr.attrs['nx_cam'] = 13500
    hdr.attrs['ny_cam'] = 720
    
    hdr.attrs['ipxgrid'] = 'cartesian'
    hdr.attrs['nx_ipx'] = 250
    hdr.attrs['ny_ipx'] = 167
    
    hdr.attrs['niter'] = niter
    hdr.attrs['chisq'] = chisq
    hdr.attrs['npoints'] = npoints
    hdr.attrs['npars'] = npars
    
    grp = f.create_group('data')
    grp.create_dataset('staridx', data=staridx)
    grp.create_dataset('m', data=m)
    grp.create_dataset('sigma', data=sigma)
    
    grp.create_dataset('camtransidx', data=camtransidx)
    grp.create_dataset('z', data=z)
    
    grp.create_dataset('intrapixidx', data=intrapixidx)
    grp.create_dataset('a', data=a)
    grp.create_dataset('b', data=b)
    grp.create_dataset('c', data=c)
    grp.create_dataset('d', data=d)


 
