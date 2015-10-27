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

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

filelist = glob.glob('/data2/talens/3mEast/fLC_2015060?LPE.hdf5')
filelist = np.sort(filelist)
print filelist

ipx1 = np.array([])
ipx2 = np.array([])
mz = np.array([])
ma = np.array([])
mb = np.array([])
mc = np.array([])
md = np.array([])

for ind in range(401, 481):

    mascc = np.array([])
    mra = np.array([])
    mdec = np.array([])
    mnobs = np.array([])
    mlst = np.array([])
    mx = np.array([])
    my = np.array([])
    mflux = np.array([])
    meflux = np.array([])
    msky = np.array([])
    mflag = np.array([])

    for filename in filelist:
        f = fLCfile(filename)
        ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])

        pg = PolarGrid(13500, 720)
        decidx = pg.find_decidx(dec)
        
        here = (decidx == ind)
        ascc = ascc[here]
        ra = ra[here]
        dec = dec[here]
        nobs = nobs[here]

        lst, x, y, flux, eflux, sky, flag = f.read_data(['lst', 'x', 'y', 'flux0', 'eflux0', 'sky', 'flag'], ascc, nobs)

        mascc = np.append(mascc, ascc)
        mra = np.append(mra, ra)
        mdec = np.append(mdec, dec)
        mnobs = np.append(mnobs, nobs)
        mlst = np.append(mlst, lst)
        mx = np.append(mx, x)
        my = np.append(my, y)
        mflux = np.append(mflux, flux)
        meflux = np.append(meflux, eflux)
        msky = np.append(msky, sky)
        mflag = np.append(mflag, flag)
        
    ascc = mascc
    ra = mra
    dec = mdec
    nobs = mnobs
    lst = mlst
    x = mx
    y = my
    flux = mflux
    eflux = meflux
    sky = msky
    flag = mflag

    nobs = nobs.astype('int')

    # Build indices
    ascc = np.repeat(ascc, nobs)
    ascc, staridx = np.unique(ascc, return_inverse=True)

    ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
    pg = PolarGrid(13500, 720)
    camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))

    pg2 = PolarGrid(270, 720)
    intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec, nobs))

    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    x = x[here]
    y = y[here]
    flux = flux[here]
    eflux = eflux[here]

    staridx = staridx[here]
    camtransidx = camtransidx[here]
    intrapixidx = intrapixidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    # Get unique indices.
    staridx, staruni = np.unique(staridx, return_inverse=True)
    camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
    intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)

    # Calculate a model fit to the data.
    m, z, a, b, c, d, niter, chisq, npoints, npars = systematics_dev.trans_ipx(staruni, camtransuni, intrapixuni, mag, emag, x, y, verbose=True, use_weights=False)
    print niter, chisq, npoints, npars
        
    ipx1 = np.append(ipx1, camtransidx)
    ipx2 = np.append(ipx2, intrapixidx)
    mz = np.append(mz, z)
    ma = np.append(ma, a)
    mb = np.append(mb, b)
    mc = np.append(mc, c)
    md = np.append(md, d)
        
with h5py.File('/data2/talens/3mEast/LBtests/camip_2015060LPE.hdf5') as f:
    
    grp = f.create_group('data')
    grp.create_dataset('camtransidx', data=ipx1)
    grp.create_dataset('z', data=mz)
    
    grp.create_dataset('intrapixidx', data=ipx2)
    grp.create_dataset('a', data=ma)
    grp.create_dataset('b', data=mb)
    grp.create_dataset('c', data=mc)
    grp.create_dataset('d', data=md)
        
        
        
        
