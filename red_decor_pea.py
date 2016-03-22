#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

import matplotlib.pyplot as plt

from package import misc
from package import IO
from package.coordinates import grids
from package.plotting import viridis

from pea_grid import polar_eqarea_caps

def cdecor(idx1, idx2, idx3, idx4, mag, error, x, y, maxiter=100, dtol=1e-3, verbose=True):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars3 = np.amax(idx3) + 1
    npars4 = 4*(np.amax(idx4) + 1)
    npars = npars2 + npars3 + npars4
    
    # Create arrays.
    weights = 1./error**2
    par3 = np.zeros(npars3)
    
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
    
    ipx_x = 0
    ipx_y = 0
    b = np.zeros(npars4/4)
    d = np.zeros(npars4/4)
    
    err1 = np.zeros(npars1)
    err3 = np.zeros(npars3)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par2 = np.bincount(idx2, weights*(mag - par3[idx3] - ipx_x - ipx_y))/np.bincount(idx2, weights)
        par3 = np.bincount(idx3, weights*(mag - par2[idx2] - ipx_x - ipx_y))/np.bincount(idx3, weights)
    
        sol1 = par2[idx2] + par3[idx3]
    
        a = np.bincount(idx4, weights*(mag - sol1 - b[idx4]*csx - ipx_y)*snx)/np.bincount(idx4, weights*snx**2)
        b = np.bincount(idx4, weights*(mag - sol1 - a[idx4]*snx - ipx_y)*csx)/np.bincount(idx4, weights*csx**2)
        
        ipx_x = a[idx4]*snx + b[idx4]*csx
        
        c = np.bincount(idx4, weights*(mag - sol1 - ipx_x - d[idx4]*csy)*sny)/np.bincount(idx4, weights*sny**2)
        d = np.bincount(idx4, weights*(mag - sol1 - ipx_x - c[idx4]*sny)*csy)/np.bincount(idx4, weights*csy**2)
        
        ipx_y = c[idx4]*sny + d[idx4]*csy
    
        par4 = np.vstack([a, b, c, d]).T
        
        if (niter > 4):
        
            err1 = find_sigma(idx1, mag - sol1 - ipx_x - ipx_y, np.sqrt(error**2 + err3[idx3]**2))
            err3 = find_sigma(idx3, mag - sol1 - ipx_x - ipx_y, np.sqrt(error**2 + err1[idx1]**2))
            weights = 1/(error**2 + (err1**2)[idx1] + (err3**2)[idx3])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit = np.nanmax(np.abs(sol1 - sol1_old))
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            dcrit3 = np.nanmax(np.abs(par3 - par3_old))
            dcrit4 = np.nanmax(np.abs(par4 - par4_old))
            
            print dcrit
            
            if (dcrit < dtol) & (dcrit4 < dtol):
                print 'Solution has converged, ending the iterations.'
                break
        
        # Check if the solution is oscillating?
        if (niter > 1):
            
            dcrit = np.nanmax(np.abs(sol1 - sol1_older))
            dcrit2 = np.nanmax(np.abs(par2 - par2_older))
            dcrit3 = np.nanmax(np.abs(par3 - par3_older))
            dcrit4 = np.nanmax(np.abs(par4 - par4_older))
            
            if (dcrit < dtol) & (dcrit4 < dtol):
                print 'Solution is oscillating, ending the iterations.'
                break
        
        if (niter > 0):
            sol1_older = np.copy(sol1_old)
            par2_older = np.copy(par2_old)
            par3_older = np.copy(par3_old)
            par4_older = np.copy(par4_old)
        
        sol1_old = np.copy(sol1)
        par2_old = np.copy(par2)
        par3_old = np.copy(par3)
        par4_old = np.copy(par4)
    
    # Compute the chi-square of the fit.
    chisq = weights*(mag - sol1 - ipx_x - ipx_y)**2        
    chisq = np.sum(chisq)
    
    return par2, par3, par4, err1, err3, Quality(niter, chisq, npoints, npars)
    

def sigma_function(idx, ressq, errsq, err):
    
    weights = 1/(errsq + (err**2)[idx])
    term = ressq*weights**2 - weights
    term = np.bincount(idx, term)
    
    return term

def find_sigma(idx, residuals, error, maxiter=10):
    
    # Search for a solution between 0 and 2.
    N = np.amax(idx) + 1
    err1 = np.zeros(N)
    err2 = np.full(N, 2)
    
    ressq = residuals*residuals
    errsq = error*error
    
    # Compute the value of the function at the beginning the interval.
    diff1 = sigma_function(idx, ressq, errsq, err1)
    args1, = np.where(diff1 < 1e-10)

    # Compute the value of the function at the end the interval.
    diff2 = sigma_function(idx, ressq, errsq, err2)
    args2, = np.where(diff2 > 1e-10)

    # Find the solution.
    for niter in range(maxiter):
        
        err3 = (err2 + err1)/2.
        diff3 = sigma_function(idx, ressq, errsq, err3)
        
        err1 = np.where(diff3 > 1e-10, err3, err1)
        err2 = np.where(diff3 > 1e-10, err2, err3)
    
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.

    return err3

def main():
    
    filename = '/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5'
    
    f = IO.fLCfile(filename)
    ascc, ra, dec, vmag, nobs = f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
    nobs = nobs.astype('int')
    staridx = np.arange(len(ascc))
    
    ring, cell, N = polar_eqarea_caps(ra, dec, 23)
    
    print np.unique(ring)
    print N[15]
    
    select = (ring == 15)
    ascc = ascc[select]
    ra = ra[select]
    dec = dec[select]
    vmag = vmag[select]
    nobs = nobs[select]
    cell = cell[select]
    staridx = staridx[select]

    lstseq, lst, flux, eflux, sky, flag, x, y = f.read_data(['lstseq', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], ascc, nobs)
    lstseq = lstseq.astype('int')
    
    ra = np.repeat(ra, nobs)
    dec = np.repeat(dec, nobs)
    cell = np.repeat(cell, nobs)  
    vmag = np.repeat(vmag, nobs)
    staridx = np.repeat(staridx, nobs)
    
    ha = np.mod(lst*15. - ra, 360.)
    camgrid = grids.PolarGrid(13500, 720)
    haidx1, decidx = camgrid.radec2idx(ha, dec)
    ipxgrid = grids.PolarGrid(270, 720)
    haidx2, decidx = ipxgrid.radec2idx(ha, dec)
    
    select = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    lstseq = lstseq[select]
    lst = lst[select]
    flux = flux[select]
    eflux = eflux[select]
    cell = cell[select]
    decidx = decidx[select]
    haidx1 = haidx1[select]
    haidx2 = haidx2[select]
    vmag = vmag[select]
    x = x[select]
    y = y[select]
    staridx = staridx[select]
    
    mag, emag = misc.flux2mag(flux, eflux)
    mag = mag - vmag
    lstseq = lstseq - np.amin(lstseq)
    lstlen = np.ptp(lstseq) + 1
    
    camidx = np.ravel_multi_index((haidx1, decidx), (13502, 722))
    skyidx = np.ravel_multi_index((cell, lstseq), (N[15], lstlen)) 
    ipxidx = np.ravel_multi_index((haidx2, decidx), (272, 722))

    staridx, idx1 = np.unique(staridx, return_inverse=True)
    camidx, idx2 = np.unique(camidx, return_inverse=True)
    skyidx, idx3 = np.unique(skyidx, return_inverse=True)
    ipxidx, idx4 = np.unique(ipxidx, return_inverse=True)

    cam = np.full((13502*722,), fill_value=np.nan) 
    sky = np.full((N[15]*lstlen,), fill_value=np.nan)
    err3 = np.full((N[15]*lstlen,), fill_value=np.nan)
    A = np.full((272*722, 4), fill_value=np.nan)

    cam[camidx], sky[skyidx], A[ipxidx], err1, err3[skyidx], quality = cdecor(idx1, idx2, idx3, idx4, mag, emag, x, y, maxiter=25)
    
    cam = cam.reshape((13502, 722))
    sky = sky.reshape((N[15], lstlen))
    err3 = err3.reshape((N[15], lstlen))
    A = A.reshape((272, 722, 4))

    plt.imshow(cam.T, aspect='auto', interpolation='nearest', cmap=viridis)
    plt.colorbar()
    plt.show()
    
    ax = plt.subplot(221)
    plt.imshow(A[:,:,0].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    plt.colorbar()
    
    plt.subplot(222, sharex=ax, sharey=ax)
    plt.imshow(A[:,:,1].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    plt.colorbar()
    
    plt.subplot(223, sharex=ax, sharey=ax)
    plt.imshow(A[:,:,2].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    plt.colorbar()
    
    plt.subplot(224, sharex=ax, sharey=ax)
    plt.imshow(A[:,:,3].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    plt.colorbar()
    
    plt.show()
    
    ax = plt.subplot(211)
    plt.imshow(sky, aspect='auto', interpolation='nearest', vmin=-.5, vmax=.5, cmap=viridis)
    plt.colorbar()

    plt.subplot(212, sharex=ax, sharey=ax)
    plt.imshow(err3, aspect='auto', interpolation='nearest', vmin=0, vmax=.2, cmap=viridis)
    plt.colorbar()
    plt.show()

    return

if __name__ == '__main__':
    main()
