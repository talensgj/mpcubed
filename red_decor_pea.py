#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
from scipy import linalg

from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

import matplotlib.pyplot as plt

from package import misc
from package import IO
from package.coordinates import grids
from package.plotting import viridis

from pea_grid import polar_eqarea_caps

def _sigma_function_vmag(idx1, idx2, value, error, par1, sigma1, err):
    
    weights = 1/(error**2 + (sigma1**2)[idx1] + (err**2)[idx2])
    diff = np.bincount(idx2, weights**2*(value - par1[idx1])**2 - weights)
    
    return diff
    
def _find_sigma_vmag(idx1, idx2, value, error, par1, sigma1, maxiter = 10):
    
    # Search for a solution between 0 and 2.
    N = np.amax(idx2) + 1
    err1 = np.zeros(N)
    err2 = 2*np.ones(N)
    
    # Compute the value of the function at the beginning the interval.
    diff1 = _sigma_function_vmag(idx1, idx2, value, error, par1, sigma1, err1)
    args1, = np.where(diff1 < 1e-10)
    
    # Compute the value of the function at the end the interval.
    diff2 = _sigma_function_vmag(idx1, idx2, value, error, par1, sigma1, err2)
    args2, = np.where(diff2 > 1e-10)
    
    # Find the solution.
    for niter in range(maxiter):
        
        err3 = (err1 + err2)/2.
        diff3 = _sigma_function_vmag(idx1, idx2, value, error, par1, sigma1, err3)
    
        err1 = np.where(diff3 > 1e-10, err3, err1)
        err2 = np.where(diff3 > 1e-10, err2, err3)
        
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.
    
    return err3
    
def _sigma_function(idx1, idx2, value, error, par1, sigma1, err):
    
    weights = 1/(error**2 + (sigma1**2)[idx1] + (err**2)[idx2])
    par2 = np.bincount(idx2, weights*(value - par1[idx1]))/np.bincount(idx2, weights)
    diff = np.bincount(idx2, weights**2*(value - par1[idx1] - par2[idx2])**2 - weights)
    
    return par2, diff
    
def _find_sigma(idx1, idx2, value, error, par1, sigma1, maxiter = 10):
    
    # Search for a solution between 0 and 2.
    N = np.amax(idx2) + 1
    err1 = np.zeros(N)
    err2 = 2*np.ones(N)
    
    # Compute the value of the function at the beginning the interval.
    par2, diff1 = _sigma_function(idx1, idx2, value, error, par1, sigma1, err1)
    args1, = np.where(diff1 < 1e-10)
    
    # Compute the value of the function at the end the interval.
    par2, diff2 = _sigma_function(idx1, idx2, value, error, par1, sigma1, err2)
    args2, = np.where(diff2 > 1e-10)
    
    # Find the solution.
    for niter in range(maxiter):
        
        err3 = (err1 + err2)/2.
        par2, diff3 = _sigma_function(idx1, idx2, value, error, par1, sigma1, err3)
    
        err1 = np.where(diff3 > 1e-10, err3, err1)
        err2 = np.where(diff3 > 1e-10, err2, err3)
        
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.
    
    par2, _ = _sigma_function(idx1, idx2, value, error, par1, sigma1, err3)
    
    return par2, err3

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
    
def fast_sigma_function(ressq, errsq, err):
    
    weights = 1/(errsq + err**2)
    term = ressq*weights**2 - weights
    term = np.sum(term)
    
    return term

def fast_find_sigma(idx, residuals, error, sort, maxiter=10):
    
    idx = idx[sort]
    residuals = residuals[sort]
    error = error[sort]
    
    strides = np.cumsum(np.bincount(idx))
    strides = np.append(0, strides)
    
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

    args, = np.where(diff1*diff2 < 0)

    err3 = np.zeros(N)
    for i in args:

        # Find the solution.
        for niter in range(maxiter):
            
            err3[i] = (err2[i] + err1[i])/2.
            diff3 = fast_sigma_function(ressq[strides[i]:strides[i+1]], errsq[strides[i]:strides[i+1]], err3[i])
            
            if (diff3 > 1e-10):
                err1[i] = err3[i]
            else:
                err2[i] = err3[i]
    
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.

    return err3

def cdecor(idx1, idx2, idx3, idx4, mag, error, x, y, maxiter=100, dtol=1e-3, verbose=True):
    
    sort = np.argsort(idx4)
    idx1 = idx1[sort]
    idx2 = idx2[sort]
    idx3 = idx3[sort]
    idx4 = idx4[sort]
    mag = mag[sort]
    error = error[sort]
    x = x[sort]
    y = y[sort]
    
    strides = np.cumsum(np.bincount(idx4))
    strides = np.append(0, strides)
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars3 = np.amax(idx3) + 1
    npars4 = 4*(np.amax(idx4) + 1)
    npars = npars2 + npars3 + npars4
    
    # Create arrays.
    weights = 1./error**2
    
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
    mat = np.vstack([snx, csx, sny, csy]).T

    ipx = np.zeros(len(mag))

    par1 = np.zeros(npars1)
    err1 = np.zeros(npars1)
    par3 = np.zeros(npars3)
    err3 = np.zeros(npars3)
    par4 = np.zeros((npars4/4, 4))
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par2 = np.bincount(idx2, weights*(mag - par3[idx3] - ipx))/np.bincount(idx2, weights)
        par3 = np.bincount(idx3, weights*(mag - par2[idx2] - ipx))/np.bincount(idx3, weights)
    
        err1 = _find_sigma_vmag(idx3, idx1, mag - par2[idx2] - ipx, error, par3, err3)
        par3, err3 = _find_sigma(idx1, idx3, mag - par2[idx2] - ipx, error, par1, err1)
        weights = 1/(error**2 + (err1**2)[idx1] + (err3**2)[idx3])
    
        sol1 = par2[idx2] + par3[idx3]
        
        res = mag - sol1
        for i in range(npars4/4):
            par4[i] = linalg.lstsq(mat[strides[i]:strides[i+1],:], res[strides[i]:strides[i+1]])[0]
            ipx[strides[i]:strides[i+1]] = np.dot(mat[strides[i]:strides[i+1],:], par4[i])
        
        #err1 = find_sigma(idx1, mag - sol1 - ipx, np.sqrt(error**2 + err3[idx3]**2))
        #err3 = find_sigma(idx3, mag - sol1 - ipx, np.sqrt(error**2 + err1[idx1]**2))
        #weights = 1/(error**2 + (err1**2)[idx1] + (err3**2)[idx3])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            tmp2 = np.abs(par2 - par2_old)
            tmp3 = np.abs(par3 - par3_old)
            tmp4 = np.abs(par4 - par4_old) 
            
            dcrit2 = np.nanmax(tmp2)
            dcrit3 = np.nanmax(tmp3)
            dcrit4 = np.nanmax(tmp4)
            
            if verbose:
                print '{}/{}, {:.3f}'.format(np.sum(tmp2<dtol), npars2, dcrit2)
                print '{}/{}, {:.3f}'.format(np.sum(tmp3<dtol), npars3, dcrit3)
                print '{}/{}, {:.3f}'.format(np.sum(tmp4<dtol), npars4, dcrit4)
            
            if (dcrit2 < dtol) & (dcrit3 < dtol) & (dcrit4 < dtol):
                print 'Solution has converged, ending the iterations.'
                break
        
        # Check if the solution is oscillating?
        if (niter > 1):
            
            tmp2 = np.abs(par2 - par2_older)
            tmp3 = np.abs(par3 - par3_older)
            tmp4 = np.abs(par4 - par4_older) 
            
            dcrit2 = np.nanmax(tmp2)
            dcrit3 = np.nanmax(tmp3)
            dcrit4 = np.nanmax(tmp4)
            
            if (dcrit2 < dtol) & (dcrit3 < dtol) & (dcrit4 < dtol):
                print 'Solution is oscillating, ending the iterations.'
                break
        
        if (niter > 0):

            par2_older = np.copy(par2_old)
            par3_older = np.copy(par3_old)
            par4_older = np.copy(par4_old)
        
        par2_old = np.copy(par2)
        par3_old = np.copy(par3)
        par4_old = np.copy(par4)
    
    # Compute the chi-square of the fit.
    chisq = weights*(mag - sol1 - ipx)**2        
    chisq = np.sum(chisq)
    
    return par2, par3, par4, err1, err3, Quality(niter, chisq, npoints, npars)
    
def main():
    
    filename = '/data2/talens/2015Q2/LPN/fLC_201506ALPN.hdf5'
    
    with h5py.File(filename, 'r') as f:
        
        lstmin = f['global'].attrs['lstmin'].astype('int')
        lstmax = f['global'].attrs['lstmax'].astype('int')
        
    f = IO.fLCfile(filename)
    ascc, ra, dec, vmag, nobs = f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
    nobs = nobs.astype('int')
    staridx = np.arange(len(ascc))
    
    lstlen = lstmax - lstmin + 1
    
    camgrid = grids.PolarGrid(13500, 720)
    ipxgrid = grids.PolarGrid(270, 720)
    ring, cell, N = polar_eqarea_caps(ra, dec, 23)
    
    skypatch = (np.cumsum(N) - N)[ring] + cell 
    Ntot = np.sum(N)
    
    err1 = np.full((len(ascc),), fill_value=np.nan)
    cam = np.full((13502*722,), fill_value=np.nan) 
    clouds = np.full((Ntot*lstlen,), fill_value=np.nan)
    err3 = np.full((Ntot*lstlen,), fill_value=np.nan)
    A = np.full((272*722, 4), fill_value=np.nan)
    
    nobs_cam = np.full((13502*722,), fill_value=np.nan) 
    nobs_A = np.full((272*722,), fill_value=np.nan) 
    nobs_clouds = np.full((Ntot*lstlen,), fill_value=np.nan)
    
    for i in range(25):

        # Select stars in this ring.
        select = (ring == i)
        ascc_ = ascc[select]
        
        if (len(ascc_) == 0):
            continue
        
        ra_ = ra[select]
        dec_ = dec[select]
        vmag_ = vmag[select]
        nobs_ = nobs[select]
        staridx_ = staridx[select]
        skypatch_ = skypatch[select]

        # Raed the data.
        lstseq, lst, flux, eflux, sky, flag, x, y = f.read_data(['lstseq', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], ascc_, nobs_)
        lstseq = lstseq.astype('int') - lstmin
        
        ra_ = np.repeat(ra_, nobs_)
        dec_ = np.repeat(dec_, nobs_)
        vmag_ = np.repeat(vmag_, nobs_)
        staridx_ = np.repeat(staridx_, nobs_)
        skypatch_ = np.repeat(skypatch_, nobs_)  
        
        ha = np.mod(lst*15. - ra_, 360.)
        haidx1, decidx = camgrid.radec2idx(ha, dec_)
        haidx2, decidx = ipxgrid.radec2idx(ha, dec_)
        
        # Remove bad data points.
        select = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
        lstseq = lstseq[select]
        lst = lst[select]
        flux = flux[select]
        eflux = eflux[select]
        x = x[select]
        y = y[select]
        vmag_ = vmag_[select]
        staridx_ = staridx_[select]
        skypatch_ = skypatch_[select]
        haidx1 = haidx1[select]
        haidx2 = haidx2[select]
        decidx = decidx[select]

        # Convert flux to magnitudes.
        mag, emag = misc.flux2mag(flux, eflux)
        mag = mag - vmag_
        
        camidx = np.ravel_multi_index((haidx1, decidx), (13502, 722))
        skyidx = np.ravel_multi_index((skypatch_, lstseq), (Ntot, lstlen)) 
        ipxidx = np.ravel_multi_index((haidx2, decidx), (272, 722))

        staridx_, idx1 = np.unique(staridx_, return_inverse=True)
        camidx, idx2 = np.unique(camidx, return_inverse=True)
        skyidx, idx3 = np.unique(skyidx, return_inverse=True)
        ipxidx, idx4 = np.unique(ipxidx, return_inverse=True)

        nobs_cam[camidx] = np.bincount(idx2)
        nobs_A[ipxidx] = np.bincount(idx4)
        nobs_clouds[skyidx] = np.bincount(idx3)

        print 'Computing solution for ring', i

        print len(mag)
        cam[camidx], clouds[skyidx], A[ipxidx], err1[staridx_], err3[skyidx], quality = cdecor(idx1, idx2, idx3, idx4, mag, emag, x, y, maxiter=50, verbose=False)
    
        print quality.niter, quality.chisq/(quality.npoints - quality.npars)
    
    nobs_cam = nobs_cam.reshape((13502, 722))
    cam = cam.reshape((13502, 722))
    nobs_clouds = nobs_clouds.reshape((Ntot, lstlen))
    clouds = clouds.reshape((Ntot, lstlen))
    err3 = err3.reshape((Ntot, lstlen))
    nobs_A = nobs_A.reshape((272, 722))
    A = A.reshape((272, 722, 4))
    
    with h5py.File('/data2/talens/pea4.hdf5') as f:

        # Write the header.
        hdr = f.create_group('header')
        
        #hdr.attrs['station'] = station
        #hdr.attrs['camera'] = camera
        
        #hdr.attrs['alt0'] = np.mean(alt0)
        #hdr.attrs['az0'] = np.mean(az0)
        #hdr.attrs['th0'] = np.mean(th0)
        #hdr.attrs['x0'] = np.mean(x0)
        #hdr.attrs['y0'] = np.mean(y0)
        
        #hdr.attrs['outer_maxiter'] = self.outer_maxiter
        #hdr.attrs['inner_maxiter'] = self.inner_maxiter
        #hdr.attrs['sigmas'] = self.sigmas
        #hdr.attrs['dtol'] = self.dtol
       
        #hdr.create_dataset('spatial/niter', data = self.spatial_niter, dtype = 'uint32')
        #hdr.create_dataset('spatial/chisq', data = self.spatial_chisq, dtype = 'float64')
        #hdr.create_dataset('spatial/npoints', data = self.spatial_npoints, dtype = 'uint32')
        #hdr.create_dataset('spatial/npars', data = self.spatial_npars, dtype = 'uint32')
        
        #hdr.create_dataset('temporal/niter', data = self.temporal_niter, dtype = 'uint32')
        #hdr.create_dataset('temporal/chisq', data = self.temporal_chisq, dtype = 'float64')
        #hdr.create_dataset('temporal/npoints', data = self.temporal_npoints, dtype = 'uint32')
        #hdr.create_dataset('temporal/npars', data = self.temporal_npars, dtype = 'uint32')
        
        # Write the data.
        grp = f.create_group('data')
        
        # Write the magnitudes.
        grp.create_dataset('magnitudes/ascc', data = ascc)
        grp.create_dataset('magnitudes/vmag', data = vmag, dtype = 'float32')
        #grp.create_dataset('magnitudes/nobs', data = self.nobs_m, dtype = 'uint32')
        grp.create_dataset('magnitudes/mag', data = vmag, dtype = 'float32')
        grp.create_dataset('magnitudes/sigma', data = err1, dtype = 'float32')
        
        # Write the camera transmission.
        idx1, idx2 = np.where(~np.isnan(cam))
        grp.create_dataset('trans/idx1', data = idx1, dtype = 'uint32')
        grp.create_dataset('trans/idx2', data = idx2, dtype = 'uint32')
        grp.create_dataset('trans/nobs', data = nobs_cam[idx1, idx2], dtype = 'uint32')
        grp.create_dataset('trans/trans', data = cam[idx1, idx2], dtype = 'float32')
        
        grp['trans'].attrs['grid'] = 'polar'
        grp['trans'].attrs['nx'] = 13500
        grp['trans'].attrs['ny'] = 722
        
        # Write the intrapixel variations.
        idx1, idx2 = np.where(~np.isnan(A[:,:,0]))
        grp.create_dataset('intrapix/idx1', data = idx1, dtype = 'uint32')
        grp.create_dataset('intrapix/idx2', data = idx2, dtype = 'uint32')
        grp.create_dataset('intrapix/nobs', data = nobs_A[idx1, idx2], dtype = 'uint32')
        grp.create_dataset('intrapix/sinx', data = A[idx1, idx2, 0], dtype = 'float32')
        grp.create_dataset('intrapix/cosx', data = A[idx1, idx2, 1], dtype = 'float32')
        grp.create_dataset('intrapix/siny', data = A[idx1, idx2, 2], dtype = 'float32')
        grp.create_dataset('intrapix/cosy', data = A[idx1, idx2, 3], dtype = 'float32')
        
        grp['intrapix'].attrs['grid'] = 'polar'
        grp['intrapix'].attrs['nx'] = 270
        grp['intrapix'].attrs['ny'] = 722
        
        # Write the sky transmission.
        idx, lstseq = np.where(~np.isnan(clouds))
        grp.create_dataset('clouds/idx', data = idx, dtype = 'uint32')
        grp.create_dataset('clouds/lstseq', data = lstseq + lstmin, dtype = 'uint32')
        grp.create_dataset('clouds/nobs', data = nobs_clouds[idx, lstseq], dtype = 'uint32')
        grp.create_dataset('clouds/clouds', data = clouds[idx, lstseq], dtype = 'float32')
        grp.create_dataset('clouds/sigma', data = err3[idx, lstseq], dtype = 'float32')
        
        grp['clouds'].attrs['grid'] = 'healpix'
        grp['clouds'].attrs['nx'] = 23
        grp['clouds'].attrs['lstmin'] = lstmin
        grp['clouds'].attrs['lstmax'] = lstmax
        grp['clouds'].attrs['lstlen'] = lstlen
    
    #cam = cam.reshape((722, 13502))
    #cam = np.roll(cam, 13502/2, axis=1)
    #clouds = clouds.reshape((Ntot, lstlen))
    #err3 = err3.reshape((Ntot, lstlen))
    #A = A.reshape((272, 722, 4))

    #for i in range(722):
    
        #if np.any(np.isfinite(cam[i])):
    
            #plt.plot(cam[i], '.', c='k')
            #plt.show()

    #vmin = np.nanpercentile(cam, 1)
    #vmax = np.nanpercentile(cam, 99)

    #plt.imshow(cam, aspect='auto', interpolation='nearest', cmap=viridis, vmin=vmin, vmax=vmax)
    #plt.colorbar()
    #plt.show()
    
    #ax = plt.subplot(221)
    #plt.imshow(A[:,:,0].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    #plt.colorbar()
    
    #plt.subplot(222, sharex=ax, sharey=ax)
    #plt.imshow(A[:,:,1].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    #plt.colorbar()
    
    #plt.subplot(223, sharex=ax, sharey=ax)
    #plt.imshow(A[:,:,2].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    #plt.colorbar()
    
    #plt.subplot(224, sharex=ax, sharey=ax)
    #plt.imshow(A[:,:,3].T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1, cmap=viridis)
    #plt.colorbar()
    
    #plt.show()
    
    #ax = plt.subplot(211)
    #plt.imshow(clouds, aspect='auto', interpolation='nearest', vmin=-.5, vmax=.5, cmap=viridis)
    #plt.colorbar()

    #plt.subplot(212, sharex=ax, sharey=ax)
    #plt.imshow(err3, aspect='auto', interpolation='nearest', vmin=0, vmax=.2, cmap=viridis)
    #plt.colorbar()
    #plt.show()

    return

if __name__ == '__main__':
    main()
