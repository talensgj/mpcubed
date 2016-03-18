#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package import misc
from package import IO
from package.coordinates import grids

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
    
    skygrid = grids.PolarGrid(36, 36)
    raidx, decidx = skygrid.radec2idx(ra, dec)
    print np.unique(decidx)
    select = (decidx == 25)
    ascc = ascc[select]
    ra = ra[select]
    dec = dec[select]
    vmag = vmag[select]
    nobs = nobs[select]
    raidx = raidx[select]
    staridx = staridx[select]

    lstseq, lst, flux, eflux, sky, flag, x, y = f.read_data(['lstseq', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], ascc, nobs)
    lstseq = lstseq.astype('int')
    
    ra = np.repeat(ra, nobs)
    dec = np.repeat(dec, nobs)
    raidx = np.repeat(raidx, nobs)  
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
    raidx = raidx[select]
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
    
    skyidx = np.ravel_multi_index((raidx, lstseq), (38, lstlen)) 
    camidx = np.ravel_multi_index((haidx1, decidx), (13502, 722))
    idx3 = np.ravel_multi_index((haidx2, decidx), (272, 722))

    sky = np.zeros(np.amax(skyidx) + 1)
    ipx_x = 0
    ipx_y = 0
    b = np.zeros(np.amax(idx3) + 1)
    d = np.zeros(np.amax(idx3) + 1)
    err2 = np.zeros(np.amax(staridx) + 1)
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
    weights = 1/emag**2
    for niter in range(50):
        
        cam = np.bincount(camidx, weights*(mag - sky[skyidx] - ipx_x - ipx_y))/np.bincount(camidx, weights)
            
        sky = np.bincount(skyidx, weights*(mag - cam[camidx] - ipx_x - ipx_y))/np.bincount(skyidx, weights)
    
        a = np.bincount(idx3, weights*(mag - cam[camidx] - sky[skyidx] - b[idx3]*csx - ipx_y)*snx)/np.bincount(idx3, weights*snx**2)
        b = np.bincount(idx3, weights*(mag - cam[camidx] - sky[skyidx] - a[idx3]*snx - ipx_y)*csx)/np.bincount(idx3, weights*csx**2)
        
        ipx_x = a[idx3]*snx + b[idx3]*csx
        
        c = np.bincount(idx3, weights*(mag - cam[camidx] - sky[skyidx] - ipx_x - d[idx3]*csy)*sny)/np.bincount(idx3, weights*sny**2)
        d = np.bincount(idx3, weights*(mag - cam[camidx] - sky[skyidx] - ipx_x - c[idx3]*sny)*csy)/np.bincount(idx3, weights*csy**2)
        
        ipx_y = c[idx3]*sny + d[idx3]*csy
    
        
        #err1 = find_sigma(skyidx, mag - cam[camidx] - sky[skyidx] - ipx_x - ipx_y, emag)
        #weights = 1/(emag**2 + (err1**2)[skyidx])
    
        err1 = find_sigma(skyidx, mag - cam[camidx] - sky[skyidx] - ipx_x - ipx_y, np.sqrt(emag**2 + err2[staridx]**2))
        err2 = find_sigma(staridx, mag - cam[camidx] - sky[skyidx] - ipx_x - ipx_y, np.sqrt(emag**2 + err1[skyidx]**2))
        
        weights = 1/(emag**2 + (err1**2)[skyidx] + (err2**2)[staridx])
    
    haidx, decidx = np.unravel_index(np.arange(len(cam)), (13502, 722))
    with h5py.File('/data2/talens/polar_eqarea.hdf5') as f:
        grp = f.create_group('data')
        grp.create_dataset('haidx', data=haidx)
        grp.create_dataset('decidx', data=decidx)
        grp.create_dataset('value', data=cam)
    
    #fit = (cam[camidx] + sky[skyidx] + ipx_x + ipx_y)
    #for i in range(0, 100000, 10000):
        #plt.plot(mag[i:i+10000], '.')
        #plt.plot(fit[i:i+10000])
        #plt.show()
        
    #array = np.full((13502, 722), fill_value=np.nan)
    #array[haidx1, decidx] = cam[camidx]
    #plt.imshow(array.T, aspect='auto', interpolation='nearest')
    #plt.colorbar()
    #plt.show()
    
    #array = np.full((272, 722), fill_value=np.nan)
    #array[haidx2, decidx] = c[idx3]
    #plt.imshow(array.T, aspect='auto', interpolation='nearest', vmin=-.1, vmax=.1)
    #plt.colorbar()
    #plt.show()
    
    #ax = plt.subplot(121)
    #array = np.full((38, lstlen), fill_value=np.nan)
    #array[raidx, lstseq] = sky[skyidx]
    #plt.imshow(array.T, aspect='auto', interpolation='nearest', vmin=-.5, vmax=.5)
    #plt.colorbar()

    #plt.subplot(122, sharex=ax, sharey=ax)
    #array = np.full((38, lstlen), fill_value=np.nan)
    #array[raidx, lstseq] = err1[skyidx]
    #plt.imshow(array.T, aspect='auto', interpolation='nearest', vmin=0, vmax=.2)
    #plt.colorbar()
    #plt.show()

    return

if __name__ == '__main__':
    main()
