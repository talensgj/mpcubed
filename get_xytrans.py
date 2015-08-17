#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob

import h5py
import numpy as np

import healpy

import matplotlib.pyplot as plt
from viridis import viridis

from index_functions import index_statistics
from coordinate_grids import HealpixGrid, PolarGrid, CartesianGrid
from sysrem import sysrem

def transmission(filename):

    with h5py.File(filename) as f:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        vmag = hdr['vmag'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        
        select = np.append(0, np.cumsum(Nobs))
    
        lst = np.zeros(np.sum(Nobs))
        flux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        sky = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        x = np.zeros(np.sum(Nobs))
        y = np.zeros(np.sum(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            lst[select[i]:select[i+1]] = data[ascc[i]]['lst']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            #eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            eflux0[select[i]:select[i+1]] = index_statistics(data[ascc[i]]['lstidx']//50, data[ascc[i]]['flux0'], statistic='std', keeplength=True)
            sky[select[i]:select[i+1]] = data[ascc[i]]['sky']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
            x[select[i]:select[i+1]] = data[ascc[i]]['x']
            y[select[i]:select[i+1]] = data[ascc[i]]['y']
    
    dec_1 = np.copy(dec)
    
    print 'Done reading.'

    # Sky coordinates of data.
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    # Create the indices.
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    
    # Remove bad data.
    here, = np.where((flags < 1)&(flux0>0)&(sky>0)&(eflux0>0))
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    ha = ha[here]
    dec = dec[here]
    x = x[here]
    y = y[here]
    star_id = star_id[here]

    #hg = PolarGrid(13500, 720)
    #bins, binnum = hg.find_gridpoint(ha, dec, compact=True)
    
    hg = CartesianGrid(800, 600, margin=50)
    bins, binnum = hg.find_gridpoint(x, y, compact=True)
    
    #hg = HealpixGrid(512)
    #bins, binnum = hg.find_gridpoint(ha, dec, compact=True)
    
    count = index_statistics(binnum, flux0, statistic='count', keeplength=False)
    print np.percentile(count, 5), np.percentile(count, 95)
    
    # Compute the transmission using sysrem.
    a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(binnum, star_id, flux0, eflux0, a2=1e7*10**(vmag/-2.5))

    with h5py.File('/data2/talens/Jul2015/Trans0716LPS_cg800x600m50.hdf5.hdf5') as f:
        f.create_dataset('binnum', data=bins)
        f.create_dataset('count', data=count)
        f.create_dataset('trans', data=a2)
    exit()
    trans = hg.put_values_on_grid(a2, ind=bins, fill_value=np.nan)
    count = hg.put_values_on_grid(count, ind=bins)
    chisq_map = hg.put_values_on_grid(chisq_pbin2, ind=bins, fill_value=np.nan)
    stars = np.unique(star_id)
    
    plt.subplot(111)
    plt.title(r'niter = %i, $\chi^2_\nu$ = %.2f, npoints = %i, npars = %i'%(niter, chisq, npoints, npars))
    plt.imshow(trans[1:-1,1:-1].T, interpolation='None', origin='lower', aspect='auto', extent=(0,360,-90,90), vmin=0, vmax=1.5, cmap=viridis)
    plt.colorbar().set_label('Transmission')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec_1), np.amax(dec_1))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    
    #plt.subplot(212)
    #arr = a1[stars]/(1e7*10**(vmag[stars]/-2.5))
    #plt.plot(dec_1[stars], arr, '.', alpha=.1)
    #plt.xlim(np.amin(dec_1), np.amax(dec_1))
    #plt.ylim(.8, 1.2)
    #plt.xlabel('Dec [deg]')
    #plt.ylabel('F/V')
    
    plt.show()
    
    exit()
    median = index_statistic(idx[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), statistic='median', keeplength=False)
    
    plt.subplot(211)
    plt.semilogy(dec_1[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), '.', alpha=.1)
    plt.xlim(np.amin(dec), np.amax(dec))
    plt.xlabel('Dec [deg]')
    plt.ylabel('F/V')
    
    a1[stars] = a1[stars]/median[idx[stars]]
    
    plt.subplot(212)
    plt.semilogy(dec_1[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), '.', alpha=.1)
    plt.xlim(np.amin(dec), np.amax(dec))
    plt.xlabel('Dec [deg]')
    plt.ylabel('F/V')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    trans[:,np.unique(idx[stars])] = trans[:,np.unique(idx[stars])]*median[np.unique(idx[stars])]
    
    ax = plt.subplot(221)
    plt.imshow(trans[1:-1,1:-1].T, interpolation='None', aspect='auto', origin='lower', vmin=0, vmax=1.5, extent=(0,360,-90,90), cmap=test_cm)
    plt.colorbar().set_label('Transmission')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    
    plt.subplot(222, sharex=ax, sharey=ax)
    plt.imshow(count[1:-1,1:-1].T, interpolation='None', aspect='auto', origin='lower', extent=(0,360,-90,90), cmap=test_cm)
    plt.colorbar().set_label('# points')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    
    plt.subplot(223, sharex=ax, sharey=ax)
    arr = chisq_map[1:-1,1:-1].T/count[1:-1,1:-1].T
    plt.imshow(arr, interpolation='None', aspect='auto', origin='lower', vmin=0, vmax=np.percentile(arr[np.isfinite(arr)], 98), extent=(0,360,-90,90), cmap=test_cm)
    plt.colorbar().set_label(r'$\chi^2/N$')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    
    plt.subplot(224)
    arr = chisq_pbin1[stars]/np.bincount(star_id)[stars]
    plt.plot(vmag[stars], arr, '.', alpha=.1)
    plt.xlim(8.4, 2)
    plt.ylim(0, np.percentile(arr[np.isfinite(arr)], 98))
    plt.xlabel('V [mag]')
    plt.ylabel(r'$\chi^2/N$')
    
    plt.tight_layout()
    #plt.savefig('0203LPE_maps.png')
    plt.show()
    plt.close()
    
    cg = CartesianGrid(200, 150, margin=50)
    binnum = cg.find_gridpoint(x, y)
    a2 = np.bincount(binnum, flux0*a1[star_id]/eflux0**2)/np.bincount(binnum, (a1**2)[star_id]/eflux0**2)
    trans = cg.put_values_on_grid(a2)
    
    count = np.bincount(binnum)
    count = cg.put_values_on_grid(count)
    
    chisq_map = np.bincount(binnum, (flux0-a1[star_id]*a2[binnum])**2/eflux0**2)
    chisq_map = cg.put_values_on_grid(chisq_map)
    
    chisq_stars = np.bincount(binnum, (flux0-a1[star_id]*a2[binnum])**2/eflux0**2)
    
    ax = plt.subplot(221)
    trans = trans/np.percentile(trans[1:-1,1:-1], 99)
    plt.imshow(trans[1:-1,1:-1].T, interpolation='None', origin='lower', aspect='equal', vmin=0, vmax=1, extent=(50,4008-50,50,2672-50), cmap=test_cm)
    plt.colorbar().set_label('Transmission')
    plt.xlim(0,4008)
    plt.ylim(0,2672)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')
    
    plt.subplot(222, sharex=ax, sharey=ax)
    plt.imshow(count[1:-1,1:-1].T, interpolation='None', origin='lower', aspect='equal', extent=(50,4008-50,50,2672-50), cmap=test_cm)
    plt.colorbar().set_label('# points')
    plt.xlim(0,4008)
    plt.ylim(0,2672)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')
    
    plt.subplot(223, sharex=ax, sharey=ax)
    arr = chisq_map[1:-1,1:-1].T/count[1:-1,1:-1].T
    plt.imshow(arr, interpolation='None', origin='lower', aspect='equal', vmin=0, vmax=np.percentile(arr[np.isfinite(arr)], 98), extent=(50,4008-50,50,2672-50), cmap=test_cm)
    plt.colorbar().set_label(r'$\chi^2/N$')
    plt.xlim(0,4008)
    plt.ylim(0,2672)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')
    
    plt.subplot(224)
    arr = chisq_stars[stars]/np.bincount(star_id)[stars]
    plt.plot(vmag[stars], arr, '.', alpha=.1)
    plt.xlim(8.4, 2)
    plt.ylim(0, np.percentile(arr[np.isfinite(arr)], 98))
    plt.xlabel('V [mag]')
    plt.ylabel(r'$\chi^2/N$')
    
    plt.tight_layout()
    #plt.savefig('0203LPE_maps.png')
    plt.show()
    
    
    
    
    return 0

if __name__ == '__main__':
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Compute the transmission map for a particular night and camera.')
    
    #parser.add_argument('path', type=str, help='the global path for search')
    #parser.add_argument('-n', '--night', default='', type=str, help='the night to reduce')
    #parser.add_argument('-c', '--camera',default='', type=str, help='the camera to reduce')
    
    #args = parser.parse_args()
    
    filelist = glob.glob('/data2/talens/Jul2015/fLC_20150716LPS.hdf5')
    #filelist = np.sort(filelist)

    #filelist = ['/data2/talens/fLC_20150203LPC.hdf5']

    for filename in filelist:
        print 'Data:', filename
        transmission(filename)
