#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from viridis import viridis

from coordinate_grids import PolarGrid, CartesianGrid, HealpixGrid

def plot_transmission(filename):
    
    with h5py.File(filename, 'r') as f:
        
        gridname = f['header'].attrs['grid']
        
        if gridname == 'polar':
            grid = PolarGrid(f['header'].attrs['nx'], f['header'].attrs['ny'])
        elif gridname == 'cartesian':
            #margin = f['Header'].attrs['margin']
            margin = 50
            grid = CartesianGrid(f['header'].attrs['nx'], f['header'].attrs['ny'], margin=50)
            
        elif gridname == 'healpix':
            grid = HealpixGrid(f['header'].attrs['nx'])
        else:
            print 'Unkown grid: %s'%gridname
            exit()
        
        binnum = f['data/binnum'].value
        trans = f['data/trans'].value
    
    trans = grid.put_values_on_grid(trans, binnum, fill_value=np.nan)
    
    if gridname == 'polar':
        
        trans = trans[1:-1,1:-1].T
        
        im = plt.imshow(trans, origin='lower', aspect='auto', interpolation='None', extent=(0,360,-90,90), cmap=viridis, vmin=0, vmax=1.5)
        cbar = plt.colorbar(im)
        cbar.set_label('Transmission', size='x-large')
        cbar.ax.tick_params(labelsize='large')
        
        plt.xticks(size='large')
        plt.yticks(size='large')
        plt.xlabel('HA [deg]', size='x-large')
        plt.ylabel('Dec [deg]', size='x-large')
        plt.tight_layout()
    
    elif gridname == 'cartesian':
    
        trans = trans[1:-1,1:-1].T
    
        im = plt.imshow(trans, origin='lower', interpolation='None', extent=(margin,4008-margin,margin,2672-margin), cmap=viridis, vmin=0, vmax=1.5)
        plt.xlim(0,4008)
        plt.ylim(0,2672)
        cbar = plt.colorbar(im)
        cbar.set_label('Transmission', size='x-large')
        cbar.ax.tick_params(labelsize='large')
    
        plt.xticks(size='large')
        plt.yticks(size='large')
        plt.xlabel('HA [deg]', size='x-large')
        plt.ylabel('Dec [deg]', size='x-large')
        plt.tight_layout()
    
    else:
        
        import healpy
        
        healpy.mollview(trans, min=0, max=1.5, cmap=viridis, unit='Transmission', xsize=3200)
        healpy.graticule()
        
    plt.show()
    
plot_transmission('/data2/talens/Jul2015/cam_20150716LPC_pg2700x720.hdf5')

with h5py.File('/data2/talens/Jul2015/cam_20150716LPC_pg2700x720.hdf5', 'r') as f:
    
    plt.subplot(211)
    plt.semilogy(f['data/vmag'], f['data/flux'], '.', alpha=.2, c='k')
    x = np.array([2,8.4])
    plt.semilogy(x, 1e7*10**(x/-2.5), c='r') 
    plt.xticks(size='large')
    plt.yticks(size='large')
    plt.xlabel('V [mag]', size='x-large')
    plt.ylabel('Flux [counts]', size='x-large') 
    plt.xlim(8.4,2)
    
    plt.subplot(212)
    plt.semilogy(f['data/dec'], f['data/flux'].value/(1e7*10**(f['data/vmag'].value/-2.5)), '.', alpha=.2, c='k')
    plt.xticks(size='large')
    plt.yticks(size='large')
    plt.xlabel('Dec [deg]', size='x-large')
    plt.ylabel('F/V', size='x-large')

plt.tight_layout()
plt.show()

#jd_ref = np.floor(lc['jdmid'])
   
#with h5py.File('/data2/talens/Jul2015/Trans0716LPN_pg2700x720.hdf5', 'r') as f:
    #bins = f['Data/binnum'].value
    #trans = f['Data/trans'].value

#pg = PolarGrid(2700, 720)
#trans_pg = pg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
#trans_pg = np.ravel(trans_pg)

#with h5py.File('/data2/talens/Jul2015/Trans0716LPN_cg400x300m50.hdf5', 'r') as f:
    #bins = f['Data/binnum'].value
    #trans = f['Data/trans'].value

#cg = CartesianGrid(400, 300, margin=50)
#trans_cg = cg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
#trans_cg = np.ravel(trans_cg)

#with h5py.File('/data2/talens/Jul2015/Trans0716LPN_hg256.hdf5', 'r') as f:
    #bins = f['Data/binnum'].value
    #trans = f['Data/trans'].value

#hg = HealpixGrid(256)
#trans_hg = hg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
        
#with h5py.File('/data2/talens/Jul2015/fLC_20150716LPN.hdf5', 'r') as f:
    #ascc = f['table_header/ascc'].value    
    #dec = f['table_header/dec'].value
    
    #here = dec>85
    #ascc = ascc[here]
    
    #for sid in ascc:
        #si = f['header/'+sid]
        #lc = f['data/'+sid]
    
        #ha = np.mod(lc['lst']*15.-np.repeat(si['ra'],si['nobs'].astype('int')), 360.)
        #dec = np.repeat(si['dec'], si['nobs'].astype('int'))
        #binnum_pg = pg.find_gridpoint(ha, dec)
        #binnum_hg = hg.find_gridpoint(ha, dec)
        #binnum_cg = cg.find_gridpoint(lc['x'], lc['y'])
    
        #jd_ref = np.floor(lc['jdmid'])
        
        #ax = plt.subplot(311)
        #plt.title('Dec = %.2f'%si['dec'])
        #plt.plot(lc['jdmid']-jd_ref, lc['flux0'], '.')
        #plt.xlabel('JD')
        #plt.ylabel('Flux')

        #plt.subplot(312, sharex=ax)
        #plt.plot(lc['jdmid']-jd_ref, trans_pg[binnum_pg], '.')
        #plt.plot(lc['jdmid']-jd_ref, trans_cg[binnum_cg], '.')
        #plt.plot(lc['jdmid']-jd_ref, trans_hg[binnum_hg], '.')
        #plt.ylim(0,1.5)
        #plt.xlabel('JD')
        #plt.ylabel('Trans')

        #plt.subplot(313, sharex=ax)
        #plt.plot(lc['jdmid']-jd_ref, lc['flux0']/trans_pg[binnum_pg], '.')
        #plt.plot(lc['jdmid']-jd_ref, lc['flux0']/trans_cg[binnum_cg], '.')
        #plt.plot(lc['jdmid']-jd_ref, lc['flux0']/trans_hg[binnum_hg], '.')
        #plt.xlim(.25,.75)
        #plt.xlabel('JD')
        #plt.ylabel('Trans')

        #plt.tight_layout()
        #plt.show()
        #plt.close()
