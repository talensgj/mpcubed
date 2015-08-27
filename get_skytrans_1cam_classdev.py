#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from index_functions import index_statistics
from coordinate_grids import HealpixGrid, PolarGrid, CartesianGrid
from sysrem import sysrem

import healpy

class SkyTransmission():
    
    def __init__(self, grid='healpix', nx=8, ny=0, margin=0):
        
        self.grid = grid
        self.nx = nx
        self.ny = ny
        self.margin = margin

        return 


    def calculate(self, fLCfile, redfile, skyfile=None):
        
        self.fLCfile = fLCfile
        self.redfile = redfile
        
        # Filename of output file.
        if skyfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'sky_'+tail.rsplit('_')[-1]
            self.skyfile = os.path.join(head, tail)
        else:
            self.skyfile = skyfile
        
        # Read the stellar header information.
        with h5py.File(self.fLCfile, 'r') as f:
            
            hdr = f['table_header']
            self.ascc = hdr['ascc'].value
            self.vmag = hdr['vmag'].value
            self.ra = hdr['ra'].value
            self.dec = hdr['dec'].value
            self.nobs = hdr['nobs'].value.astype('int')
        
        # Calculate the transmission.
        if self.grid == 'healpix':
            self._calculate_healpix()
          
        elif self.grid == 'cartesian':
            self._calculate_cartesian()
             
        elif self.grid in ['polar']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid

        return

    def _read_data(self, ascc, nobs, readxy=False):
        nstars = len(ascc)
        ndata = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        lstidx = np.zeros(ndata)
        sky = np.zeros(ndata)
        flags1 = np.zeros(ndata)
        
        if readxy:
            x = np.zeros(ndata)
            y = np.zeros(ndata)
        
        cflux0 = np.zeros(ndata)
        ecflux0 = np.zeros(ndata)
        flags2 = np.zeros(ndata) 
        
        with h5py.File(self.fLCfile, 'r') as f, h5py.File(self.redfile, 'r') as g:
        
            lc = f['data']
            rc = g['data']
        
            for i in range(nstars):
                
                lstidx[select[i]:select[i+1]] = lc[ascc[i]]['lstidx']
                sky[select[i]:select[i+1]] = lc[ascc[i]]['sky']
                flags1[select[i]:select[i+1]] = lc[ascc[i]]['flag']
                
                if readxy:
                    x[select[i]:select[i+1]] = lc[ascc[i]]['x']
                    y[select[i]:select[i+1]] = lc[ascc[i]]['y']
                
                cflux0[select[i]:select[i+1]] = rc[ascc[i]]['cflux0']
                ecflux0[select[i]:select[i+1]] = rc[ascc[i]]['ecflux0']
                flags2[select[i]:select[i+1]] = rc[ascc[i]]['flags']
        
        lstidx = lstidx.astype('int')
   
        if readxy:
            return lstidx, sky, flags1, cflux0, ecflux0, flags2, x, y
        else:
            return lstidx, sky, flags1, cflux0, ecflux0, flags2

    def _calculate_healpix(self):
        
        # Healpix grid instance.
        hg = HealpixGrid(self.nx)
        
        # Assign stars to sky bins and count the number of stars in each bin.
        skyidx, skyuni = hg.find_gridpoint(self.ra, self.dec, compact=True)
        starcount = index_statistics(skyuni, skyuni, statistic='count')
        
        # Create arrays.
        niter = np.zeros(len(skyidx), dtype='int')
        chisq = np.zeros(len(skyidx), dtype='float')
        npoints = np.zeros(len(skyidx), dtype='int')
        npars = np.zeros(len(skyidx), dtype='int')
        
        for ind in range(len(skyidx)):
            
            # Select stars in the current sky bin.
            here = (skyuni == ind)
            ascc = self.ascc[here]
            vmag = self.vmag[here]
            nobs = self.nobs[here]
            
            # Read data for these stars.
            lstidx, sky, flags1, cflux0, ecflux0, flags2 = self._read_data(ascc, nobs)
            
            # Create the staridx
            staridx = np.repeat(np.arange(len(ascc)), nobs)
            
            # Remove bad datapoints.
            here = (cflux0 > 0)&(ecflux0 > 0)&(sky > 0)&(flags1 < 1)&(flags2 < 1)
            cflux0 = cflux0[here]
            ecflux0 = ecflux0[here]
            lstidx = lstidx[here]
            staridx = staridx[here]
            
            # If no good data skip this bin. I don't like skipping bins, but if I do I may as well do it better than this.
            if len(cflux0) == 0:
                print 'No good data points in this bin.'
                continue
            
            # Make the lstidx ascending from 0 and count the number of datapoints at each lstidx.
            lstidx, lstuni = np.unique(lstidx, return_inverse=True)
            pointcount = index_statistics(lstuni, lstuni, statistic='count')
            
            # Compute the sky transmission curve.
            skytrans, flux, niter[ind], chisq[ind], chisq_sky, chisq_flux, npoints[ind], npars[ind] = sysrem(lstuni, staridx, cflux0, ecflux0, a2 = (1e7)*10**(vmag/-2.5))
        
            with h5py.File(self.skyfile) as f:
                
                grp = f.create_group('data/%i'%skyidx[ind])
                grp.create_dataset('lstidx', data=lstidx)
                grp.create_dataset('pointcount', data=pointcount)
                grp.create_dataset('skytrans', data=skytrans)
                grp.create_dataset('chisq_sky', data=chisq_sky)
        
        with h5py.File(self.skyfile) as f:
            
            grp = f.create_group('header')
            
            grp.attrs['grid'] = self.grid
            grp.attrs['nx'] = self.nx
            grp.attrs['ny'] = self.ny
            grp.attrs['margin'] = self.margin
            
            grp.create_dataset('skyidx', data = skyidx)
            grp.create_dataset('starcount', data = starcount)
            grp.create_dataset('niter', data = niter)
            grp.create_dataset('chisq', data = chisq)
            grp.create_dataset('npoints', data = npoints)
            grp.create_dataset('npars', data = npars)
            
            ### NEEDS TO BE INSIDE THE LOOP?
            #stars = np.unique(star_id)
            #grp.create_dataset('ascc', data = ascc[stars])
            #grp.create_dataset('vmag', data = vmag[stars])
            #grp.create_dataset('flux', data = a1[stars])
            #grp.create_dataset('dec', data = dec[stars])
            #grp.create_dataset('chisq_flux', data = chisq_pbin1[stars])

    def _calculate_cartesian(self):

        # Cartesian grid instance.
        cg = CartesianGrid(self.nx, self.ny, margin=self.margin)
        
        # Read data for these stars.
        lstidx, sky, flags1, cflux0, ecflux0, flags2, x, y = self._read_data(self.ascc, self.nobs, readxy=True)
            
        # Create the staridx
        staridx = np.repeat(np.arange(len(self.ascc)), self.nobs)
            
        # Remove bad datapoints.
        here = (cflux0 > 0)&(ecflux0 > 0)&(sky > 0)&(flags1 < 1)&(flags2 < 1)
        cflux0 = cflux0[here]
        ecflux0 = ecflux0[here]
        lstidx = lstidx[here]
        staridx = staridx[here]
        x = x[here]
        y = y[here]
        
        skyidx, skyuni = cg.find_gridpoint(x, y, compact=True)
        lstidx, lstuni = np.unique(lstidx, return_inverse=True)
            
        idx = np.ravel_multi_index([skyuni, lstuni], (len(skyidx), len(lstidx)))
        pointcount = np.bincount(idx)
            
        idx, idxuni = np.unique(idx, return_inverse=True)
            
        # Compute the sky transmission curve.
        skytrans, flux, niter, chisq, chisq_sky, chisq_flux, npoints, npars = sysrem(idxuni, staridx, cflux0, ecflux0, a2 = (1e7)*10**(self.vmag/-2.5))
        
        array = np.full((len(skyidx), len(lstidx)), fill_value=np.nan)
        array[np.unravel_index(idx, (len(skyidx), len(lstidx)))] = skytrans
        
        for ind in range(len(skyidx)):
        
            with h5py.File(self.skyfile) as f:
                
                grp = f.create_group('data/%i'%skyidx[ind])
                grp.create_dataset('lstidx', data=lstidx)
                #grp.create_dataset('pointcount', data=pointcount)
                grp.create_dataset('skytrans', data=array[ind])
                #grp.create_dataset('chisq_sky', data=chisq_sky)
        
        with h5py.File(self.skyfile) as f:
            
            grp = f.create_group('header')
            
            grp.attrs['grid'] = self.grid
            grp.attrs['nx'] = self.nx
            grp.attrs['ny'] = self.ny
            grp.attrs['margin'] = self.margin
            
            grp.create_dataset('skyidx', data = skyidx)
            #grp.create_dataset('starcount', data = starcount)
            grp.create_dataset('niter', data = niter)
            grp.create_dataset('chisq', data = chisq)
            grp.create_dataset('npoints', data = npoints)
            grp.create_dataset('npars', data = npars)
            
            ### NEEDS TO BE INSIDE THE LOOP?
            #stars = np.unique(star_id)
            #grp.create_dataset('ascc', data = ascc[stars])
            #grp.create_dataset('vmag', data = vmag[stars])
            #grp.create_dataset('flux', data = a1[stars])
            #grp.create_dataset('dec', data = dec[stars])
            #grp.create_dataset('chisq_flux', data = chisq_pbin1[stars])

        return
        
class SkyFile():
    
    def __init__(self, skyfile):
        
        self.skyfile = skyfile
        
        return
        
    def _read_skyfile(self):
        
        with h5py.File(self.skyfile, 'r') as f:
            
            self.grid = f['header'].attrs['grid']
            self.nx = f['header'].attrs['nx']
            self.ny = f['header'].attrs['ny']
            self.margin = f['header'].attrs['margin']
            
            npix = healpy.nside2npix(self.nx)
            
            self.starcount = np.full(npix, fill_value=np.nan)
            self.skytrans = np.full((npix, 13500), fill_value=np.nan)
            self.pointcount = np.full((npix, 13500), fill_value=np.nan)
            
            skyidx = f['header/skyidx'].value
            self.starcount[skyidx] = f['header/starcount'].value
            
            for ind in skyidx:
                
                data = f['data/%i'%ind]
                lstidx = data['lstidx'].value
                
                self.skytrans[ind, lstidx] = data['skytrans'].value
                self.pointcount[ind, lstidx] = data['pointcount'].value
            
        return
        
    def visualize(self):
        
        self._read_skyfile()
            
        if self.grid == 'healpix':
            self._visualize_healpix()
            
        elif self.grid in ['polar', 'cartesian']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid
            exit()
            
        return

    def _visualize_healpix(self):
        
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from viridis import viridis 
        
        rcParams['xtick.labelsize'] = 'large'
        rcParams['ytick.labelsize'] = 'large'
        rcParams['axes.labelsize'] = 'x-large'
        rcParams['image.interpolation'] = 'none'
        rcParams['image.origin'] = 'lower'
        
        xlim, ylim = np.where(np.isfinite(self.skytrans))
        
        array = self.skytrans
        array[self.pointcount<=10] = np.nan
        array = array[np.unique(xlim)]
        array = array[:,np.unique(ylim)]
        
        # Figure showing the transmission map.
        plt.imshow(array/np.nanmedian(array, axis=1, keepdims=True), aspect='auto', cmap=viridis, vmin=.8, vmax=1.2)
        cb = plt.colorbar()
        plt.xlabel('LST [idx]')
        plt.ylabel('Sky [idx]')
        cb.set_label('Sky')
        plt.show()
        
        hg = HealpixGrid(self.nx)
        
        for i in np.unique(ylim):
            print i
            array = self.skytrans
            array = array/np.nanmedian(array, axis=1, keepdims=True)
        
            array = hg.put_values_on_grid(array[:,i])
            
            healpy.mollview(array, min=.8, max=1.2, cmap=viridis, title='LSTIDX = %i'%i, unit='Sky')
            healpy.graticule()
            plt.savefig('/data2/talens/Jul2015/Sky/sky_%.5i.png'%i)
            plt.close()
        
        return
        
#test = SkyTransmission()
#test.calculate('/data2/talens/Jul2015/fLC_20150714LPE.hdf5', '/data2/talens/Jul2015/redtest.hdf5', skyfile='/data2/talens/Jul2015/skytest.hdf5')

test = SkyTransmission(grid='cartesian', nx=15, ny=10, margin=45)
test.calculate('/data2/talens/Jul2015/fLC_20150714LPE.hdf5', '/data2/talens/Jul2015/redtest.hdf5', skyfile='/data2/talens/Jul2015/skytest_cart.hdf5')

#test = SkyFile('/data2/talens/Jul2015/skytest.hdf5')
#test.visualize()










