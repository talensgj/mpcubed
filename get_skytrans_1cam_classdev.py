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
    
    def __init__(self):
        self.grid = 'healpix'
        self.nx = 8
        self.ny = 0
        self.margin = 0

        return 


    def get_gridparameters(self):
        
        print 'grid = %s'%self.grid
        print 'nx = %i'%self.nx
        print 'ny = %i'%self.ny
        print 'margin = %i'%self.margin
        
        return 0
        
        
    def set_gridparameters(self, grid, nx, ny, margin):
        
        self.grid = grid
        self.nx = nx
        self.ny = ny
        self.margin = margin
        
        return 0


    def calculate(self, fLC, red):
        
        self.fLC = fLC
        self.red = red
        
        # Read the stellar header information.
        with h5py.File(self.fLC, 'r') as f:
            
            hdr = f['table_header']
            self.ascc = hdr['ascc'].value
            self.vmag = hdr['vmag'].value
            self.ra = hdr['ra'].value
            self.dec = hdr['dec'].value
            self.nobs = hdr['nobs'].value.astype('int')
        
        if self.grid == 'healpix':
            self.process_healpix()
             
        elif self.grid in ['polar', 'cartesian']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid

        return 0

    def read_data(self, ascc, nobs):
        nstars = len(ascc)
        ndata = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        lstidx = np.zeros(ndata)
        sky = np.zeros(ndata)
        flags1 = np.zeros(ndata)
        
        cflux0 = np.zeros(ndata)
        ecflux0 = np.zeros(ndata)
        flags2 = np.zeros(ndata) 
        
        with h5py.File(self.fLC, 'r') as f, h5py.File(self.red, 'r') as g:
        
            lc = f['data']
            rc = g['data']
        
            for i in range(nstars):
                
                lstidx[select[i]:select[i+1]] = lc[ascc[i]]['lstidx']
                sky[select[i]:select[i+1]] = lc[ascc[i]]['sky']
                flags1[select[i]:select[i+1]] = lc[ascc[i]]['flag']
                
                cflux0[select[i]:select[i+1]] = rc[ascc[i]]['cflux0']
                ecflux0[select[i]:select[i+1]] = rc[ascc[i]]['ecflux0']
                flags2[select[i]:select[i+1]] = rc[ascc[i]]['flags']
        
        lstidx = lstidx.astype('int')
   
        return lstidx, sky, flags1, cflux0, ecflux0, flags2


    def process_healpix(self):
        
        # Create the skygrid and assign the stars to the grid.
        grid = HealpixGrid(self.nx)
        bins, binnum = grid.find_gridpoint(self.ra, self.dec, compact=True)
        
        # Count the number of stars in each sky bin.
        starcount = index_statistics(binnum, binnum, statistic='count')
        
        niter = np.zeros(len(bins), dtype='int')
        chisq = np.zeros(len(bins), dtype='float')
        npoints = np.zeros(len(bins), dtype='int')
        npars = np.zeros(len(bins), dtype='int')
        
        head, tail = os.path.split(self.fLC)
        tail = 'sky_'+tail.rsplit('_')[-1]
        skyf = os.path.join(head, tail)
        
        for ind in range(len(bins)):
            
            # Select stars in the current sky bin.
            here = (binnum == ind)
            ascc = self.ascc[here]
            vmag = self.vmag[here]
            nobs = self.nobs[here]
            
            # Read data for these stars.
            lstidx, sky, flags1, cflux0, ecflux0, flags2 = self.read_data(ascc, nobs)
            
            # Create a stellar index.
            star_id = np.repeat(np.arange(len(ascc)), nobs)
            
            # Remove bad datapoints.
            here = (flags1<1)&(flags2<1)&(cflux0>0)&(sky>0)&(ecflux0>0)
            lstidx = lstidx[here]
            cflux0 = cflux0[here]
            ecflux0 = ecflux0[here]
            star_id = star_id[here]
            
            # If no good data skip this bin. I don't like skipping bins, but if I do I may as well do it better than this.
            if len(cflux0) == 0:
                continue
            
            # Make the lstidx ascending from 0.
            times, lstidx = np.unique(lstidx, return_inverse=True)
            
            # Count the number of good datapoints at each time stamp for this sky patch.
            count = index_statistics(lstidx, lstidx, statistic='count')
            
            # Compute the sky transmission curve.
            a2, a1, niter[ind], chisq[ind], chisq_pbin2, chisq_pbin1, npoints[ind], npars[ind] = sysrem(lstidx, star_id, cflux0, ecflux0, a2 = (1e7)*10**(vmag/-2.5))
        
            with h5py.File(skyf) as f:
                
                grp = f.create_group('data/%i'%bins[ind])
                grp.create_dataset('lstidx', data=times)
                grp.create_dataset('count', data=count)
                grp.create_dataset('sky', data=a2)
                grp.create_dataset('chisq_sky', data=chisq_pbin2)
        
        with h5py.File(skyf) as f:
            
            grp = f.create_group('header')
            
            grp.attrs['grid'] = self.grid
            grp.attrs['nx'] = self.nx
            grp.attrs['ny'] = self.ny
            grp.attrs['margin'] = self.margin
            
            grp.create_dataset('binnum', data = bins)
            grp.create_dataset('nstars', data = starcount)
            grp.create_dataset('niter', data = niter)
            grp.create_dataset('chisq', data = chisq)
            grp.create_dataset('npoints', data = npoints)
            grp.create_dataset('npars', data = npars)
            
            #grp.create_dataset('ascc', data=ascc[np.unique(star_id)])
            #grp.create_dataset('flux', data=a1)
            #grp.create_dataset('chisq_flux', data=chisq_pbin1)

        return 0

test = SkyTransmission()
test.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/red_20150716LPC.hdf5')
