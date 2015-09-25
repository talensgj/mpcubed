#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from index_functions import index_statistics
from coordinate_grids import HealpixGrid, PolarGrid, CartesianGrid
from coarse_dev import coarse_decorrelation

import healpy

class SkyTransmission():
    
    def __init__(self, grid='healpix', nx=8, ny=0, margin=0):
        
        self.grid = grid
        self.nx = nx
        self.ny = ny
        self.margin = margin

        return 

    def calculate(self, fLCfile, redfile=None, skyfile=None):
        
        self.fLCfile = fLCfile
        
        # Filename of reduced data file.
        if redfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'coarsered_'+tail.rsplit('_')[-1]
            self.redfile = os.path.join(head, tail)
        else:
            self.redfile = redfile
        
        # Filename of output file.
        if skyfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'coarsesky_'+tail.rsplit('_')[-1]
            self.skyfile = os.path.join(head, tail)
        else:
            self.skyfile = skyfile
        
        # Read the stellar header information.
        with h5py.File(self.fLCfile, 'r') as f:
            
            hdr = f['header_table']
            self.ascc = hdr['ascc'].value
            self.vmag = hdr['vmag'].value
            self.ra = hdr['ra'].value
            self.dec = hdr['dec'].value
            self.nobs = hdr['nobs'].value.astype('int')
        
        # Calculate the transmission.
        if self.grid == 'healpix':
            self._calculate_healpix()
          
        elif self.grid == 'cartesian':
            print 'Grid %s not implemented yet.'%self.grid
            exit()
             
        elif self.grid in ['polar']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid

        return

    def _read_data(self, ascc, nobs):
        nstars = len(ascc)
        ndata = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        lstidx = np.zeros(ndata)
        sky = np.zeros(ndata)
        flags1 = np.zeros(ndata)

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
                
                cflux0[select[i]:select[i+1]] = rc[ascc[i]]['ipcflux0']
                ecflux0[select[i]:select[i+1]] = 2.5/np.log(10)*lc[ascc[i]]['eflux0']/lc[ascc[i]]['flux0']
                flags2[select[i]:select[i+1]] = rc[ascc[i]]['flags']
        
        lstidx = lstidx.astype('int')
   
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
            here = (ecflux0 > 0)&(sky > 0)&(flags1 < 1)&(flags2 < 1)
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
            magnitude, skytrans, sigma2, niter[ind], chisq[ind], npoints[ind], npars[ind] = coarse_decorrelation(staridx, lstuni, cflux0, ecflux0, verbose=True)
            
            with h5py.File(self.skyfile) as f:
                
                grp = f.create_group('data/%i'%skyidx[ind])
                grp.create_dataset('lstidx', data=lstidx)
                grp.create_dataset('pointcount', data=pointcount)
                grp.create_dataset('skytrans', data=skytrans)
                grp.create_dataset('sigma2', data=sigma2)
        
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
            
            self.starcount = np.full(npix, fill_value=0)
            self.skytrans = np.full((npix, 13500), fill_value=np.nan)
            self.clouds = np.full((npix, 13500), fill_value=np.nan)
            self.pointcount = np.full((npix, 13500), fill_value=0)
            #self.chisq_sky = np.full((npix, 13500), fill_value=np.nan)
            
            self.skyidx = f['header/skyidx'].value
            self.starcount[self.skyidx] = f['header/starcount'].value
            
            for idx in self.skyidx:
                
                try:
                    data = f['data/%i'%idx]
                except:
                    pass
                else:
                    lstidx = data['lstidx'].value
                
                    self.clouds[idx, lstidx] = data['sigma2'].value
                    self.skytrans[idx, lstidx] = data['skytrans'].value
                    self.pointcount[idx, lstidx] = data['pointcount'].value
                    #self.chisq_sky[idx, lstidx] = data['chisq_sky'].value
            
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
        
        self.skytrans[self.clouds>.05] = np.nan
        
        # Figure showing the transmission map.
        plt.imshow(self.skytrans, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
        cb = plt.colorbar()
        plt.xlabel('LST [idx]')
        plt.ylabel('Sky [idx]')
        plt.ylim(np.amin(self.skyidx)-.5,np.amax(self.skyidx)+.5) 
        cb.set_label('Sky')
        plt.show()
        
        return
        
    def correct(self, fLCfile=None, redfile=None):
        
        if fLCfile is None:
            head, tail = os.path.split(self.skyfile)
            tail = 'fLC_'+tail.rsplit('_')[-1]
            self.fLCfile = os.path.join(head, tail)
        else:
            self.fLCfile = fLCfile
            
        if redfile is None:
            head, tail = os.path.split(self.skyfile)
            tail = 'coarsered_'+tail.rsplit('_')[-1]
            self.redfile = os.path.join(head, tail)
        else:
            self.redfile = redfile
        
        self._read_skyfile()
            
        if self.grid == 'healpix':
            self._correct_healpix()
            
        elif self.grid in ['polar', 'cartesian']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid
            exit()
            
        return


    def _correct_healpix(self):
    
        # Healpix grid instance.
        hg = HealpixGrid(self.nx)
        
        with h5py.File(self.fLCfile, 'r') as f, h5py.File(self.redfile, 'r+') as g:
            
            ascc = f['header_table/ascc'].value
            ra = f['header_table/ra'].value
            dec = f['header_table/dec'].value
            
            skyidx = hg.find_gridpoint(ra, dec)
            
            # Add header to resulting file.
            g['header'].attrs['skyfile'] = self.skyfile
            try: del g['data2']
            except: pass
            for i in range(len(ascc)):
                
                lc = f['data/'+ascc[i]]
                rc = g['data/'+ascc[i]]
                
                lstidx = lc['lstidx'].astype('int')
            
                # Find stellar transmission curve and correct flux and eflux.
                skytrans0 = self.skytrans[skyidx[i], lstidx]
                scflux0 = rc['ipcflux0']-skytrans0
        
                # Add flags.
                flags = np.where(np.isnan(skytrans0), 1, 0)
                flags = flags + np.where(self.pointcount[skyidx[i], lstidx]<=5, 2, 0) 
                if self.starcount[skyidx[i]] <= 5:
                    flags += 4
                flags = flags + np.where(self.clouds[skyidx[i], lstidx]>.05, 8, 0)
        
                # Combine the reduced data in a record array.
                record = np.rec.fromarrays([skytrans0, scflux0, flags], names=['skytrans0', 'scflux0', 'flags'])
                
                # Write the result to file.
                g.create_dataset('data2/'+ascc[i], data=record)
    
        return
