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


    def calculate(self, fLCfile, redfile=None, skyfile=None):
        
        self.fLCfile = fLCfile
        
        # Filename of reduced data file.
        if redfile is None:
            head, tail = os.path.split(self.fLCfile)
            tail = 'red_'+tail.rsplit('_')[-1]
            self.redfile = os.path.join(head, tail)
        else:
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
        pointcount = np.bincount(idxuni)
            
        # Compute the sky transmission curve.
        skytrans, flux, niter, chisq, chisq_sky, chisq_flux, npoints, npars = sysrem(idxuni, staridx, cflux0, ecflux0, a2 = (1e7)*10**(self.vmag/-2.5), maxiter=250)
        
        array = np.full((len(skyidx), len(lstidx)), fill_value=np.nan)
        array[np.unravel_index(idx, (len(skyidx), len(lstidx)))] = skytrans
        
        array1 = np.full((len(skyidx), len(lstidx)), fill_value=np.nan)
        array1[np.unravel_index(idx, (len(skyidx), len(lstidx)))] = pointcount
        
        array2 = np.full((len(skyidx), len(lstidx)), fill_value=np.nan)
        array2[np.unravel_index(idx, (len(skyidx), len(lstidx)))] = chisq_sky
        
        for ind in range(len(skyidx)):
        
            with h5py.File(self.skyfile) as f:
                
                grp = f.create_group('data/%i'%skyidx[ind])
                grp.create_dataset('lstidx', data=lstidx)
                grp.create_dataset('pointcount', data=array1[ind])
                grp.create_dataset('skytrans', data=array[ind])
                grp.create_dataset('chisq_sky', data=array2[ind])
        
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
            
            self.starcount = np.full(npix, fill_value=0)
            self.skytrans = np.full((npix, 13500), fill_value=np.nan)
            self.pointcount = np.full((npix, 13500), fill_value=0)
            self.chisq_sky = np.full((npix, 13500), fill_value=np.nan)
            
            self.skyidx = f['header/skyidx'].value
            self.starcount[self.skyidx] = f['header/starcount'].value
            
            for idx in self.skyidx:
                
                try:
                    data = f['data/%i'%idx]
                except:
                    pass
                else:
                    lstidx = data['lstidx'].value
                    
                    self.skytrans[idx, lstidx] = data['skytrans'].value
                    self.pointcount[idx, lstidx] = data['pointcount'].value
                    self.chisq_sky[idx, lstidx] = data['chisq_sky'].value
            
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
        
        # Figure showing the transmission map.
        plt.imshow(self.skytrans, aspect='auto', cmap=viridis, vmin=0, vmax=1.5)
        cb = plt.colorbar()
        plt.xlabel('LST [idx]')
        plt.ylabel('Sky [idx]')
        plt.ylim(np.amin(self.skyidx)-.5,np.amax(self.skyidx)+.5) 
        cb.set_label('Sky')
        plt.show()
        
        for idx in self.skyidx:
            
            xlim, = np.where(np.isfinite(self.skytrans[idx]))
            
            plt.figure(figsize=(16,8))
            
            ax = plt.subplot(311)
            plt.title('skyidx = %i'%idx)
            plt.plot(self.skytrans[idx], '.')
            plt.ylim(0,1.5)
            plt.ylabel('Sky')
            
            plt.subplot(312, sharex=ax)
            plt.plot(self.pointcount[idx], '.')
            plt.ylabel('# points')
            
            plt.subplot(313, sharex=ax)
            plt.plot(self.chisq_sky[idx], '.')
            plt.xlim(np.amin(xlim)-.5, np.amax(xlim+.5))
            plt.ylabel(r'$\chi^2$')
            plt.xlabel('LST [idx]')
            
            plt.tight_layout()
            plt.show()
            plt.close()
        
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
            tail = 'red_'+tail.rsplit('_')[-1]
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
            
            ascc = f['table_header/ascc'].value
            
            # Add header to resulting file.
            g['header'].attrs['skyfile'] = self.skyfile
            
            for sid in ascc:
                
                si = f['header/'+sid]
                lc = f['data/'+sid]
                rc = g['data/'+sid]
                
                # Find skyidx for the lightcurve.
                skyidx = hg.find_gridpoint(si['ra'], si['dec'])
                lstidx = lc['lstidx'].astype('int')
            
                # Find stellar transmission curve and correct flux and eflux.
                skytrans0 = self.skytrans[skyidx, lstidx]
                scflux0 = rc['cflux0']/skytrans0
                escflux0 = rc['ecflux0']/skytrans0
        
                # Add flags.
                flags = np.where(np.isnan(skytrans0), 1, 0)
                flags = flags + np.where(self.pointcount[skyidx, lstidx]<=5, 2, 0) 
                if self.starcount[skyidx] <= 5:
                    flags += 4
        
                # Combine the reduced data in a record array.
                record = np.rec.fromarrays([skytrans0, scflux0, escflux0, flags], names=['skytrans0', 'scflux0', 'escflux0', 'flags'])
                
                # Write the result to file.
                g.create_dataset('data2/'+sid, data=record)
        
        return
        
        
#test = SkyTransmission()
#test.calculate('/data2/talens/Jul2015/fLC_20150715LPN.hdf5')
#test.calculate('/data2/talens/Jul2015/fLC_20150715LPE.hdf5')
#test.calculate('/data2/talens/Jul2015/fLC_20150715LPS.hdf5')
#test.calculate('/data2/talens/Jul2015/fLC_20150715LPW.hdf5')
#test.calculate('/data2/talens/Jul2015/fLC_20150715LPC.hdf5')

#test = SkyFile('/data2/talens/Jul2015/sky_20150716LPN.hdf5')
#test.visualize()
#test.correct()









