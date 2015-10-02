#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from index_functions import index_statistics
from coordinate_grids import HealpixGrid
from coarse_decor import coarse_decorrelation

class SkyTransmission():
    
    def __init__(self, grid='healpix', nx=8, ny=0, margin=0):
        
        self.grid = grid
        self.nx = nx
        self.ny = ny
        self.margin = margin

        return 

    def _read_data(self, ascc, nobs):
        nstars = len(ascc)
        ndata = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        lstidx = np.zeros(ndata)
        sky = np.zeros(ndata)
        flags1 = np.zeros(ndata)
        
        cmag0 = np.zeros(ndata)
        ecmag0 = np.zeros(ndata)
        flags2 = np.zeros(ndata) 
        
        with h5py.File(self.fLCfile, 'r') as f, h5py.File(self.redfile, 'r') as g:
        
            lc = f['data']
            rc = g['data']
        
            for i in range(nstars):
                
                lstidx[select[i]:select[i+1]] = lc[ascc[i]]['lstidx']
                sky[select[i]:select[i+1]] = lc[ascc[i]]['sky']
                flags1[select[i]:select[i+1]] = lc[ascc[i]]['flag']
                
                cmag0[select[i]:select[i+1]] = rc[ascc[i]]['ipcmag0']
                ecmag0[select[i]:select[i+1]] = rc[ascc[i]]['smag0']
                flags2[select[i]:select[i+1]] = rc[ascc[i]]['flags']
        
        lstidx = lstidx.astype('int')

        return lstidx, sky, flags1, cmag0, ecmag0, flags2

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
            
            hdr = f['header_table']
            self.ascc = hdr['ascc'].value
            self.vmag = hdr['vmag'].value
            self.ra = hdr['ra'].value
            self.dec = hdr['dec'].value
            self.nobs = hdr['nobs'].value.astype('int')
        
        # Calculate the transmission.
        if self.grid == 'healpix':
            self._calculate_healpix()
             
        elif self.grid in ['cartesian', 'polar']:
            print 'Grid %s not implemented yet.'%self.grid
            exit()
            
        else:
            print 'Unknown grid: %s'%self.grid

        return

    def _calculate_healpix(self):
        
        # Healpix grid instance.
        hg = HealpixGrid(self.nx)
        
        # Assign stars to sky bins and count the number of stars in each bin.
        skyidx, skyuni = hg.find_gridpoint(self.ra, self.dec, compact=True)
        starcount = np.bincount(skyuni)
        
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
            lstidx, sky, flags1, cmag0, ecmag0, flags2 = self._read_data(ascc, nobs)
            
            # Create the staridx
            staridx = np.repeat(np.arange(len(ascc)), nobs)
            
            # Remove bad datapoints.
            here = np.isfinite(cmag0)&(ecmag0 > 0)&(sky > 0)&(flags1 < 1)&(flags2 < 1)
            cmag0 = cmag0[here]
            ecmag0 = ecmag0[here]
            lstidx = lstidx[here]
            staridx = staridx[here]
            
            if len(cmag0) == 0: continue
            
            # Make the lstidx ascending from 0 and count the number of datapoints at each lstidx.
            lstidx, lstuni = np.unique(lstidx, return_inverse=True)
            pointcount = np.bincount(lstuni)
            
            # Compute the sky transmission curve.
            magnitude, skytrans, sigma1, sigma2, niter[ind], chisq[ind], npoints[ind], npars[ind] = coarse_decorrelation(staridx, lstuni, cmag0, ecmag0)
        
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
            
            import healpy
            npix = healpy.nside2npix(self.nx)
            
            self.starcount = np.full(npix, fill_value=0)
            self.skytrans = np.full((13500, npix), fill_value=np.nan)
            self.pointcount = np.full((13500, npix), fill_value=0)
            
            self.skyidx = f['header/skyidx'].value
            self.starcount[self.skyidx] = f['header/starcount'].value
            
            for idx in self.skyidx:
                
                try:
                    data = f['data/%i'%idx]
                except:
                    pass
                else:
                    lstidx = data['lstidx'].value
                    self.skytrans[lstidx, idx] = data['skytrans'].value
                    self.pointcount[lstidx, idx] = data['pointcount'].value

        return
        
    def visualize(self, wrap=False):
        
        self.wrap = wrap
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
        import healpy
        
        rcParams['xtick.labelsize'] = 'large'
        rcParams['ytick.labelsize'] = 'large'
        rcParams['axes.labelsize'] = 'x-large'
        rcParams['image.interpolation'] = 'none'
        rcParams['image.origin'] = 'lower'
        
        array = self.skytrans
        #array = array/np.nanmedian(array, axis=0, keepdims=True)
        if self.wrap: array = np.roll(array, 13500/2, axis=0)
        
        xlim, ylim = np.where(np.isfinite(array))
        
        ax = plt.subplot(111)
        plt.imshow(array.T, aspect='auto', cmap=viridis, vmin=-.2, vmax=.2)
        cb = plt.colorbar()
        plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
        plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
        plt.xlabel('LST [idx]')
        plt.ylabel('sky [idx]')
        cb.set_label('Sky')
        plt.show()
        
        for idx in np.unique(xlim):
        
            healpy.mollview(array[idx], cmap=viridis, min=-.2, max=.2, unit='Sky', title='')
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
                skytrans0 = self.skytrans[lstidx, skyidx[i]]
                scmag0 = rc['ipcmag0']-skytrans0
        
                # Add flags.
                flags = np.where(np.isnan(skytrans0), 1, 0)
                flags = flags + np.where(self.pointcount[lstidx, skyidx[i]]<=5, 2, 0) 
        
                # Combine the reduced data in a record array.
                arlist = [skytrans0, scmag0, flags]
                names = ['skytrans0', 'scmag0', 'flags']
                record = np.rec.fromarrays(arlist, names=names)
                
                # Write the result to file.
                g.create_dataset('data2/'+ascc[i], data=record)
    
        return
