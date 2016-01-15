#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from progressbar import ProgressBar

import misc
import IO
from coordinates import grids
from systematics import cdecor
from statistics import statistics

def create_baseline(date, camera, mode, filepath, outpath):
    
    if (mode == 0):
        part = 'A'
        dates = [date + '%.2i'%i + camera for i in range(1, 16)]
    elif (mode == 1):
        part = 'B'
        dates = [date + '%.2i'%i + camera for i in range(16, 32)]
    elif (mode == 2):
        part = ''
        dates = [date + '%.2i'%i + camera for i in range(1, 32)]
    else:
        print 'Unknown value for mode.'
        print 'exiting...'
        exit()
        
    outfile = os.path.join(outpath, 'fLC_%s%s%s.hdf5'%(date, part, camera))
    filelist = [os.path.join(filepath, '%s/fLC/fLC_%s.hdf5'%(date, date)) for date in dates]
    
    IO.make_baseline(filelist, outfile)
    
    return outfile

class CoarseDecorrelation():
    
    def __init__(self, LBfile, aperture, sysfile = None, **kwargs):
        """
            Performs a coarse decorrelation on all data in a given Long Baseline
            fLC file. Removes camera transmission, intrapixel variations and
            sky transmission.
        """
        
        # fLC file and aperture to work on.
        self.LBfile = LBfile
        self.aper = aperture
        
        if not os.path.isfile(self.LBfile):
            print 'File not found:', self.LBfile
            print 'exiting...'
            exit()
        else:
            print 'Calculating corrections for aperture %i of file:'%self.aper, self.LBfile
        
        # The systematics file.
        if sysfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'sys%i_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            sysfile = os.path.join(head, tail)
        
        self.sysfile = sysfile
        
        if os.path.isfile(self.sysfile):
            print 'Systematics file already exists:', self.sysfile
            print 'exiting...'
            exit()
        else:
            print 'Writing results to:', self.sysfile
        
        # Initialize with defualt parameters unless arguments were given.
        self.sigmas = kwargs.pop('sigmas', True)
        self.outer_maxiter = kwargs.pop('outer_maxiter', 5)
        
        self.camgrid = 'polar'
        self.camnx = kwargs.pop('camnx', 13500)
        self.camny = kwargs.pop('camny', 720)
        
        self.ipxgrid = 'polar'
        self.ipxnx = kwargs.pop('ipxnx', 270)
        self.ipxny = self.camny
        
        self.skygrid = 'healpix'
        self.skynx = kwargs.pop('skynx', 8)
    
        # Options to coarse_decor functions.
        self.inner_maxiter = kwargs.pop('inner_maxiter', 100)
        self.dtol = kwargs.pop('dtol', 1e-3)
        self.verbose = kwargs.pop('verbose', False)
    
        return
    
    def read_data(self, here):
        """
            Reads portions of the data from file, creates appropriate indices,
            removes flagged data and converts fluxes to magnitudes.
        """
        
        ascc = self.ascc[here]
        ra = self.ra[here]
        dec = self.dec[here]
        nobs = self.nobs[here]
        
        staridx = self.staridx[here]
        skyidx = self.skyidx[here]
        
        # Read data.
        fields = ['flux%i'%self.aper, 'eflux%i'%self.aper, 'sky', 'x', 'y', 'lst', 'lstseq', 'flag']
        flux, eflux, sky, x, y, lst, lstseq, flags = self.f.read_data(fields, ascc, nobs)
        lstseq = lstseq.astype('int') - self.lstmin
        
        ra = np.repeat(ra, nobs)
        dec = np.repeat(dec, nobs)
        ha = np.mod(lst*15. - ra, 360.)
        
        # Create indices.
        staridx = np.repeat(staridx, nobs)
        camtransidx, decidx = self.camgrid.radec2idx(ha, dec)
        intrapixidx, decidx = self.ipxgrid.radec2idx(ha, dec)
        skyidx = np.repeat(skyidx, nobs)
        
        # Remove bad data.
        here = (flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1)
        
        flux = flux[here]
        eflux = eflux[here]
        x = x[here]
        y = y[here]
        lstseq = lstseq[here]
        
        staridx = staridx[here]
        decidx = decidx[here]
        camtransidx = camtransidx[here]
        intrapixidx = intrapixidx[here]
        skyidx = skyidx[here]
        
        # Convert flux to magnitudes:
        mag, emag = misc.flux2mag(flux, eflux)
        
        return mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, skyidx, lstseq
    
    def spatial(self):
        """
            Solves for the camera transmission and intrapixel variations which
            are assumed to be time independent.
        """
        
        decidx, decuni = np.unique(self.decidx, return_inverse=True)
        
        nbins = len(decidx)
        niter = np.zeros(nbins)
        chisq = np.zeros(nbins)
        npoints = np.zeros(nbins)
        npars = np.zeros(nbins)
        
        pbar = ProgressBar(maxval = nbins).start()
        for idx in range(nbins):
            
            # Read data.
            here = (decuni == idx)
            mag, emag, x, y, staridx, _, camtransidx, intrapixidx, skyidx, lstseq = self.read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known temporal correction.
            if self.got_sky == True:
                mag = mag - self.s[skyidx, lstseq]
                
                if self.sigmas:
                    emag = np.sqrt(emag**2 + self.sigma1[staridx]**2 + self.sigma2[skyidx, lstseq]**2)
            
            # Create unique indices.
            staridx, ind1, m_nobs = np.unique(staridx, return_inverse=True, return_counts=True)
            camtransidx, ind2, z_nobs = np.unique(camtransidx, return_inverse=True, return_counts=True)
            intrapixidx, ind3, A_nobs = np.unique(intrapixidx, return_inverse=True, return_counts=True)
            
            # Calculate new spatial correction.
            m, z, A, niter[idx], chisq[idx], npoints[idx], npars[idx] = cdecor.cdecor_intrapix(ind1, ind2, ind3, mag, emag, x, y, maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            
            # Simple magnitude calibration.
            offset = np.nanmedian(m - self.vmag[staridx])
            
            self.m[staridx] = m - offset
            self.z[camtransidx, decidx[idx]] = z + offset
            self.A[intrapixidx, decidx[idx]] = A
            
            self.m_nobs[staridx] = m_nobs
            self.z_nobs[camtransidx, decidx[idx]] = z_nobs
            self.A_nobs[intrapixidx, decidx[idx]] = A_nobs
            
            pbar.update(idx + 1)
            
        pbar.finish()
            
        self.niter_spatial = niter
        self.chisq_spatial = chisq
        self.npoints_spatial = npoints
        self.npars_spatial = npars 
            
        return
        
    def temporal(self):
        """
            Solves for the sky transmission as a function of position and time.
        """
        
        skyidx, skyuni = np.unique(self.skyidx, return_inverse=True)
        
        nbins = len(skyidx)
        niter = np.zeros(nbins)
        chisq = np.zeros(nbins)
        npoints = np.zeros(nbins)
        npars = np.zeros(nbins)
        
        pbar = ProgressBar(maxval = nbins).start()
        for idx in range(nbins):
            
            # Read data.
            here = (skyuni == idx)
            mag, emag, x, y, staridx, decidx, camtransidx, intrapixidx, _, lstseq = self.read_data(here)
            
            if (len(mag) == 0): continue
            
            # Apply known spatial correction.
            mag = mag - self.z[camtransidx, decidx]
            mag = mag - np.sum(self.A[intrapixidx, decidx]*np.array([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T, axis=1)
            
            # Create unique indices.
            staridx, ind1, m_nobs = np.unique(staridx, return_inverse=True, return_counts=True)
            lstseq, ind2, s_nobs = np.unique(lstseq, return_inverse=True, return_counts=True)
            
            # Calculate new temporal correction.
            if self.sigmas:
                m, s, sigma1, sigma2, niter[idx], chisq[idx], npoints[idx], npars[idx] = cdecor.cdecor_sigmas(ind1, ind2, mag, emag, self.sigma1[staridx], self.sigma2[skyidx[idx], lstseq], maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            else:
                m, s, niter[idx], chisq[idx], npoints[idx], npars[idx] = cdecor.cdecor(ind1, ind2, mag, emag, maxiter = self.inner_maxiter, dtol = self.dtol, verbose = self.verbose)
            
            # Simple magnitude calibration.
            offset = np.nanmedian(m - self.vmag[staridx])
            
            self.m[staridx] = m - offset
            self.s[skyidx[idx], lstseq] = s + offset
            
            if self.sigmas:
                self.sigma1[staridx] = sigma1
                self.sigma2[skyidx[idx], lstseq] = sigma2
            
            self.m_nobs[staridx] = m_nobs
            self.s_nobs[skyidx[idx], lstseq] = s_nobs
            
            pbar.update(idx + 1)
           
        pbar.finish()
           
        self.niter_temporal = niter
        self.chisq_temporal = chisq
        self.npoints_temporal = npoints
        self.npars_temporal = npars 
           
        self.got_sky = True
            
        return
    
    def run(self):
        """
            Performs the coarse decorrelation on the given Long Baseline fLC file
            and writes the result to the given output location.
        """

        # Set up the IO and coordinate grids.
        self.f = IO.fLCfile(self.LBfile)
        self.camgrid = grids.PolarGrid(self.camnx, self.camny)
        self.ipxgrid = grids.PolarGrid(self.ipxnx, self.ipxny)
        self.skygrid = grids.HealpixGrid(self.skynx)
        
        # Read the required header data.
        self.ascc, self.ra, self.dec, self.nobs, self.vmag = self.f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
        self.nobs = self.nobs.astype('int')
        
        # Create indices.
        self.staridx = np.arange(len(self.ascc))
        _, self.decidx = self.camgrid.radec2idx(self.ra, self.dec)
        self.skyidx = self.skygrid.radec2idx(self.ra, self.dec)
        
        # Read global information.
        with h5py.File(self.LBfile, 'r') as f:
            
            station = f['global'].attrs['station']
            camera = f['global'].attrs['camera']
            
            alt0 = f['global/alt0'].value
            az0 = f['global/az0'].value
            th0 = f['global/th0'].value
            x0 = f['global/x0'].value
            y0 = f['global/y0'].value
            
            lstmin = f['global'].attrs['lstmin'].astype('int')
            lstmax = f['global'].attrs['lstmax'].astype('int')
        
        self.lstmin = lstmin
        lstlen = lstmax - lstmin + 1
        
        # Create arrays to hold the results.
        self.m = np.full(len(self.ascc), fill_value=np.nan)
        self.z = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        self.A = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2, 4), fill_value=np.nan)
        self.s = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        if self.sigmas:
            self.sigma1 = np.full(len(self.ascc), fill_value=0)
            self.sigma2 = np.full((self.skygrid.npix, lstlen), fill_value=0)
        
        self.m_nobs = np.full(len(self.ascc), fill_value=np.nan)
        self.z_nobs = np.full((self.camgrid.nx+2, self.camgrid.ny+2), fill_value=np.nan)
        self.A_nobs = np.full((self.ipxgrid.nx+2, self.ipxgrid.ny+2), fill_value=np.nan)
        self.s_nobs = np.full((self.skygrid.npix, lstlen), fill_value=np.nan)
        
        self.got_sky = False
        
        # Perform the coarse decorrelation.
        print 'Performing coarse decorrelation for:', self.LBfile
        for niter in range(self.outer_maxiter):
        
            print 'Iteration %i out of %i:'%(niter + 1, self.outer_maxiter)
            
            print '    Calculating camera systematics...'
            self.spatial()
            
            print '    Calculating atmospheric systematics...'
            self.temporal()
        
        # Write the results to file.
        with h5py.File(self.sysfile) as f:
    
            # Write the header.
            hdr = f.create_group('header')
            
            hdr.attrs['station'] = station
            hdr.attrs['camera'] = camera
            
            hdr.attrs['alt0'] = np.mean(alt0)
            hdr.attrs['az0'] = np.mean(az0)
            hdr.attrs['th0'] = np.mean(th0)
            hdr.attrs['x0'] = np.mean(x0)
            hdr.attrs['y0'] = np.mean(y0)
            
            hdr.attrs['outer_maxiter'] = self.outer_maxiter
            hdr.attrs['inner_maxiter'] = self.inner_maxiter
            hdr.attrs['sigmas'] = self.sigmas
            hdr.attrs['dtol'] = self.dtol
           
            hdr.create_dataset('spatial/niter', data = self.niter_spatial, dtype = 'uint32')
            hdr.create_dataset('spatial/chisq', data = self.chisq_spatial, dtype = 'float64')
            hdr.create_dataset('spatial/npoints', data = self.npoints_spatial, dtype = 'uint32')
            hdr.create_dataset('spatial/npars', data = self.npars_spatial, dtype = 'uint32')
            
            hdr.create_dataset('temporal/niter', data = self.niter_temporal, dtype = 'uint32')
            hdr.create_dataset('temporal/chisq', data = self.chisq_temporal, dtype = 'float64')
            hdr.create_dataset('temporal/npoints', data = self.npoints_temporal, dtype = 'uint32')
            hdr.create_dataset('temporal/npars', data = self.npars_temporal, dtype = 'uint32')
            
            # Write the data.
            grp = f.create_group('data')
            
            # Write the magnitudes.
            grp.create_dataset('magnitudes/ascc', data = self.ascc)
            grp.create_dataset('magnitudes/vmag', data = self.vmag, dtype = 'float32')
            grp.create_dataset('magnitudes/nobs', data = self.m_nobs, dtype = 'uint32')
            grp.create_dataset('magnitudes/mag', data = self.m, dtype = 'float32')
            if self.sigmas:
                grp.create_dataset('magnitudes/sigma', data = self.sigma1, dtype = 'float32')
            
            # Write the camera transmission.
            self.z = np.ravel(self.z)
            self.z_nobs = np.ravel(self.z_nobs)
            idx, = np.where(~np.isnan(self.z))
            grp.create_dataset('trans/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('trans/nobs', data = self.z_nobs[idx], dtype = 'uint32')
            grp.create_dataset('trans/trans', data = self.z[idx], dtype = 'float32')
            
            grp['trans'].attrs['grid'] = 'polar'
            grp['trans'].attrs['nx'] = self.camnx
            grp['trans'].attrs['ny'] = self.camny
            
            # Write the intrapixel variations.
            self.A = np.reshape(self.A, ((self.ipxgrid.nx+2)*(self.ipxgrid.ny+2), 4))
            self.A_nobs = np.ravel(self.A_nobs)
            idx, = np.where(~np.isnan(self.A[:,0]))
            grp.create_dataset('intrapix/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('intrapix/nobs', data = self.A_nobs[idx], dtype = 'uint32')
            grp.create_dataset('intrapix/sinx', data = self.A[idx,0], dtype = 'float32')
            grp.create_dataset('intrapix/cosx', data = self.A[idx,1], dtype = 'float32')
            grp.create_dataset('intrapix/siny', data = self.A[idx,2], dtype = 'float32')
            grp.create_dataset('intrapix/cosy', data = self.A[idx,3], dtype = 'float32')
            
            grp['intrapix'].attrs['grid'] = 'polar'
            grp['intrapix'].attrs['nx'] = self.ipxnx
            grp['intrapix'].attrs['ny'] = self.ipxny
            
            # Write the sky transmission.
            idx, lstseq = np.where(~np.isnan(self.s))
            grp.create_dataset('clouds/idx', data = idx, dtype = 'uint32')
            grp.create_dataset('clouds/lstseq', data = lstseq + self.lstmin, dtype = 'uint32')
            grp.create_dataset('clouds/nobs', data = self.s_nobs[idx, lstseq], dtype = 'uint32')
            grp.create_dataset('clouds/clouds', data = self.s[idx, lstseq], dtype = 'float32')
            if self.sigmas:
                grp.create_dataset('clouds/sigma', data = self.sigma2[idx, lstseq], dtype = 'float32')
            
            grp['clouds'].attrs['grid'] = 'healpix'
            grp['clouds'].attrs['nx'] = self.skynx
            grp['clouds'].attrs['lstmin'] = self.lstmin
            grp['clouds'].attrs['lstmax'] = lstmax
            grp['clouds'].attrs['lstlen'] = lstlen
        
        return


class SysCorr():
    
    def __init__(self, LBfile, aperture, sysfile = None, outfile = None):
        """
        Given an fLC file and a systematics file it computes the systematics
        corrected lightcurve, bins it and writes the result to a temporary file.
        """
        
        # fLC file and aperture to work on.
        self.LBfile = LBfile
        self.aper = aperture
        
        if not os.path.isfile(self.LBfile):
            print 'File not found:', self.LBfile
            print 'exiting...'
            exit()
        else:
            print 'Applying corrections to aperture %i of file:'%self.aper, self.LBfile
        
        # The systematics file.
        if sysfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'sys%i_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            sysfile = os.path.join(head, tail)
        
        self.sysfile = sysfile
        
        if not os.path.isfile(self.sysfile):
            print 'Systematics file not found:', self.sysfile
            print 'exiting...'
            exit()
        else:
            print 'Reading corrections from:', self.sysfile
        
        # The output file.
        if outfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'tmp%i_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            outfile = os.path.join(head, tail)
        
        self.outfile = outfile
        
        if os.path.isfile(self.outfile):
            print 'Output file already exists:', self.outfile
            print 'exiting...'
            exit()
        else:
            print 'Writing results to:', self.outfile
        
        return
        
    def correct(self):

        # Open the fLCfile.
        f = IO.fLCfile(self.LBfile)
        
        nstars = len(self.ascc)
        pbar = ProgressBar(maxval = nstars).start()
        for i in range(nstars):
            
            # Get the header information for this star.
            ascc = self.ascc[i]
            ra = self.ra[i]
            dec = self.dec[i]
            nobs = self.nobs[i]
            skyidx = self.skyidx[i]
            
            # Read data.
            fields = ['flux%i'%self.aper, 'eflux%i'%self.aper, 'sky', 'x', 'y', 'jdmid', 'lst', 'lstseq', 'flag']
            flux, eflux, sky, x, y, jdmid, lst, lstseq, flags = f.read_data(fields, [ascc], [nobs])
            lstseq = lstseq.astype('int') - self.lstmin
            
            # Convert flux to magnitudes.
            mag, emag = misc.flux2mag(flux, eflux)
            
            # Create indices.    
            ra = np.repeat(ra, nobs)
            dec = np.repeat(dec, nobs)
            ha = np.mod(lst*15. - ra, 360.)
            
            camtransidx = self.pgcam.find_gridpoint(ha, dec)
            intrapixidx = self.pgipx.find_gridpoint(ha, dec)
            
            # Get the correction terms.
            trans = self.trans[camtransidx]
            intrapix = self.a[intrapixidx]*np.sin(2*np.pi*x) + self.b[intrapixidx]*np.cos(2*np.pi*x) + self.c[intrapixidx]*np.sin(2*np.pi*y) + self.d[intrapixidx]*np.cos(2*np.pi*y)
            clouds = self.clouds[skyidx, lstseq]
            correction = trans + intrapix + clouds

            # Get new flags.
            flags = np.where((flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1), 0, 1)
            flags = flags + np.where(np.isnan(correction), 2, 0)
            flags = flags + np.where((self.nobs_trans[camtransidx] < 25) | (self.nobs_ipx[intrapixidx] < 25) | (self.nobs_clouds[skyidx, lstseq] < 25), 4, 0)
        
            # Determine the bins.
            binidx = (lstseq + self.lstmin) // 50
            
            # Remove bad data.
            here = (flags < 1)
            binidx = binidx[here]
            
            if (len(binidx) < 1):
                continue
            
            lst = lst[here]
            jdmid = jdmid[here]
            x = x[here]
            y = y[here]
            sky = sky[here]
            
            mag = mag[here]
            trans = trans[here]
            clouds = clouds[here]
            correction = correction[here]
            
            # Bin the data.
            lstseq = np.unique(binidx)
            nobs = statistics.idxstats(binidx, None, statistic = 'count')
            lst = statistics.idxstats(binidx, lst, statistic = 'mean')
            jdmid = statistics.idxstats(binidx, jdmid, statistic = 'mean')
            x = statistics.idxstats(binidx, x, statistic = 'mean')
            y = statistics.idxstats(binidx, y, statistic = 'mean')
            bsky = statistics.idxstats(binidx, sky, statistic = 'mean')
            esky = statistics.idxstats(binidx, sky, statistic = 'std')
            
            bmag = statistics.idxstats(binidx, mag - correction, statistic = 'mean')
            emag = statistics.idxstats(binidx, mag - correction, statistic = 'std')
            btrans = statistics.idxstats(binidx, trans, statistic = 'mean')
            etrans = statistics.idxstats(binidx, trans, statistic = 'std')
            bclouds = statistics.idxstats(binidx, clouds, statistic = 'mean')
            eclouds = statistics.idxstats(binidx, clouds, statistic = 'std')
            
            # Create a record array.
            arlist = [lstseq, nobs, lst, jdmid, x, y, bsky, esky, bmag, emag, btrans, etrans, bclouds, eclouds]
            names = ['lstseq', 'nobs', 'lst', 'jdmid', 'x', 'y', 'sky', 'esky', 'mag%i'%self.aper, 'emag%i'%self.aper, 'trans%i'%self.aper, 'etrans%i'%self.aper, 'clouds%i'%self.aper, 'eclouds%i'%self.aper]
            formats = ['uint32', 'uint8', 'float64', 'float64', 'float32', 'float32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']
            record = np.rec.fromarrays(arlist, names = names, formats = formats)
    
            # Write the lightcurve to file.
            with h5py.File(self.outfile) as g:
                g.create_dataset(ascc, data = record)
            
            pbar.update(i + 1)
           
        pbar.finish()
            
        return

    def run(self):
        
        # Read the required header data.
        f = IO.fLCfile(self.LBfile)
        self.ascc, self.ra, self.dec, self.vmag, self.nobs = f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
        self.nobs = self.nobs.astype('int')
        
        # Read the correction terms.
        sys = IO.SysFile(self.sysfile)
        self.pgcam, self.trans, self.nobs_trans = sys.read_trans(ravel = True)
        self.pgipx, self.a, self.b, self.c, self.d, self.nobs_ipx = sys.read_intrapix(ravel = True)
        self.hg, self.clouds, self.nobs_clouds, self.lstmin, lstmax = sys.read_clouds()
        
        # Create indices.
        self.skyidx = self.hg.find_gridpoint(self.ra, self.dec)
        
        # Apply the corrections.
        self.correct()
        
        return
        
def make_quarterfile(filelist, redfile):
    """ Merge the temporary files created by SysCorr."""
    ### NEEDS HEADER ###

    nfiles = len(filelist)
    
    # Create a list of unique stars in the files.
    ascc = np.array([])
    for i in range(nfiles):
        with h5py.File(filelist[i], 'r') as f:
            ascc = np.append(ascc, f.keys())
            
    ascc = np.unique(ascc)
    
    # Read the data.
    for sID in ascc:
        first = True
        for i in range(nfiles):
        
            with h5py.File(filelist[i], 'r') as f:
                
                # Try to read the star.
                try:
                    tmp = f[sID].value
                except:
                    continue
                
                # Add the data to the lightcurve.
                if first:
                    lc = tmp
                    first = False
                else:
                    lc = stack_arrays((lc, tmp), asrecarray=True)
        
        # Write the data to the redfile.
        with h5py.File(redfile) as f:
            for key in lc.dtype.names:
                f.create_dataset('data/' + sID + '/' + key, data = lc[key])

    return
