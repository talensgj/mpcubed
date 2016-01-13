#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import h5py
import numpy as np

import IO
import misc
from core import statistics

class SysCorr():
    
    def __init__(self, LBfile, aperture, sysfile = None, outfile = None):
        
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
        
        # Read the required header data.
        self.f = IO.fLCfile(self.LBfile)
        self.ascc, self.ra, self.dec, self.vmag, self.nobs = self.f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
        self.nobs = self.nobs.astype('int')
        
        # Read the correction terms.
        sys = IO.SysFile(self.sysfile)
        self.pgcam, self.trans, self.nobs_trans = sys.read_trans()
        self.pgipx, self.a, self.b, self.c, self.d, self.nobs_ipx = sys.read_intrapix()
        self.hg, self.clouds, self.sigma, self.nobs_clouds, self.lstmin, lstmax = sys.read_clouds()
        
        # Create indices.
        self.skyidx = self.hg.radec2idx(self.ra, self.dec)
        
        return
        
    def get_correction(self, ascc):

        # Get the header information for this star.
        i, = np.where(self.ascc == ascc)
        ra = self.ra[i]
        dec = self.dec[i]
        nobs = self.nobs[i]
        skyidx = self.skyidx[i]
        
        # Read data.
        fields = ['x', 'y', 'lst', 'lstseq']
        x, y, lst, lstseq = self.f.read_data(fields, [ascc], [nobs])
        lstseq = lstseq.astype('int') - self.lstmin
        
        # Create indices.    
        ra = np.repeat(ra, nobs)
        dec = np.repeat(dec, nobs)
        ha = np.mod(lst*15. - ra, 360.)
        
        camidx, decidx = self.pgcam.radec2idx(ha, dec)
        ipxidx, decidx = self.pgipx.radec2idx(ha, dec)
        
        # Get the correction terms.
        trans = self.trans[camidx, decidx]
        intrapix = self.a[ipxidx, decidx]*np.sin(2*np.pi*x) + self.b[ipxidx, decidx]*np.cos(2*np.pi*x) + self.c[ipxidx, decidx]*np.sin(2*np.pi*y) + self.d[ipxidx, decidx]*np.cos(2*np.pi*y)
        clouds = self.clouds[skyidx, lstseq]
        
        correction = trans + intrapix + clouds
        
        # Flag data where the correction may not be good.
        flags = np.where(np.isnan(correction), 1, 0)
        flags = flags + np.where((self.nobs_trans[camidx, decidx] < 25) | (self.nobs_ipx[ipxidx, decidx] < 25) | (self.nobs_clouds[skyidx, lstseq] < 25), 2, 0)
        flags = flags + np.where((self.sigma[skyidx, lstseq] > .05), 4, 0)
        
        return trans, intrapix, clouds, flags

    def binned_corrected_lightcurve(self, ascc):
        
        # Get the header information for this star.
        i, = np.where(self.ascc == ascc)
        ra = self.ra[i]
        dec = self.dec[i]
        nobs = self.nobs[i]
        skyidx = self.skyidx[i]
        
        # Read data.
        fields = ['flux%i'%self.aper, 'eflux%i'%self.aper, 'sky', 'x', 'y', 'jdmid', 'lst', 'lstseq', 'flag']
        flux, eflux, sky, x, y, jdmid, lst, lstseq, flag = self.f.read_data(fields, [ascc], [nobs])
        lstseq = lstseq.astype('int')
        
        # Compute the corrected magnitudes.
        mag, emag = misc.flux2mag(flux, eflux)
        trans, intrapix, clouds, flags = self.get_correction(ascc)
        mag = mag - trans - intrapix - clouds
        
        # Remove flagged data.
        mask = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1) & (flags < 1)

        mag = mag[mask]
        emag = emag[mask]
        sky = sky[mask]
        x = x[mask]
        y = y[mask]
        jdmid = jdmid[mask]
        lst = lst[mask]
        lstseq = lstseq[mask]
        
        trans = trans[mask]
        clouds = clouds[mask]
        
        # Compute the final binned lightcurve.
        binidx = (lstseq // 50)
        lstseq = np.unique(binidx)
        nobs = statistics.idxstats(binidx, None, statistic = 'count')
        
        bmag = statistics.idxstats(binidx, mag, statistic = 'mean')
        emag = statistics.idxstats(binidx, mag, statistic = 'std')
        bsky = statistics.idxstats(binidx, sky, statistic = 'mean')
        esky = statistics.idxstats(binidx, sky, statistic = 'std')
        x = statistics.idxstats(binidx, x, statistic = 'mean')
        y = statistics.idxstats(binidx, y, statistic = 'mean')
        jdmid = statistics.idxstats(binidx, jdmid, statistic = 'mean')
        lst = statistics.idxstats(binidx, lst, statistic = 'mean')
        
        btrans = statistics.idxstats(binidx, trans, statistic = 'mean')
        etrans = statistics.idxstats(binidx, trans, statistic = 'std')
        bclouds = statistics.idxstats(binidx, clouds, statistic = 'mean')
        eclouds = statistics.idxstats(binidx, clouds, statistic = 'std')
        
        # Create a record array.
        arlist = [lstseq, nobs, lst, jdmid, x, y, bsky, esky, bmag, emag, btrans, etrans, bclouds, eclouds]
        names = ['lstseq', 'nobs', 'lst', 'jdmid', 'x', 'y', 'sky', 'esky', 'mag%i'%self.aper, 'emag%i'%self.aper, 'trans%i'%self.aper, 'etrans%i'%self.aper, 'clouds%i'%self.aper, 'eclouds%i'%self.aper]
        formats = ['uint32', 'uint8', 'float64', 'float64', 'float32', 'float32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']
        record = np.rec.fromarrays(arlist, names = names, formats = formats)

        return record

    def write_corrected_lightcurves(self, outfile=None):
        
        # The output file.
        if outfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'red%i_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            outfile = os.path.join(head, tail)
        
        if os.path.isfile(outfile):
            print 'Output file already exists:', outfile
            print 'exiting...'
            exit()
        else:
            print 'Writing results to:', outfile
        
        # Write the global.
        data = self.f.read_global()
        data = dict(data)
        with h5py.File(outfile) as f:
            grp = f.create_group('global')
            grp.attrs['station'] = data['station']
            grp.attrs['camera'] = data['camera']
        
        # Write the header_table.
        fields = ['ascc', 'ra', 'dec', 'vmag', 'bmag', 'spectype']
        data = self.f.read_header(fields)
        with h5py.File(outfile) as f:
            grp = f.create_group('header_table')
            grp.create_dataset('ascc', data = data[0])
            grp.create_dataset('ra', data = data[1])
            grp.create_dataset('dec', data = data[2])
            grp.create_dataset('vmag', data = data[3])
            grp.create_dataset('bmag', data = data[4])
            grp.create_dataset('spectype', data = data[5])
        
        # Write the corrected lightcurves.
        with h5py.File(outfile) as f:
            grp = f.create_group('data')
            
            nstars = len(self.ascc)
            for i in range(nstars):
                data = self.binned_corrected_lightcurve(self.ascc[i])
                grp.create_dataset(self.ascc[i], data = data)
            
        return outfile