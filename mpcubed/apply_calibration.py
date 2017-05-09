#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import h5py
import numpy as np

from . import io, misc
from .statistics import statistics

class CorrectLC():
    
    def __init__(self, LBfile, aperture, sysfile = None):
        
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
            prefix = 'sys%i_vmag_'%self.aper
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
        self.f = io.PhotFile(self.LBfile)
        self.stars = self.f.read_stars(['ascc', 'ra', 'dec', 'nobs'])
        self.stars['nobs'] = self.stars['nobs'].astype('int')
        
        # Read the correction terms.
        sys = io.SysFile(self.sysfile)
        ascc, vmag, self.mag, sigma, nobs = sys.read_magnitudes()
        self.pgcam, self.trans, self.nobs_trans = sys.read_trans()
        self.pgipx, self.a, self.b, self.c, self.d, self.nobs_ipx = sys.read_intrapix()
        self.hg, self.clouds, self.sigma, self.nobs_clouds, self.lstmin, lstmax = sys.read_clouds()
        
        # Create indices.
        self.skyidx = self.hg.radec2idx(self.stars['ra'], self.stars['dec'])
        
        return
        
    def get_correction(self, ascc):

        # Get the header information for this star.
        i, = np.where(self.stars['ascc'] == ascc)
        ra = self.stars['ra'][i]
        dec = self.stars['dec'][i]
        nobs = self.stars['nobs'][i]
        skyidx = self.skyidx[i]
        
        # Read data.
        fields = ['x', 'y', 'lst', 'lstseq']
        lc = self.f.read_lightcurves(ascc=ascc, fields=fields)
        
        x = lc['x']
        y = lc['y']
        lst = lc['lst']
        lstseq = lc['lstseq'].astype('int') - self.lstmin
        
        # Create indices.    
        ra = np.repeat(ra, nobs)
        dec = np.repeat(dec, nobs)
        ha = np.mod(lst*15. - ra, 360.)
        
        camidx, decidx = self.pgcam.radec2idx(ha, dec)
        ipxidx, decidx = self.pgipx.radec2idx(ha, dec)
        
        # Get the correction terms.
        trans = self.trans[camidx, decidx] + self.mag[i]
        intrapix = self.a[ipxidx, decidx]*np.sin(2*np.pi*x) + self.b[ipxidx, decidx]*np.cos(2*np.pi*x) + self.c[ipxidx, decidx]*np.sin(2*np.pi*y) + self.d[ipxidx, decidx]*np.cos(2*np.pi*y)
        clouds = self.clouds[skyidx, lstseq]
        
        correction = trans + intrapix + clouds
        
        # Flag data where the correction may not be good.
        flags = np.where(np.isnan(correction), 1, 0)
        flags = flags + np.where((self.nobs_trans[camidx, decidx] < 25) | (self.nobs_ipx[ipxidx, decidx] < 25) | (self.nobs_clouds[skyidx, lstseq] < 25), 2, 0)
        if self.sigma is not None:
            flags = flags + np.where((self.sigma[skyidx, lstseq] > .05), 4, 0)
        
        return trans, intrapix, clouds, flags

    def binned_corrected_lightcurve(self, ascc):
        
        names = ['lstseq', 'nobs', 'lst', 'jdmid', 'x', 'y', 'sky', 'esky', 'mag%i'%self.aper, 'emag%i'%self.aper, 'trans%i'%self.aper, 'etrans%i'%self.aper, 'clouds%i'%self.aper, 'eclouds%i'%self.aper]
        formats = ['uint32', 'uint8', 'float64', 'float64', 'float32', 'float32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']        
        
        # Get the header information for this star.
        i, = np.where(self.stars['ascc'] == ascc)
        nobs = self.stars['nobs'][i]

        # Read data.
        fields = ['flux%i'%self.aper, 'eflux%i'%self.aper, 'sky', 'x', 'y', 'jdmid', 'lst', 'lstseq', 'flag']
        lc = self.f.read_lightcurves(ascc=ascc, fields=fields)
        
        flux = lc['flux%i'%self.aper]
        eflux = lc['eflux%i'%self.aper]
        sky = lc['sky']
        x = lc['x']
        y = lc['y']
        jdmid = lc['jdmid']
        lst = lc['lst']
        lstseq = lc['lstseq'].astype('int')
        flag = lc['flag']
        
        # Compute the corrected magnitudes.
        mag, emag = misc.flux2mag(flux, eflux)
        trans, intrapix, clouds, flags = self.get_correction(ascc)
        mag = mag - trans - intrapix - clouds
        
        # Remove flagged data.
        mask = (flux > 0) & (eflux > 0) & (flag < 1) & (flags < 1) & (sky > 0) 

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
        lstseq, binidx, nobs = np.unique(lstseq//50, return_inverse=True, return_counts=True)        
        
        lc_bin = np.recarray(len(lstseq), names=names, formats=formats)        
        
        lc_bin['lstseq'] = lstseq
        lc_bin['nobs'] = nobs # Number of raw points used for each binned point.        
        
        lc_bin['jdmid'] = statistics.idxstats(binidx, jdmid, statistic='mean')
        lc_bin['lst'] = statistics.idxstats(binidx, lst, statistic='mean')
        #lc_bin['exptime'] = idxstats(binidx, exptime, statistic='sum')           
        
        lc_bin['x'] = statistics.idxstats(binidx, x, statistic='mean')
        lc_bin['y'] = statistics.idxstats(binidx, y, statistic='mean')        
        
        lc_bin['mag{}'.format(self.aper)] = statistics.idxstats(binidx, mag, statistic='mean')
        lc_bin['emag{}'.format(self.aper)] = statistics.idxstats(binidx, mag, statistic='std')#/np.sqrt(nobs)
        lc_bin['sky'] = statistics.idxstats(binidx, sky, statistic='mean')
        lc_bin['esky'] = statistics.idxstats(binidx, sky, statistic='std')#/np.sqrt(nobs)        
        
        lc_bin['trans{}'.format(self.aper)] = statistics.idxstats(binidx, trans, statistic='mean')
        lc_bin['etrans{}'.format(self.aper)] = statistics.idxstats(binidx, trans, statistic='std')#/np.sqrt(nobs)
        lc_bin['clouds{}'.format(self.aper)] = statistics.idxstats(binidx, clouds, statistic='mean')
        lc_bin['eclouds{}'.format(self.aper)] = statistics.idxstats(binidx, clouds, statistic='std')#/np.sqrt(nobs)

        return lc_bin

    def make_redfile(self, outfile=None):
        
        # The output file.
        if outfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'red%i_vmag_'%self.aper
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
        with h5py.File(self.LBfile, 'r') as f:
            grp = f['global']
            filelist = grp['filelist'].value
            aversion = grp['aversion'].value
            rversion = grp['rversion'].value
            cversion = grp['cversion'].value
        
        with h5py.File(outfile) as f:
            grp = f.create_group('global')
            grp.attrs['station'] = data['station']
            grp.attrs['camera'] = data['camera']
            grp.attrs['exptime'] = 6.4 # Hardcoded ...
            grp.attrs['naper'] = data['naper']
            grp.attrs['aper0'] = data['aper0']
            grp.attrs['aper1'] = data['aper1']
            grp.attrs['skyrad0'] = data['skyrad0']
            grp.attrs['skyrad1'] = data['skyrad1']
            
            grp.create_dataset('filelist', data=filelist)
            grp.create_dataset('aversion', data=aversion)
            grp.create_dataset('rversion', data=rversion)
            grp.create_dataset('cversion', data=cversion)
            grp.attrs['pversion'] = '1.0.0' # Hardcoded ...
        
        # Write the header_table.
        fields = ['ascc', 'ra', 'dec', 'vmag', 'bmag', 'spectype']
        stars = self.f.read_stars(fields)
        with h5py.File(outfile) as f:
            grp = f.create_group('header_table')
            for key in stars.keys():
                grp.create_dataset(key, data=stars[key])
                
        # Write the corrected lightcurves.
        nstars = len(self.stars['ascc'])
        nobs = np.zeros(nstars)
        lstseqmin = np.zeros(nstars)
        lstseqmax = np.zeros(nstars)
        for i in range(nstars):
            
            # Get the binned corrected lightcurve.
            lc = self.binned_corrected_lightcurve(self.stars['ascc'][i])
            
            if (len(lc['jdmid']) > 0):
                nobs[i] = len(lc['jdmid'])
                lstseqmin[i] = np.amin(lc['lstseq'])
                lstseqmax[i] = np.amax(lc['lstseq'])
            
            # Write the lightcurve.
            with h5py.File(outfile) as f:
                f.create_dataset('data/' + self.stars['ascc'][i], data = lc)
            
        with h5py.File(outfile) as f:
            grp = f['header_table']
            grp.create_dataset('nobs', data = nobs)
            grp.create_dataset('lstseqmin', data = lstseqmin)
            grp.create_dataset('lstseqmax', data = lstseqmax)
            
        return outfile
