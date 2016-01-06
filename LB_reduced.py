#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from progressbar import ProgressBar

from fLCfile import fLCfile
from sysfile import SysFile

from core.index_functions import index_statistics
from usefull_functions_dev import flux2mag

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
        f = fLCfile(self.LBfile)
        
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
            mag, emag = flux2mag(flux, eflux)
            
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
            nobs = index_statistics(binidx, None, statistic = 'count')
            lst = index_statistics(binidx, lst, statistic = 'mean')
            jdmid = index_statistics(binidx, jdmid, statistic = 'mean')
            x = index_statistics(binidx, x, statistic = 'mean')
            y = index_statistics(binidx, y, statistic = 'mean')
            bsky = index_statistics(binidx, sky, statistic = 'mean')
            esky = index_statistics(binidx, sky, statistic = 'std')
            
            bmag = index_statistics(binidx, mag - correction, statistic = 'mean')
            emag = index_statistics(binidx, mag - correction, statistic = 'std')
            btrans = index_statistics(binidx, trans, statistic = 'mean')
            etrans = index_statistics(binidx, trans, statistic = 'std')
            bclouds = index_statistics(binidx, clouds, statistic = 'mean')
            eclouds = index_statistics(binidx, clouds, statistic = 'std')
            
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
        f = fLCfile(self.LBfile)
        self.ascc, self.ra, self.dec, self.vmag, self.nobs = f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
        self.nobs = self.nobs.astype('int')
        
        # Read the correction terms.
        sys = SysFile(self.sysfile)
        self.pgcam, self.trans, self.nobs_trans = sys.read_trans(ravel = True)
        self.pgipx, self.a, self.b, self.c, self.d, self.nobs_ipx = sys.read_intrapix(ravel = True)
        self.hg, self.clouds, self.nobs_clouds, self.lstmin, lstmax = sys.read_clouds()
        
        # Create indices.
        self.skyidx = self.hg.find_gridpoint(self.ra, self.dec)
        
        # Apply the corrections.
        self.correct()
        
        return

if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Apply the systematics corrections to a given fLC file.')
    parser.add_argument('LBfile', help = 'The input fLC file.')
    parser.add_argument('aperture', type = int, help = 'The aperture to be reduced.')
    parser.add_argument('-s', '--sysfile', help = 'The systematics file to use.', default = None)
    parser.add_argument('-o', '--outfile', help = 'The output file.', default = None)
    args = parser.parse_args()
    
    obj = SysCorr(args.LBfile, args.aperture, args.sysfile, args.outfile)
    obj.run()
    
    
