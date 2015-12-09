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
from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics
from usefull_functions_dev import flux2mag

class SysCorr():
    
    def __init__(self, LBfile, aperture, sysfile = None, outfile = None):
        
        self.LBfile = LBfile
        self.aper = aperture
        
        if sysfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'sys%i_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            sysfile = os.path.join(head, tail)
        
        self.sysfile = sysfile
        
        if outfile is None:
            head, tail = os.path.split(self.LBfile)
            prefix = 'tmp%i_'%self.aper
            tail = prefix + tail.rsplit('_')[-1]
            outfile = os.path.join(head, tail)
        
        self.outfile = outfile
        
        return
        
    def correct(self):

        nstars = len(self.ascc)
        
        f = fLCfile(self.LBfile)
        
        pbar = ProgressBar(maxval = nstars).start()
        for i in range(nstars):
            
            ascc = self.ascc[i]
            ra = self.ra[i]
            dec = self.dec[i]
            nobs = self.nobs[i]
            skyidx = self.skyidx[i]
            
            # Read data.
            flux, eflux, sky, x, y, jdmid, lst, lstseq, flags = f.read_data(['flux%i'%self.aper, 'eflux%i'%self.aper, 'sky', 'x', 'y', 'jdmid', 'lst', 'lstseq', 'flag'], [ascc], [nobs])
            lstseq = lstseq.astype('int') - self.lstmin
                
            ra = np.repeat(ra, nobs)
            dec = np.repeat(dec, nobs)
            ha = np.mod(lst*15. - ra, 360.)
            
            # Create indices.
            camtransidx = self.pgcam.find_gridpoint(ha, dec)
            intrapixidx = self.pgipx.find_gridpoint(ha, dec)
            
            # Convert flux to magnitudes.
            mag, emag = flux2mag(flux, eflux)

            # Get correction
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
        
        # Read the correction terms.
        sys = SysFile(self.sysfile)
        pgcam, trans, nobs_trans = sys.read_trans(ravel = True)
        pgipx, a, b, c, d, nobs_ipx = sys.read_intrapix(ravel = True)
        hg, clouds, nobs_clouds, lstmin, lstmax = sys.read_clouds()
    
        self.pgcam = pgcam
        self.trans = trans
        self.nobs_trans = nobs_trans
        
        self.pgipx = pgipx
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.nobs_ipx = nobs_ipx
        
        self.hg = hg
        self.clouds = clouds
        self.nobs_clouds = nobs_clouds
        self.lstmin = lstmin
        
        # Put in self
        f = fLCfile(self.LBfile)
        ascc, ra, dec, vmag, nobs = f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
        nobs = nobs.astype('int')
        
        self.ascc = ascc
        self.ra = ra
        self.dec = dec
        self.vmag = vmag
        self.nobs = nobs
        self.skyidx = self.hg.find_gridpoint(ra, dec)
        
        self.correct()
        
        return

def merge(filelist):
    
    nfiles = len(filelist)
    
    ascc = np.array([])
    for i in range(nfiles):
        
        with h5py.File(filelist[i], 'r') as f:
            ascc = np.append(ascc, f.keys())
            
    ascc = np.unique(ascc)
    
    for sID in ascc:
        print sID
        first = True
        for i in range(nfiles):
        
            with h5py.File(filelist[i], 'r') as f:
                
                try:
                    tmp = f[sID].value
                except:
                    pass
                else:
                    if first:
                        lc = tmp
                        first = False
                    else:
                        lc = stack_arrays((lc, tmp), asrecarray=True)
        
        with h5py.File('test.hdf5') as f:
            for key in lc.dtype.names:
                f.create_dataset('data/' + sID + '/' + key, data = lc[key])
        
        
        
    return

if __name__ == '__main__':
    
    #obj = SysCorr('/data2/talens/2015Q2/LPE/fLC_201506BLPE.hdf5', 0)
    #obj.run()
    
    filelist = glob.glob('/data2/talens/2015Q2/LPE/tmp0_*')
    filelist = np.sort(filelist)
    merge(filelist)
