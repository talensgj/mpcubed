#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from progressbar import ProgressBar

from fLCfile import fLCfile
from sysfile import SysFile
from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics
from usefull_functions_dev import flux2mag

class SysCorr():
    
    def __init__(self, LBfile, sysfile, outfile):
        
        self.LBfile = LBfile
        self.sysfile = sysfile
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
            flux, eflux, sky, x, y, jdmid, lst, lstseq, flags = f.read_data(['flux0', 'eflux0', 'sky', 'x', 'y', 'jdmid', 'lst', 'lstseq', 'flag'], [ascc], [nobs])
            lstseq = lstseq.astype('int') - self.lstmin
                
            ra = np.repeat(ra, nobs)
            dec = np.repeat(dec, nobs)
            ha = np.mod(lst*15. - ra, 360.)
            
            # Create indices.
            camtransidx = self.camgrid.find_gridpoint(ha, dec)
            intrapixidx = self.ipxgrid.find_gridpoint(ha, dec)
            
            # Convert flux to magnitudes.
            mag, emag = flux2mag(flux, eflux)

            # Get correction
            camtrans = self.z[camtransidx]
            intrapix = self.a[intrapixidx]*np.sin(2*np.pi*x) + self.b[intrapixidx]*np.cos(2*np.pi*x) + self.c[intrapixidx]*np.sin(2*np.pi*y) + self.d[intrapixidx]*np.cos(2*np.pi*y)
            skytrans = self.s[skyidx, lstseq]
            correction = camtrans + intrapix + skytrans

            # Get new flags.
            flags = np.where((flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1), 0, 1)
            flags = flags + np.where(np.isnan(correction), 2, 0)
            flags = flags + np.where((self.nobs_z[camtransidx] < 25) | (self.nobs_A[intrapixidx] < 25) | (self.nobs_s[skyidx, lstseq] < 25), 4, 0)
        
            binidx = (lstseq + self.lstmin) // 50
            here = (flags < 1)
            
            binidx = binidx[here]
            mag = mag[here]
            correction = correction[here]
            jdmid = jdmid[here]
            lst = lst[here]
            x = x[here]
            y = y[here]
            sky = sky[here]
            camtrans = camtrans[here]
            skytrans = skytrans[here]
            
            lstseq = np.unique(binidx)
            nobs = index_statistics(binidx, None, statistic = 'count')
            bmag = index_statistics(binidx, mag - correction, statistic = 'mean')
            emag = index_statistics(binidx, mag - correction, statistic = 'std')
            jdmid = index_statistics(binidx, jdmid, statistic = 'mean')
            lst = index_statistics(binidx, lst, statistic = 'mean')
            x = index_statistics(binidx, x, statistic = 'mean')
            y = index_statistics(binidx, y, statistic = 'mean')
            bsky = index_statistics(binidx, sky, statistic = 'mean')
            esky = index_statistics(binidx, sky, statistic = 'std')
            bcamtrans = index_statistics(binidx, camtrans, statistic = 'mean')
            ecamtrans = index_statistics(binidx, camtrans, statistic = 'std')
            bskytrans = index_statistics(binidx, skytrans, statistic = 'mean')
            eskytrans = index_statistics(binidx, skytrans, statistic = 'std')
        
            with h5py.File(self.outfile) as g:
                g.create_dataset('data/' + ascc + '/lstseq', data = lstseq)
                g.create_dataset('data/' + ascc + '/nobs', data = nobs)
                g.create_dataset('data/' + ascc + '/mag0', data = bmag)
                g.create_dataset('data/' + ascc + '/emag0', data = emag)
                g.create_dataset('data/' + ascc + '/jdmid', data = jdmid)
                g.create_dataset('data/' + ascc + '/lst', data = lst)
                g.create_dataset('data/' + ascc + '/x', data = x)
                g.create_dataset('data/' + ascc + '/y', data = y)
                g.create_dataset('data/' + ascc + '/sky', data = bsky)
                g.create_dataset('data/' + ascc + '/esky', data = esky)
                g.create_dataset('data/' + ascc + '/camtrans', data = bcamtrans)
                g.create_dataset('data/' + ascc + '/ecamtrans', data = ecamtrans)
                g.create_dataset('data/' + ascc + '/skytrans', data = bskytrans)
                g.create_dataset('data/' + ascc + '/eskytrans', data = eskytrans)
            
            pbar.update(i + 1)
           
        pbar.finish()
            
        return

    def run(self):
        
        # Read the correction terms.
        sys = SysFile(self.sysfile)
        z, nobs_z = sys.read_camtrans()
        a, b, c, d, nobs_A = sys.read_intrapix()
        s, nobs_s, lstmin = sys.read_skytrans()
    
        # Should come from SysFile???
        self.z = z.ravel()
        self.nobs_z = nobs_z.ravel()
        self.a = a.ravel()
        self.b = b.ravel()
        self.c = c.ravel()
        self.d = d.ravel()
        self.s = s
        self.nobs_s = nobs_s
        self.nobs_A = nobs_A.ravel()
        self.camgrid = PolarGrid(13500, 720)
        self.ipxgrid = PolarGrid(270, 720)
        self.lstmin = lstmin
        
        # Put in self
        f = fLCfile('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')
        ascc, ra, dec, vmag, nobs = f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs'])
        nobs = nobs.astype('int')
        
        self.ascc = ascc
        self.ra = ra
        self.dec = dec
        self.vmag = vmag
        self.nobs = nobs
        
        skygrid = HealpixGrid(8)
        self.skyidx = skygrid.find_gridpoint(ra, dec)
        
        self.correct()
        
        return

if __name__ == '__main__':
    
    obj = SysCorr('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', '/data2/talens/2015Q2/LPE/sys_201506ALPE.hdf5', '/data2/talens/2015Q2/LPE/cor_201506ALPE.hdf5')
    obj.run()
