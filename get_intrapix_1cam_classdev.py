#!/usr/bin/env python
# -*- coding: utf-8 -*-

class IntraPixel():
    
    def __init__(self):
        
        
    def calculate(self, fLCfile)
    
        self.fLCfile = fLCfile
        
        # Read the stellar header information.
        with h5py.File(self.fLCfile, 'r') as f:
            
            hdr = f['table_header']
            self.ascc = hdr['ascc'].value
            self.vmag = hdr['vmag'].value
            self.ra = hdr['ra'].value
            self.dec = hdr['dec'].value
            self.nobs = hdr['nobs'].value.astype('int')
            
        # Calculate the intrapixel variations.
        pg = PolarGrid(270, 720)
        decidx, decuni = pg.find_decidx(self.dec, compact=True)
        starcount = index_statistics(decuni, decuni, statistic='count')
        
        # Create arrays.
        niter = np.zeros(len(decidx), dtype='int')
        chisq = np.zeros(len(decidx), dtype='float')
        npoints = np.zeros(len(decidx), dtype='int')
        npars = np.zeros(len(decidx), dtype='int')
    
        for ind in range(len(decidx)):
            
            # Select stars in the current sky bin.
            here = (decuni == ind)
            ascc = self.ascc[here]
            ra = self.ra[here]
            vmag = self.vmag[here]
            nobs = self.nobs[here]
            
            # Read data for these stars.
            lst, flux0, eflux0, sky, flags = self._read_data(ascc, nobs)
            
            # Create the haidx and staridx. 
            ha = np.mod(lst*15.-np.repeat(ra,nobs), 360.)
            haidx = pg.find_raidx(ha)
            staridx = np.repeat(np.arange(len(ascc)), nobs)
            
            # Remove bad datapoints.
            here = (flux0 > 0)&(eflux0 > 0)&(sky > 0)&(flags < 1)
            flux0 = flux0[here]
            eflux0 = eflux0[here]
            haidx = haidx[here]
            staridx = staridx[here]

            if len(flux0) == 0: continue
            
            # Make the haidx ascending from 0 and count the number of datapoints at each haidx.
            haidx, hauni = np.unique(haidx, return_inverse=True)
            pointcount = index_statistics(hauni, hauni, statistic='count')
            
            a, b, niter[ind], chisq[ind] = intrapix(hauni, y, flux0, eflux0)
            
            A = np.sqrt(a**2+b**2)
            phi = np.arctan(a/b)
            
            with h5py.File(self.camfile) as f:
                
                grp = f.create_group('data/%i'%decidx[ind])
                grp.create_dataset('haidx', data=haidx)
                grp.create_dataset('pointcount', data=pointcount)
                grp.create_dataset('amplitude', data=A)
                grp.create_dataset('phase', data=phi)
            
         with h5py.File(self.camfile) as f:
            
            grp = f.create_group('header')
            
            grp.create_dataset('decidx', data = decidx)
            grp.create_dataset('niter', data = niter)
            grp.create_dataset('chisq', data = chisq)
        
            
            
            
            
            
            
            
