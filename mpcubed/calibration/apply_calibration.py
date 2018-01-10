#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

from .. import io, misc, statistics

class ApplyDecorrelation(object):

    def __init__(self, sysfile):
        
        self.sysfile = sysfile
        
        f = io.SysFile(self.sysfile)
    
        ascc, vmag, mag, sigma, nobs = f.read_magnitudes()    
        self.magnitudes = dict()
        self.magnitudes['ascc'] = ascc
        self.magnitudes['nobs'] = nobs
        self.magnitudes['mag'] = mag
        self.magnitudes['sigma'] = sigma
    
        camgrid, trans, nobs = f.read_trans()
        self.camgrid = camgrid
        self.transmission = dict()
        self.transmission['nobs'] = nobs
        self.transmission['trans'] = trans
        
        ipxgrid, sinx, cosx, siny, cosy, nobs = f.read_intrapix()  
        self.ipxgrid = ipxgrid
        self.intrapix = dict()
        self.intrapix['nobs'] = nobs
        self.intrapix['sinx'] = sinx
        self.intrapix['cosx'] = cosx
        self.intrapix['siny'] = siny
        self.intrapix['cosy'] = cosy
        
        skygrid, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()        
        self.skygrid = skygrid
        self.clouds = dict()
        self.clouds['nobs'] = nobs
        self.clouds['clouds'] = clouds
        self.clouds['sigma'] = sigma
        self.lstmin = lstmin        
        
        return
        
    def get_magnitudes(self, ascc):

        mask = self.magnitudes['ascc'] == ascc
        mag = self.magnitudes['mag'][mask]
        nobs = self.magnitudes['nobs'][mask]
        sigma = self.magnitudes['sigma'][mask]

        return mag, nobs, sigma        
        
    def get_transmission(self, ra, dec, lst):
        
        ha = np.mod(lst*15 - ra, 360.)
        dec = np.repeat(dec, len(lst))

        k, n = self.camgrid.radec2idx(ha, dec)

        trans = self.transmission['trans'][k,n]
        nobs = self.transmission['nobs'][k,n]

        return trans, nobs
        
    def get_intrapix(self, ra, dec, lst, x, y):
        
        ha = np.mod(lst*15 - ra, 360.)
        dec = np.repeat(dec, len(lst))

        l, n = self.ipxgrid.radec2idx(ha, dec)
        
        ipx_x = self.intrapix['sinx'][l,n]*np.sin(2*np.pi*x) + self.intrapix['cosx'][l,n]*np.cos(2*np.pi*x)
        ipx_y = self.intrapix['siny'][l,n]*np.sin(2*np.pi*y) + self.intrapix['cosy'][l,n]*np.cos(2*np.pi*y)        
        nobs = self.intrapix['nobs'][l,n]        
        
        return ipx_x + ipx_y, nobs
        
    def get_clouds(self, ra, dec, lstseq):
      
        q = self.skygrid.radec2idx(ra, dec)
        t = lstseq - self.lstmin
        
        clouds = self.clouds['clouds'][q,t]
        nobs = self.clouds['nobs'][q,t]
        sigma = self.clouds['sigma'][q,t]
      
        return clouds, nobs, sigma
        
    def get_systematics(self, ascc, ra, dec, lstseq, lst, x, y):

        flag = np.zeros(len(lstseq), dtype='uint8')
        
        mag, nobs, sigma = self.get_magnitudes(ascc)        
        
        trans, nobs = self.get_transmission(ra, dec, lst)
        flag = np.where(nobs < 25, flag+2, flag)
        
        ipx, nobs = self.get_intrapix(ra, dec, lst, x, y)
        flag = np.where(nobs < 25, flag+4, flag)
        
        clouds, nobs, sigma = self.get_clouds(ra, dec, lstseq)
        flag = np.where(nobs < 25, flag+8, flag)
        flag = np.where(sigma > .05, flag+16, flag)
        
        systematics = trans + ipx + clouds
        flag = np.where(np.isnan(systematics), flag + 1, flag)
        
        return mag, trans, ipx, clouds, flag 

    def __call__(self, photfile, aper, redfile=None):
        
        # Check if the photometry file exists.
        if not os.path.isfile(photfile):
            raise IOError('Photometry file not found: {}'.format(photfile))
         
        # The reduced lightcurves file.
        if 'vmag' in self.sysfile:
            prefix = 'red{}_vmag_'.format(aper)
        else:
            prefix = 'red{}_'.format(aper)
        
        if redfile is None:
            head, tail = os.path.split(photfile)
            tail = prefix + tail.rsplit('_')[-1]
            redfile = os.path.join(head, tail)
        
        # Check if the reduced lightcurves file exists.
        if os.path.isfile(redfile):
            raise IOError('Reduced lightcurves file already exists: {}'.format(redfile))
        
        print 'Reading corrections from: {}'.format(self.sysfile) 
        print 'Applying corrections to aperture {} of file: {}'.format(aper, photfile)
        print 'Writing results to: {}'.format(redfile)
        
        # Read the raw photometry.
        f = io.PhotFile(photfile)
    
        settings = f.read_global()
    
        fields = ['ascc', 'ra', 'dec', 'vmag', 'bmag', 'spectype']
        stars = f.read_stars(fields)
        stars['nobs'] = np.zeros(len(stars['ascc']), dtype='uint32')  
        stars['lstseqmin'] = np.zeros(len(stars['ascc']), dtype='uint32')     
        stars['lstseqmax'] = np.zeros(len(stars['ascc']), dtype='uint32')         
           
        # Fields and datatypes for the binned lightcurves.
        lightcurves = dict()
        names = ['lstseq', 'nobs', 'lst', 'jdmid', 'x', 'y', 'sky', 'esky', 'mag%i'%aper, 'emag%i'%aper, 'trans%i'%aper, 'etrans%i'%aper, 'clouds%i'%aper, 'eclouds%i'%aper]
        formats = ['uint32', 'uint8', 'float64', 'float64', 'float32', 'float32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32']            
                
        for i in range(len(stars['ascc'])):                    
                
            # Read the lightcurve.
            lc = f.read_lightcurves(ascc=stars['ascc'][i])
    
            mask = (lc['flux{}'.format(aper)] > 0) & (lc['eflux{}'.format(aper)] > 0) & (lc['flag'] < 1) & (lc['sky'] > 0)
            lc = lc[mask]
    
            if (len(lc) == 0):
                stars['nobs'][i] = 0
                continue   
    
            # Compute the corrected lightcurve.
            x, y = lc['x'].astype('float64'), lc['y'].astype('float64')
            mag, emag = misc.flux2mag(lc['flux{}'.format(aper)], lc['eflux{}'.format(aper)])
            mag0, trans, ipx, clouds, cflag = self.get_systematics(stars['ascc'][i], stars['ra'][i], stars['dec'][i], lc['lstseq'], lc['lst'], x, y)
            trans = trans + mag0
            mag = mag - trans - ipx - clouds          
                
            # Mask values.
            mask = (cflag < 1)         
        
            lc = lc[mask]
            
            if (len(lc) == 0):
                stars['nobs'][i] = 0
                continue 
            
            x = x[mask]
            y = y[mask]
            mag = mag[mask]
            emag = emag[mask]
            trans = trans[mask]
            ipx = ipx[mask]
            clouds = clouds[mask]            
                
            # Compute the final binned lightcurve.
            lstseq, binidx, nobs = np.unique(lc['lstseq']//50, return_inverse=True, return_counts=True)        
            
            lc_bin = np.recarray(len(lstseq), names=names, formats=formats)        
            
            lc_bin['lstseq'] = lstseq
            lc_bin['nobs'] = nobs # Number of raw points used for each binned point.        
            
            lc_bin['jdmid'] = statistics.idxstats(binidx, lc['jdmid'], statistic='mean')
            lc_bin['lst'] = statistics.idxstats(binidx, lc['lst'], statistic='mean')
            
            lc_bin['x'] = statistics.idxstats(binidx, x, statistic='mean')
            lc_bin['y'] = statistics.idxstats(binidx, y, statistic='mean')        
            
            lc_bin['mag{}'.format(aper)] = statistics.idxstats(binidx, mag, statistic='mean')
            lc_bin['emag{}'.format(aper)] = statistics.idxstats(binidx, mag, statistic='std')/np.sqrt(nobs)
            lc_bin['sky'] = statistics.idxstats(binidx, lc['sky'], statistic='mean')
            lc_bin['esky'] = statistics.idxstats(binidx, lc['sky'], statistic='std')/np.sqrt(nobs)        
            
            lc_bin['trans{}'.format(aper)] = statistics.idxstats(binidx, trans, statistic='mean')
            lc_bin['etrans{}'.format(aper)] = statistics.idxstats(binidx, trans, statistic='std')/np.sqrt(nobs)
            lc_bin['clouds{}'.format(aper)] = statistics.idxstats(binidx, clouds, statistic='mean')
            lc_bin['eclouds{}'.format(aper)] = statistics.idxstats(binidx, clouds, statistic='std')/np.sqrt(nobs)            
            
            lightcurves[stars['ascc'][i]] = lc_bin
        
            stars['nobs'][i] = len(lstseq)
            stars['lstseqmin'][i] = lstseq[0]
            stars['lstseqmax'][i] = lstseq[-1]            
               
        io.write_reduced(redfile, settings, stars, lightcurves)
            
        return redfile
