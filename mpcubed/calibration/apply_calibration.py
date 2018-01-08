#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import h5py
import numpy as np

from .. import io, misc, statistics

class GetSystematics(object):

    def __init__(self, filename):
        
        f = io.SysFile(filename)
    
        ascc, vmag, mag, sigma, nobs = f.read_magnitudes()    
        self.magnitudes = dict()
        self.magnitudes['ascc'] = ascc
        self.magnitudes['nobs'] = nobs
        self.magnitudes['mag'] = mag
        self.magnitudes['sigma'] = sigma
    
        pg, trans, nobs = f.read_trans()
        self.pgcam = pg
        self.transmission = dict()
        self.transmission['nobs'] = nobs
        self.transmission['trans'] = trans
        
        pg, sinx, cosx, siny, cosy, nobs = f.read_intrapix()  
        self.pgipx = pg
        self.intrapix = dict()
        self.intrapix['nobs'] = nobs
        self.intrapix['sinx'] = sinx
        self.intrapix['cosx'] = cosx
        self.intrapix['siny'] = siny
        self.intrapix['cosy'] = cosy
        
        hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()        
        self.hg = hg
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

        idx1, idx2 = self.pgcam.radec2idx(ha, dec)

        trans = self.transmission['trans'][idx1,idx2]
        nobs = self.transmission['nobs'][idx1,idx2]

        return trans, nobs
        
    def get_intrapix(self, ra, dec, lst, x, y):
        
        ha = np.mod(lst*15 - ra, 360.)
        dec = np.repeat(dec, len(lst))

        idx1, idx2 = self.pgipx.radec2idx(ha, dec)
        
        ipx_x = self.intrapix['sinx'][idx1,idx2]*np.sin(2*np.pi*x) + self.intrapix['cosx'][idx1,idx2]*np.cos(2*np.pi*x)
        ipx_y = self.intrapix['siny'][idx1,idx2]*np.sin(2*np.pi*y) + self.intrapix['cosy'][idx1,idx2]*np.cos(2*np.pi*y)        
        nobs = self.intrapix['nobs'][idx1,idx2]        
        
        return ipx_x + ipx_y, nobs
        
    def get_clouds(self, ra, dec, lstseq):
      
        idx1 = self.hg.radec2idx(ra, dec)
        idx2 = lstseq - self.lstmin
        
        clouds = self.clouds['clouds'][idx1,idx2]
        nobs = self.clouds['nobs'][idx1,idx2]
        sigma = self.clouds['sigma'][idx1,idx2]
      
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

def apply_calibration(LBfile, aper, sysfile=None, outfile=None):
    
    # fLC file and aperture to work on.
    if not os.path.isfile(LBfile):
        print 'File not found:', LBfile
        print 'exiting...'
        exit()
    else:
        print 'Applying corrections to aperture %i of file:'%aper, LBfile
    
    # The systematics file.
    if sysfile is None:
        head, tail = os.path.split(LBfile)
        prefix = 'sys%i_vmag_'%aper
        tail = prefix + tail.rsplit('_')[-1]
        sysfile = os.path.join(head, tail)
    
    if not os.path.isfile(sysfile):
        print 'Systematics file not found:', sysfile
        print 'exiting...'
        exit()
    else:
        print 'Reading corrections from:', sysfile    
    
    # The output file.
    if outfile is None:
        head, tail = os.path.split(LBfile)
        prefix = 'red%i_vmag_'%aper
        tail = prefix + tail.rsplit('_')[-1]
        outfile = os.path.join(head, tail)
    
    if os.path.isfile(outfile):
        print 'Output file already exists:', outfile
        print 'exiting...'
        exit()
    else:
        print 'Writing results to:', outfile
    
    # Read the stars 
    f = io.PhotFile(LBfile) 
    
    # Write the global.
    data = f.read_global()
    data = dict(data)
    
    with h5py.File(LBfile, 'r') as f:
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
    
    # Read the stars field from the photometry.
    f = io.PhotFile(LBfile)

    fields = ['ascc', 'ra', 'dec', 'vmag', 'bmag', 'spectype']
    stars = f.read_stars(fields)
    stars['nobs'] = np.zeros(len(stars['ascc']), dtype='uint32')  
    stars['lstseqmin'] = np.zeros(len(stars['ascc']), dtype='uint32')     
    stars['lstseqmax'] = np.zeros(len(stars['ascc']), dtype='uint32')         
       
    # Open the systematics file.
    sys = GetSystematics(sysfile)       
       
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
        mag0, trans, ipx, clouds, cflag = sys.get_systematics(stars['ascc'][i], stars['ra'][i], stars['dec'][i], lc['lstseq'], lc['lst'], x, y)
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
            
    with h5py.File(outfile) as f:
        
        idx, = np.where(stars['nobs'] > 0)            
        
        grp = f.create_group('header_table')
        for key in stars.keys():
            grp.create_dataset(key, data=stars[key][idx])
            
        grp = f.create_group('data')
        for key in lightcurves.keys():
            grp.create_dataset(key, data=lightcurves[key])
        
    return outfile
