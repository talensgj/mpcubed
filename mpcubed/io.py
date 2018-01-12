# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:23:14 2017

@author: talens
"""

import os

import h5py
import numpy as np

from .calibration import grids

###############################################################################
### Code for writing files.
###############################################################################

def write_calibration(filename, settings, spatial, temporal, magnitudes, trans, intrapix, clouds):
    
    with h5py.File(filename) as f:

        # Write the header.
        hdr = f.create_group('header')
        
        hdr.create_dataset('filelist', data=settings['filelist'])
        hdr.attrs['station'] = settings['station']
        hdr.attrs['camera'] = settings['camera']
        
        hdr.attrs['alt0'] = settings['alt0']
        hdr.attrs['az0'] = settings['az0']
        hdr.attrs['th0'] = settings['th0']
        hdr.attrs['x0'] = settings['x0']
        hdr.attrs['y0'] = settings['y0']
        
        hdr.attrs['outer_maxiter'] = settings['outer_maxiter']
        hdr.attrs['inner_maxiter'] = settings['inner_maxiter']
        hdr.attrs['sigmas'] = True
        hdr.attrs['dtol'] = settings['dtol']
       
        # Write the spatial quality checks.
        idx = spatial['n']
        
        subhdr = hdr.create_group('spatial')
        subhdr.create_dataset('niter', data=spatial['niter'][idx])
        subhdr.create_dataset('chisq', data=spatial['chisq'][idx])
        subhdr.create_dataset('npoints', data=spatial['npoints'][idx])
        subhdr.create_dataset('npars', data=spatial['npars'][idx])
        
        # Write the temporal quality checks.
        idx = temporal['q']
        
        subhdr = hdr.create_group('temporal')
        subhdr.create_dataset('niter', data=temporal['niter'][idx])
        subhdr.create_dataset('chisq', data=temporal['chisq'][idx])
        subhdr.create_dataset('npoints', data=temporal['npoints'][idx])
        subhdr.create_dataset('npars', data=temporal['npars'][idx])
        
        # Write the data.
        grp = f.create_group('data')
        
        # Write the magnitudes.
        subgrp = grp.create_group('magnitudes')
        subgrp.create_dataset('ascc', data=magnitudes['ascc'])
        subgrp.create_dataset('vmag', data=magnitudes['vmag'], dtype='float32')
        subgrp.create_dataset('nobs', data=magnitudes['nobs'])
        subgrp.create_dataset('mag', data=magnitudes['mag'], dtype='float32')
        subgrp.create_dataset('sigma', data=magnitudes['sigma'], dtype='float32')

        # Write the camera transmission.
        idx1, idx2 = np.where(trans['nobs'] > 0)  
        
        subgrp = grp.create_group('trans')          
        subgrp.create_dataset('idx1', data=idx1, dtype='uint32')
        subgrp.create_dataset('idx2', data=idx2, dtype='uint32')          
        subgrp.create_dataset('nobs', data=trans['nobs'][idx1,idx2])
        subgrp.create_dataset('trans', data=trans['trans'][idx1,idx2], dtype='float32')
        
        subgrp.attrs['grid'] = trans['gridtype']
        subgrp.attrs['nx'] = trans['num_k']
        subgrp.attrs['ny'] = trans['num_n']
        
        # Write the intrapixel variations.
        idx1, idx2 = np.where(intrapix['nobs'] > 0) 

        subgrp = grp.create_group('intrapix')           
        subgrp.create_dataset('idx1', data=idx1, dtype='uint32')
        subgrp.create_dataset('idx2', data=idx2, dtype='uint32')
        subgrp.create_dataset('nobs', data=intrapix['nobs'][idx1,idx2])
        subgrp.create_dataset('sinx', data=intrapix['amplitudes'][idx1,idx2,0], dtype='float32')
        subgrp.create_dataset('cosx', data=intrapix['amplitudes'][idx1,idx2,1], dtype='float32')
        subgrp.create_dataset('siny', data=intrapix['amplitudes'][idx1,idx2,2], dtype='float32')
        subgrp.create_dataset('cosy', data=intrapix['amplitudes'][idx1,idx2,3], dtype='float32')
        
        subgrp.attrs['grid'] = intrapix['gridtype']
        subgrp.attrs['nx'] = intrapix['num_l']
        subgrp.attrs['ny'] = intrapix['num_n']
        
        # Write the sky transmission.
        idx, lstseq = np.where(clouds['nobs'] > 0)

        subgrp = grp.create_group('clouds')            
        subgrp.create_dataset('idx', data=idx, dtype='uint32')
        subgrp.create_dataset('lstseq', data=clouds['lstmin'] + lstseq, dtype='uint32')
        subgrp.create_dataset('nobs', data=clouds['nobs'][idx, lstseq])
        subgrp.create_dataset('clouds', data=clouds['clouds'][idx, lstseq], dtype='float32')
        subgrp.create_dataset('sigma', data=clouds['sigma'][idx, lstseq], dtype='float32')
        
        subgrp.attrs['grid'] = clouds['gridtype']
        subgrp.attrs['nx'] = clouds['num_q']
        subgrp.attrs['lstmin'] = clouds['lstmin']
        subgrp.attrs['lstmax'] = clouds['lstmax']
        subgrp.attrs['lstlen'] = clouds['lstlen']

    return

def write_reduced(filename, settings, stars, lightcurves):
    
    with h5py.File(filename) as f:
            
        # Write the global information.
        grp = f.create_group('global')
        grp.attrs['station'] = settings['station']
        grp.attrs['camera'] = settings['camera']
        grp.attrs['exptime'] = 6.4 # Hardcoded ...
        grp.attrs['naper'] = settings['naper']
        grp.attrs['aper0'] = settings['aper0']
        grp.attrs['aper1'] = settings['aper1']
        grp.attrs['skyrad0'] = settings['skyrad0']
        grp.attrs['skyrad1'] = settings['skyrad1']
        
        grp.create_dataset('filelist', data=settings['filelist'])
        grp.create_dataset('aversion', data=settings['aversion'])
        grp.create_dataset('rversion', data=settings['rversion'])
        grp.create_dataset('cversion', data=settings['cversion'])
        grp.attrs['pversion'] = '1.0.0' # Hardcoded ...
        
        # Write the header_table.
        idx, = np.where(stars['nobs'] > 0)            
        
        grp = f.create_group('header_table')
        for key in stars.keys():
            grp.create_dataset(key, data=stars[key][idx])
            
        # Write the reduced lightcurves.
        grp = f.create_group('data')
        for key in lightcurves.keys():
            grp.create_dataset(key, data=lightcurves[key])
            
    return

def write_boxlstsq(filename, ascc, chisq0, boxpars, criteria, freq, dchisq):
    
    with h5py.File(filename) as f:
        
        # Write the header.
        grp = f.create_group('header')
        grp.create_dataset('ascc', data=ascc)
        grp.create_dataset('chisq0', data=chisq0, dtype='float32')

        for key in boxpars.dtype.names:
            grp.create_dataset(key, data=boxpars[key])
        
        for key in criteria.dtype.names:
            grp.create_dataset(key, data=criteria[key])

        # Write the data.
        grp = f.create_group('data')
        grp.create_dataset('freq', data=freq)
        grp.create_dataset('dchisq', data=dchisq, dtype='float32')
    
    return

###############################################################################
### Code for reading files.
###############################################################################

class PhotFile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
    
        return
        
    def read_global(self):
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['global']
            
            data = grp.attrs.items()
            data = dict(data)
            
            for key in grp.keys():
                data[key] = grp[key].value
        
        return data
        
    def read_stars(self, fields=None):
        
        stars = dict()        
            
        with h5py.File(self.filename, 'r') as f:
            
            if 'stars' in f.keys():
                grp = f['stars']
            elif 'header_table' in f.keys():
                grp = f['header_table']
            elif 'header' in f.keys():
                grp = f['header']
            else:
                raise IOError('No valid stars field found.')

            if fields is None:
                fields = grp.keys()
                
            for field in fields:
            
                if field in grp.keys():
                    stars[field] = grp[field].value
                else:
                    print 'Warning: skipping field {}, field not found.'.format(field)
                    
        return stars
    
    def read_lightcurves(self, ascc=None, fields=None, perstar=True, verbose=True):
        
        onestar = False        
        
        if ascc is None:
            stars = self.read_stars(['ascc'])
            ascc = stars['ascc']
            
        elif isinstance(ascc, basestring):
            onestar = True
            ascc = [ascc]
       
        nstars = len(ascc)
        curves = dict()
        nobs = np.zeros(nstars, dtype='int')
            
        # Read the data.
        with h5py.File(self.filename, 'r') as f:
            
            try:
                grp = f['data'] # La Palma
            except:
                grp = f['lightcurves'] # bRing, La Silla
            
            if not hasattr(self, 'ascc0'):
                self.ascc0 = set(grp.keys())            
            
            for i in range(nstars):
                
                if ascc[i] in self.ascc0:
                    curves[ascc[i]] = grp[ascc[i]].value
                    
                else:
                    if verbose:
                        print 'Warning: skipping star {}, star not found.'.format(ascc[i])
                    continue
                
                nobs[i] = len(curves[ascc[i]])
                    
        if not curves:
            return curves
            
        # Select specified fields.
        if fields is not None:
            
            for i in range(nstars):
                curves[ascc[i]] = curves[ascc[i]][fields]
                    
        # Combine lightcurves.
        if not perstar:
            
            strides = np.append(0, np.cumsum(nobs))
            tmp = np.recarray(strides[-1], dtype=curves[ascc[0]].dtype)
            
            for i in range(nstars):
                tmp[strides[i]:strides[i+1]] = curves[ascc[i]]
                        
            curves = tmp                        
        
        if onestar:
            return curves[ascc[0]]             
            
        return curves
        
class SysFile(object):
    """ Read data from systematics files.
    
    Attributes:
        sysfile (str): The full path to the file.
        
    """
    
    
    def __init__(self, sysfile):
        """ Initialize a reader of systematics files.
        
        Args:
            sysfile (str): The full path to the file.
            
        """
        
        self.sysfile = sysfile
        
        return
        
    def read_header(self):
        """ Read all header attributes.
        
        Returns:
            data: A list of attribute (key, value) pairs.
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            data = f['header'].attrs.items()
        
        return data
        
    def read_pointing(self):
        """ Read the pointing associated with the systematics.
        
        Returns:
            alt0 (float): The pointing altitude in degrees.
            az0 (float): The pointing azimuth in degrees. 
            th0 (float): The pointing orientation in degrees.
            x0 (float): The pointing center in pixels.
            y0 (float): The pointing center in pixels.
            
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            alt0 = f['header'].attrs['alt0']
            az0 = f['header'].attrs['az0']
            th0 = f['header'].attrs['th0']
            x0 = f['header'].attrs['x0']
            y0 = f['header'].attrs['y0']
        
        return alt0, az0, th0, x0, y0
        
    def read_spatial(self):
        """ Read the spatial quality statistics.
                
        Returns:
            spatial (dict): A dictionary containing the quality statistics.
        
        """
        
        spatial = dict()
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['header/spatial']
            
            spatial['niter'] = grp['niter'].value
            spatial['chisq'] = grp['chisq'].value
            spatial['npoints'] = grp['npoints'].value
            spatial['npars'] = grp['npars'].value
            
        return spatial
    
    def read_temporal(self):
        """ Read the temporal quality statistics.
                
        Returns:
            temporal (dict): A dictionary containing the quality statistics.
        
        """
        
        temporal = dict()
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['header/temporal']
            
            temporal['niter'] = grp['niter'].value
            temporal['chisq'] = grp['chisq'].value
            temporal['npoints'] = grp['npoints'].value
            temporal['npars'] = grp['npars'].value
            
        return temporal
        
    def read_magnitudes(self):
        """ Read the best-fit magnitudes.
        
        Returns:
            magnitudes (dict): Containing: ascc, vmag, nobs, mag and sigma.
            
        """
        
        magnitudes = dict()
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/magnitudes']
            
            magnitudes['ascc'] = grp['ascc'].value
            magnitudes['vmag'] = grp['vmag'].value
            magnitudes['nobs'] = grp['nobs'].value
            magnitudes['mag'] = grp['mag'].value
            magnitudes['sigma'] = grp['sigma'].value
            
        return magnitudes
        
    def read_trans(self):
        """ Read the fitted transmission.
        
        Returns:
            grid: A polar grid instance corresponding to the transmission map.
            trans (dict): Containing: grid, num_k, num_n, nobs, trans.
        
        """

        trans = dict()
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/trans']
            
            trans['gridtype'] = grp.attrs['grid']
            trans['num_k'] = grp.attrs['nx']
            trans['num_n'] = grp.attrs['ny']
            
            grid = grids.PolarGrid(trans['num_k'], trans['num_n'])
            
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            trans['nobs'] = grid.values2grid(idx1, idx2, grp['nobs'].value)
            trans['trans'] = grid.values2grid(idx1, idx2, grp['trans'].value)
        
        return grid, trans
        
    def read_intrapix(self):
        """ Read the fitted intrapixel variations.
        
        Returns:
            grid: A polar grid instance corresponding to the intrapixel
                variations.
            intrapix (dict): Containing grid, num_l, num_n, nobs, amplitudes/
                
        """
        
        intrapix = dict()
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/intrapix']
            
            intrapix['gridtype'] = grp.attrs['grid']
            intrapix['num_l'] = grp.attrs['nx']
            intrapix['num_n'] = grp.attrs['ny']
            
            grid = grids.PolarGrid(intrapix['num_l'], intrapix['num_n'])
            
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            
            sinx = grid.values2grid(idx1, idx2, grp['sinx'].value)
            cosx = grid.values2grid(idx1, idx2, grp['cosx'].value)
            siny = grid.values2grid(idx1, idx2, grp['siny'].value)
            cosy = grid.values2grid(idx1, idx2, grp['cosy'].value)
            
            intrapix['nobs'] = grid.values2grid(idx1, idx2, grp['nobs'].value)
            intrapix['amplitudes'] = np.stack([sinx, cosx, siny, cosy], axis=-1)
            
        return grid, intrapix
    
    def read_clouds(self):
        """ Read the fitted clouds.
        
        Returns:
            grid: A healpix grid instance corresponding to the clouds.
            clouds (dict): Containing: gridtype, num_q, lstmin, lstmax, lstlen,
                nobs, clouds, sigma.
        
        """
        
        clouds = dict()
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/clouds']
            
            clouds['gridtype'] = grp.attrs['grid']
            clouds['num_q'] = grp.attrs['nx']
            clouds['lstmin'] = grp.attrs['lstmin']
            clouds['lstmax'] = grp.attrs['lstmax']
            clouds['lstlen'] = grp.attrs['lstlen']
            
            grid = grids.HealpixGrid(clouds['num_q'])
            
            idx = grp['idx'].value
            lstseq = grp['lstseq'].value - clouds['lstmin']
            
            nobs_ = grp['nobs'].value
            clouds_ = grp['clouds'].value
            sigma_ = grp['sigma'].value
            
        clouds['nobs'] = np.full((grid.npix, clouds['lstlen']), fill_value=np.nan)
        clouds['nobs'][idx, lstseq] = nobs_
        
        clouds['clouds'] = np.full((grid.npix, clouds['lstlen']), fill_value=np.nan)
        clouds['clouds'][idx, lstseq] = clouds_
    
        clouds['sigma'] = np.full((grid.npix, clouds['lstlen']), fill_value=np.nan)
        clouds['sigma'][idx, lstseq] = sigma_
        
        return grid, clouds
        
class blsFile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
    
        return
        
    def read_header(self, fields=None):
        
        hdr = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['header']
            
            if fields is None:
                fields = grp.keys()            
            
            for field in fields:
                hdr[field] = grp[field].value
                
        return hdr
        
    def read_data(self, fields=None):
        
        data = dict()
        
        with h5py.File(self.filename, 'r') as f:
            
            grp = f['data']
            
            if fields is None:
                fields = grp.keys()
            
            for field in fields:
                data[field] = grp[field].value
        
        return data
        
###############################################################################
### Code for combining files.
###############################################################################
        
def verify_filelist(filelist):
    
    nfiles = len(filelist)
    
    args = []
    for i in range(nfiles):
        if os.path.isfile(filelist[i]):
            args.append(i)
        else:
            print filelist[i], 'does not exist.'
    
    filelist = filelist[args]
    
    if len(filelist) == 0:
        print 'No valid files.'
        print 'exiting...'
        return filelist
    
    return filelist

def _index_files(filelist):
    
    nfiles = len(filelist)
    
    idx1 = np.array([], dtype='uint16')
    
    stars = dict()
    for i in range(nfiles):
        
        filename = filelist[i]        
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['header_table']
            
            for key in grp.keys():
                
                ascc_ = grp['ascc'].value
                
                if key not in stars.keys():
                    stars[key] = grp[key].value
                else:
                    stars[key] = np.append(stars[key], grp[key].value)

            grp = f['data']
            
            if len(ascc_) > 0:
                dtype = grp[ascc_[0]].dtype

        idx1 = np.append(idx1, np.repeat(i, len(ascc_)))
    
    ascc, args, idx2 = np.unique(stars['ascc'], return_index=True, return_inverse=True)
    nstars = len(ascc)    
    nobs = np.zeros((nfiles, nstars), dtype='uint32')
    stars['nobs'] = stars['nobs'].astype('uint32')
    nobs[idx1, idx2] = stars['nobs']
    
    for key in stars.keys():
        stars[key] = stars[key][args]
    
    return stars, nobs, dtype

def _read_curves(filelist, ascc, nobs, dtype):
    
    nfiles = len(filelist)
    nstars = len(ascc)
    
    strides = np.row_stack([nstars*[0], np.cumsum(nobs, axis=0)]).astype('int')
    curves = {ascc[i]:np.recarray(strides[-1,i], dtype=dtype) for i in range(nstars)}
    
    for i in range(nfiles):
        
        filename = filelist[i]
        
        with h5py.File(filename, 'r') as f:
            
            grp = f['data']
            
            for j in range(nstars):
                
                if (nobs[i,j] > 0):
                    
                    curves[ascc[j]][strides[i,j]:strides[i+1,j]] = grp[ascc[j]].value
                    
    return curves
    
def _read_global(filelist):
    
    nfiles = len(filelist)
    
    # Create arrays.
    attrdict = {}
    
    aversion = []
    rversion = []
    cversion = []
    
    alt0 = np.zeros(nfiles)
    az0 = np.zeros(nfiles)
    th0 = np.zeros(nfiles)
    x0 = np.zeros(nfiles)
    y0 = np.zeros(nfiles)
    
    exptime = np.zeros(nfiles)
    ccdtemp = np.zeros(nfiles)
    
    # Read the global groups.
    for i in range(nfiles):
        with h5py.File(filelist[i], 'r') as f:
            grp = f['global']
            
            if (i < 1):
                
                try: attrdict['station'] = grp.attrs['station']
                except: attrdict['station'] = grp.attrs['STATION']
            
                try: attrdict['camera'] = grp.attrs['camera']
                except: attrdict['camera'] = grp.attrs['CAMERA']
            
                try: attrdict['naper'] = grp.attrs['naper']
                except: attrdict['naper'] = grp.attrs['NAPER']
                
                try: attrdict['aper0'] = grp.attrs['aper0']
                except: attrdict['aper0'] = grp.attrs['APER0']
                
                try: attrdict['aper1'] = grp.attrs['aper1']
                except: attrdict['aper1'] = grp.attrs['APER1']
                
                try: attrdict['skyrad0'] = grp.attrs['skyrad0']
                except: attrdict['skyrad0'] = grp.attrs['SKYRAD0']
                
                try: attrdict['skyrad1'] = grp.attrs['skyrad1']
                except: attrdict['skyrad1'] = grp.attrs['SKYRAD1']
            
            try: alt0[i] = grp.attrs['alt0']
            except: alt0[i] = grp.attrs['ALT0']
            
            try: az0[i] = grp.attrs['az0']
            except: az0[i] = grp.attrs['AZ0']
            
            try: th0[i] = grp.attrs['th0']
            except: th0[i] = grp.attrs['TH0']
            
            try: x0[i] = grp.attrs['x0']
            except: x0[i] = grp.attrs['X0']
            
            try: y0[i] = grp.attrs['y0']
            except: y0[i] = grp.attrs['Y0']
            
            try: aversion.append(grp.attrs['aversion'])
            except: aversion.append(grp.attrs['AVERSION'])
            
            try: rversion.append(grp.attrs['rversion'])
            except: rversion.append(grp.attrs['RVERSION'])
            
            try: cversion.append(grp.attrs['cversion'])
            except: cversion.append(grp.attrs['CVERSION'])
            
            try: exptime[i] = grp.attrs['exptime']
            except: exptime[i] = grp.attrs['EXPTIME']
            
            try: ccdtemp[i] = grp.attrs['ccdtemp']
            except: ccdtemp[i] = grp.attrs['CCDTEMP']
    
    # Put result in dictionary
    arrdict = {}
    arrdict['alt0'] = alt0
    arrdict['az0'] = az0
    arrdict['th0'] = th0
    arrdict['x0'] = x0
    arrdict['y0'] = y0
    arrdict['aversion'] = np.array(aversion)
    arrdict['rversion'] = np.array(rversion)
    arrdict['cversion'] = np.array(cversion)
    arrdict['ccdtemp'] = ccdtemp
    arrdict['exptime'] = exptime
    
    return attrdict, arrdict

def make_baseline(filename, filelist, nsteps=1000):
    
    filelist = np.sort(filelist)    
    filelist = verify_filelist(filelist)
    
    if len(filelist) == 0:
        return
        
    # Read the combined stars field and index the files.
    stars, nobs, dtype = _index_files(filelist)    
    
    stars['lstsqmin'] = np.zeros(len(stars['ascc']), dtype='uint32')
    stars['lstsqmax'] = np.zeros(len(stars['ascc']), dtype='uint32')
    
    nstars = len(stars['ascc'])
    for i in range(0, nstars, nsteps):
        
        # Read the combined lightcurves for a group of stars.
        curves = _read_curves(filelist, stars['ascc'][i:i+nsteps], nobs[:,i:i+nsteps], dtype)
             
        # Write the combined lightcurves for a group of stars.
        with h5py.File(filename) as f:
            
            for j in range(i, i+len(stars['ascc'][i:i+nsteps])):
             
                tmp = curves[stars['ascc'][j]]
                
                stars['nobs'][j] = len(tmp)
                stars['lstsqmin'][j] = tmp['lstseq'][0]
                stars['lstsqmax'][j] = tmp['lstseq'][-1]
                
                f.create_dataset('data/{}'.format(stars['ascc'][j]), data=tmp)    

    # Write the combined "header_table" field.
    with h5py.File(filename) as f:
        
        grp = f.create_group('header_table')
        for key in stars.keys():
            grp.create_dataset(key, data=stars[key])
            
    # Write the global group.
    lstmin = np.amin(stars['lstsqmin'])
    lstmax = np.amax(stars['lstsqmax'])
    
    attrdict, arrdict = _read_global(filelist)
    
    with h5py.File(filename) as f:
        
        grp = f.create_group('global')
        
        grp.attrs['lstmin'] = lstmin
        grp.attrs['lstmax'] = lstmax
        
        grp.create_dataset('filelist', data=filelist)
        
        for key, value in attrdict.iteritems():
            grp.attrs[key] = value
            
        for key, value in arrdict.iteritems():
            grp.create_dataset(key, data = value)

    return

def make_quarter(filename, filelist, nsteps=1000):
    
    filelist = np.sort(filelist)
    
    # Write the global group.
    with h5py.File(filelist[0]) as f:
        data = f['global'].attrs.items()
        
    data = dict(data)
    
    with h5py.File(filename) as f:
        grp = f.create_group('global')
        grp.attrs['station'] = data['station']
        grp.attrs['camera'] = data['camera']
    
    # Merge the headers.
    stars, nobs, dtype = _index_files(filelist)
        
    # Merge the lightcurves.
    nstars = len(stars['ascc'])
    stars['jdmin'] = np.zeros(nstars)
    stars['jdmax'] = np.zeros(nstars)
    stars['lstseqmin'] = np.zeros(nstars, dtype='uint32')
    stars['lstseqmax'] = np.zeros(nstars, dtype='uint32')
    
    for i in range(0, nstars, nsteps):
        
        # Read the combined lightcurves for a group of stars.
        curves = _read_curves(filelist, stars['ascc'][i:i+nsteps], nobs[:,i:i+nsteps], dtype)
             
        # Write the combined lightcurves for a group of stars.
        with h5py.File(filename) as f:
            
            for j in range(i, i+len(stars['ascc'][i:i+nsteps])):
             
                tmp = curves[stars['ascc'][j]]
                
                stars['nobs'][j] = len(tmp)
    
                if (len(tmp) == 0): continue                
                
                stars['jdmin'][j] = tmp['jdmid'][0]
                stars['jdmax'][j] = tmp['jdmid'][-1]
                stars['lstseqmin'][j] = tmp['lstseq'][0]
                stars['lstseqmax'][j] = tmp['lstseq'][-1]
                
                f.create_dataset('data/{}'.format(stars['ascc'][j]), data=tmp)                

    idx, = np.where(stars['nobs'] > 0)
    with h5py.File(filename) as f:
        
        grp = f.create_group('header')
        
        for key, value in stars.iteritems():
            grp.create_dataset(key, data = value[idx])

    return
    
def merge_files(filetype, filename, filelist, nsteps=1000):
    
    if filetype == 'fLC':
        make_baseline(filename, filelist, nsteps)
    elif filetype == 'red':
        make_quarter(filename, filelist, nsteps)
    else:
        raise ValueError('Filetype {}, is not suported.'.format(filetype))
    
    return

def main():
    
    import argparse

    parser = argparse.ArgumentParser(description='Combine temporary lightcurve files.')
    parser.add_argument('filetype', type=str, default=['fLC', 'red'],
                        help='the type of lightcurve files we are merging')
    parser.add_argument('filename', type=str,
                        help='the name of the resulting file')
    parser.add_argument('filelist', type=str, nargs='+',
                        help='the files to combine')
    parser.add_argument('-n', '--nsteps', type=int, default=1000,
                        help='the number of stars to read in one go', dest='nsteps')
    
    args = parser.parse_args()
    
    merge_files(args.filetype, args.filename, args.filelist, args.nsteps)
    
    return

if __name__ == '__main__':
    
    main()
    