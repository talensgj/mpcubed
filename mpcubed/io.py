# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:23:14 2017

@author: talens
"""

import os

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

from .coordinates import grids

###############################################################################
### Code for reading files.
###############################################################################

class PhotFile(object):
    
    def __init__(self, filename):
        
        self.filename = filename
    
        return
        
    def read_global(self):
        
        with h5py.File(self.filename, 'r') as f:
            
            data = f['global'].attrs.items()
        
        return data
        
    def read_stars(self, fields=None, grpname='header_table'):
        
        stars = dict()        
            
        with h5py.File(self.filename, 'r') as f:
            
            grp = f[grpname]

            if fields is None:
                fields = grp.keys()
                
            for field in fields:
            
                if field in grp.keys():
                    stars[field] = grp[field].value
                else:
                    print 'Warning: skipping field {}, field not found.'.format(field)
                    
        return stars
    
    def _dsets2recarray(self, grp):
        
        names = grp.keys()
        arrays = [grp[key].value for key in names]
        
        return np.rec.fromarrays(arrays, names=names)
    
    def read_lightcurves(self, ascc=None, fields=None, perstar=True, grpname='data'):
        
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
            
            grp = f[grpname]
            
            if not hasattr(self, 'ascc0'):
                self.ascc0 = set(grp.keys())            
            
            for i in range(nstars):
                
                if ascc[i] in self.ascc0:
                    try:
                        curves[ascc[i]] = grp[ascc[i]].value
                    except: 
                        curves[ascc[i]] = self._dsets2recarray(grp[ascc[i]])
                else:
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
        
    def read_statistics(self, mode):
        """ Read the quality statistics from the header.
        
        Args:
            mode (str): Either 'spatial' or 'temporal' depending
                on which quality parameters should be read.
                
        Returns:
            niter (int): The number of iterations performed.
            chisq (float): The chi-squared value achieved.
            npoints (int): The number of datapoints in the fit.
            npars (int): The number of model parameters used.
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['header/' + mode]
            niter = grp['niter'].value
            chisq = grp['chisq'].value
            npoints = grp['npoints'].value
            npars = grp['npars'].value
            
        return niter, chisq, npoints, npars
        
    def read_magnitudes(self):
        """ Read the fitted magnitudes.
        
        Returns:
            ascc (str): The ascc numbers of the stars.
            vmag (float): The catalogue magnitude of the stars.
            mag (float): The fitted magnitude of the stars.
            sigma (float): The extra error term. Set to None if absent.
            nobs (int): The number of datapoints used to compute each magnitude.
            
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/magnitudes']
            ascc = grp['ascc'].value
            vmag = grp['vmag'].value
            mag = grp['mag'].value
            nobs = grp['nobs'].value
            
            try: sigma = grp['sigma'].value
            except: sigma = None
        
        return ascc, vmag, mag, sigma, nobs
        
    def read_trans(self):
        """ Read the fitted transmission.
        
        Returns:
            pg: A polar grid instance corresponding to the transmision.
            trans (float): The transmission values on the grid.
            nobs (int): The number of datapoints used to compute each
                transmission.
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/trans']
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            trans = grp['trans'].value
            nobs = grp['nobs'].value
            
            nx = grp.attrs['nx']
            ny = grp.attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        trans = pg.values2grid(idx1, idx2, trans, np.nan)
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan)

        return pg, trans, nobs
        
    def read_intrapix(self):
        """ Read the fitted intrapixel variations.
        
        Returns:
            pg: A polar grid instance corresponding to the intrapixel
                variations.
            sinx (float): The amplitudes of the sinx term on the grid. 
            cosx (float): The amplitudes of the cosx term on the grid. 
            siny (float): The amplitudes of the siny term on the grid. 
            cosy (float): The amplitudes of the cosy term on the grid. 
            nobs (float): The number of datapoints used to compute the
                amplitudes.
                
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/intrapix']
            idx1 = grp['idx1'].value
            idx2 = grp['idx2'].value
            sinx = grp['sinx'].value
            cosx = grp['cosx'].value
            siny = grp['siny'].value
            cosy = grp['cosy'].value
            nobs = grp['nobs'].value
            
            nx = grp.attrs['nx']
            ny = grp.attrs['ny']
            
        pg = grids.PolarGrid(nx, ny)
        sinx = pg.values2grid(idx1, idx2, sinx, np.nan)
        cosx = pg.values2grid(idx1, idx2, cosx, np.nan)
        siny = pg.values2grid(idx1, idx2, siny, np.nan)
        cosy = pg.values2grid(idx1, idx2, cosy, np.nan)
        nobs = pg.values2grid(idx1, idx2, nobs, np.nan)
        
        return pg, sinx, cosx, siny, cosy, nobs
    
    def read_clouds(self):
        """ Read the fitted clouds.
        
        Returns:
            hg: A halpix grid instance corresponding to the clouds.
            clouds: The cloud values computed on the grid at each instant.
            sigma: The extra error term. Set to None if absent.
            nobs: The number of datapoints used to compute the clouds.
            lstmin: The lowest lstseq for which clouds were computed.
            lstmax: The highest lstseq for which clouds were computed. 
        
        """
        
        with h5py.File(self.sysfile, 'r') as f:
            
            grp = f['data/clouds']
            idx = grp['idx'].value
            lstseq = grp['lstseq'].value
            clouds = grp['clouds'].value
            nobs = grp['nobs'].value
            
            try: sigma = grp['sigma'].value
            except: sigma = None
            
            grid = grp.attrs['grid']
            nx = grp.attrs['nx']
            lstmin = grp.attrs['lstmin']
            lstmax = grp.attrs['lstmax']
            lstlen = grp.attrs['lstlen']
    
        lstseq = lstseq - lstmin
    
        if (grid == 'healpix'):
            hg = grids.HealpixGrid(nx)
        elif (grid == 'polarea'):
            hg = grids.PolarEAGrid(nx)
        else:
            print 'Unknown grid.'
            exit()
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = clouds
        clouds = tmp
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = nobs
        nobs = tmp
        
        if sigma is not None:
            tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
            tmp[idx, lstseq] = sigma
            sigma = tmp

        return hg, clouds, sigma, nobs, lstmin, lstmax
        
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
                
                grp = f.create_group('data/{}'.format(stars['ascc'][j]))                
                for key in tmp.dtype.names:
                    grp.create_dataset(key, data=tmp[key])

    idx, = np.where(stars['nobs'] > 0)
    with h5py.File(filename) as f:
        
        grp = f.create_group('header')
        
        for key, value in stars.iteritems():
            grp.create_dataset(key, data = value[idx])

    return
    