#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from .. import IO

class fLCplot(object):
    """ Plot data from an fLC file.
    
    Attributes:
        fLCfile (str): The fLC file to read the data from.
    
    """
    
    def __init__(self, fLCfile):
        """ Initialize a plotter of fLC files.
        
        Args:
            fLCfile (str): The fLC file to read the data from.
            
        """
        
        self.fLCfile = fLCfile
        
        return
        
    def plot_coverage(self):
        """ Plot the sky coverage of the data in the file.""" 
        
        f = IO.fLCfile(self.fLCfile)
        ra, dec, nobs = f.read_header(['ra', 'dec', 'nobs'])
        
        plt.figure(figsize = (16, 8))
        
        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0], xticks = np.linspace(0, 360, 13), yticks = np.linspace(-80, 80, 9))

        im = plt.scatter(ra, dec, c = nobs, cmap = viridis)
        plt.xlim(0, 360)
        plt.ylim(-90, 90)
        plt.xlabel('Hour Angle [degrees]')
        plt.ylabel('Declination [degrees]')
        
        cax = plt.subplot(gs[1,1])
        cb = plt.colorbar(im, cax = cax)
        cb.set_label('nobs')
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        return

    def plot_star(self, ascc):
        """ Plot an overview of the lightcurve for a particular star.
        
        Args:
            ascc (str): The star to plot.
        
        """
        
        f = IO.fLCfile(self.fLCfile)
        lc = f.read_star(ascc)
        
        here = (lc['flag'] > 0)
        
        offset = np.floor(lc['jdmid'][0])
        lc['jdmid'] = lc['jdmid'] - offset
        
        baseline = np.ptp(lc['jdmid'])
        xmin = np.amin(lc['jdmid']) - .05*baseline
        xmax = np.amax(lc['jdmid']) + .05*baseline
        
        plt.figure(figsize=(16,8))
        
        ax = plt.subplot(311)
        plt.title(r'ASCC %s'%(ascc))
        plt.plot(lc['jdmid'], lc['flux0'], '.', label='flux0')
        plt.plot(lc['jdmid'], lc['flux1'], '.', label='flux1')
        plt.ylabel('Flux')
        plt.legend()
        
        plt.subplot(312, sharex=ax)
        plt.plot(lc['jdmid'], lc['flux0']/lc['flux1'], '.', label='flux0/flux1')
        plt.plot(lc['jdmid'], (lc['peak'] - lc['sky'])/lc['flux0'], '.', label='peak/flux0')
        plt.plot(lc['jdmid'], (lc['peak'] - lc['sky'])/lc['flux1'], '.', label='peak/flux1')
        plt.ylim(0,1)
        plt.ylabel('Fraction')
        plt.legend()
        
        plt.subplot(313, sharex=ax)
        plt.plot(lc['jdmid'], lc['sky'], '.')
        plt.xlim(xmin, xmax)
        plt.xlabel('Time [JD-%i]'%offset)
        plt.ylabel('Sky')
        
        plt.tight_layout()
        plt.show()
        plt.close()

        return

    def plot_as_array(self, ascc, field):
        """ Plot the data of a particular star as an array.
        
        This function is intended for merged fLC files and
        will plot each night as a row in a figure where the rows
        are the lstday and the columns are the lstidx.
        
        Args:
            ascc (str): The star to plot data from.
            field (str): The lightcurve field to plot.
        
        """
        
        f = IO.fLCfile(self.fLCfile)
        lc = f.read_star(ascc)

        lstday = (lc['lstseq'] // 13500)
        lstidx = (lc['lstseq'] % 13500)
                
        lstday = lstday - np.amin(lstday)
                
        nx = np.amax(lstday) + 1
        ny = 13500
                
        image = np.full((nx, ny), fill_value = np.nan)
        image[lstday, lstidx] = lc[field]
        
        plt.figure(figsize = (16, 8))
        
        plt.suptitle('ASCC %s'%ascc, size='xx-large')
        
        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0])
        
        im = plt.imshow(image, aspect = 'auto', cmap = viridis)
        plt.xlim(np.amin(lstidx) - .5, np.amax(lstidx) + .5)
        plt.xlabel('lstidx')
        plt.ylabel('lstday')
        
        cax = plt.subplot(gs[1,1])
        cb = plt.colorbar(im, cax = cax)
        cb.set_label(field)
        
        plt.tight_layout()
        plt.show()
        plt.close()
        
        return
