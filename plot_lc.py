#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt

from coordinate_grids import PolarGrid
from index_functions import index_statistics
import os
import glob

class PlotLC():
    
    def __init__(self, path, date, cameras = ['LPN', 'LPE', 'LPS', 'LPW', 'LPC'], colors = ['b', 'r', 'g', 'y', 'c']):

        self.cameras = cameras
        self.colors = colors
        self.raw_files = [os.path.join(path, 'fLC_'+date+camera+'.hdf5') for camera in cameras]
        self.red_files = [os.path.join(path, 'red_'+date+camera+'.hdf5') for camera in cameras]
        self.figure = False
        
    def plot_lightcurve(self, ascc, raw = True, binned = None):
        
        self.ascc = ascc

        if not self.figure:
            self.figure = plt.figure(figsize=(16,4))
        
        self.ax1 = plt.subplot(111)
        
        for filename, camera, color in zip(self.raw_files, self.cameras, self.colors):
            with h5py.File(filename, 'r') as f:
                
                try:
                    si = f['header/'+ascc]
                except:
                    pass
                else:
                    lc = f['data/'+ascc]
                    jd_ref = np.floor(lc['jdmid'])
                    
                    if raw:
                        plt.plot(lc['jdmid']-jd_ref, lc['flux0'], '.', label=camera, c=color, zorder=0)
                    elif binned is None:
                        binned = 50
                    
                    if binned is not None:
                        
                        jd_bin = index_statistics(lc['lstidx']//binned, lc['jdmid']-jd_ref) 
                        flux_bin = index_statistics(lc['lstidx']//binned, lc['flux0']) 
                        eflux_bin = index_statistics(lc['lstidx']//binned, lc['flux0'], statistic='std')/np.sqrt(index_statistics(lc['lstidx']//binned, lc['flux0'], statistic='count'))
                        
                        plt.errorbar(jd_bin, flux_bin, yerr=eflux_bin, fmt='o', ecolor='k', c=color, zorder=1, label=camera) 
                        
                    
                    plt.xlim(.25, .75)
                    plt.legend(loc=2)
                    plt.xticks(size='large')
                    plt.yticks(size='large')
                    plt.xlabel('Time [JD-%i]'%jd_ref[0], size='x-large')
                    plt.ylabel('Flux', size='x-large')

    def add_xytrans(self):
        
        self.figure.set_size_inches(16,8, forward=True)
        self.ax1.change_geometry(2,1,1)
        self.ax2 = plt.subplot(212, sharex=self.ax1)
        
        for raw_file, red_file in zip(self.raw_files, self.red_files):
            
            with h5py.File(red_file, 'r') as f:
                trans = f['data/'+self.ascc+'/trans0'].value
            
            with h5py.File(raw_file, 'r') as f:
                
                try:
                    si = f['header/'+self.ascc]
                except:
                    pass
                else:
                    lc = f['data/'+self.ascc]
                    jd_ref = np.floor(lc['jdmid'])
                    plt.plot(lc['jdmid']-jd_ref, trans, '.')#, label=camera, c=color, zorder=0)
        
        plt.xlim(.25, .75)
        plt.xticks(size='large')
        plt.yticks(size='large')
        plt.xlabel('Time [JD-%i]'%jd_ref[0], size='x-large')
        plt.ylabel('Trans', size='x-large')
        
    def add_skytrans(self):
        
        self.figure.set_size_inches(16,12, forward=True)
        self.ax1.change_geometry(3,1,1)
        self.ax2.change_geometry(3,1,2)
        self.ax3 = plt.subplot(313, sharex=self.ax1)
        
        for raw_file, red_file in zip(self.raw_files, self.red_files):
            
            with h5py.File(red_file, 'r') as f:
                trans = f['data/'+self.ascc+'/cflux0'].value
            
            with h5py.File(raw_file, 'r') as f:
                
                try:
                    si = f['header/'+self.ascc]
                except:
                    pass
                else:
                    lc = f['data/'+self.ascc]
                    jd_ref = np.floor(lc['jdmid'])
                    plt.plot(lc['jdmid']-jd_ref, trans, '.')#, label=camera, c=color, zorder=0)
        
        plt.xlim(.25, .75)
        plt.xticks(size='large')
        plt.yticks(size='large')
        plt.xlabel('Time [JD-%i]'%jd_ref[0], size='x-large')
        plt.ylabel('cFlux', size='x-large')
        
    def show(self):
        plt.tight_layout()
        plt.show()

obj = PlotLC('/data2/talens/Jul2015', '20150714')
obj.plot_lightcurve('807144')
#obj.add_xytrans()
#obj.add_skytrans()
obj.show()
