#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

class sysfile():
    
    def __init__(self, sysfile):
        
        self.sysfile = sysfile
        
        return
        
    def read(self, systematic):
        
        if systematic == 'camtrans':
            
            with h5py.File(self.sysfile, 'r') as f:
                grp = f['data/'+systematic]
                grid = grp.attrs['grid']
                if grid == 'polar':
                    grid = PolarGrid(grp.attrs['nx'], grp.attrs['ny'])
                    idx = grp['idx'].value
                    value = grp['value'].value
                    
        return grid, idx, value
        
    def visualize(self, systematic):
        
        if systematic not in ['camtrans', 'intrapix', 'skytrans']:
            raise ValueError('invalid systematic %r' % (systematic,))
        
        if systematic == 'camtrans':
            
            grid, idx, value = self.read(systematic)
                
            array = grid.put_values_on_grid(value, idx, np.nan)
            
            plt.imshow(array, aspect='auto')
            plt.colorbar()
            plt.show()
        
        return

obj = sysfile('/data2/talens/3mEast/LBtests/reformat.hdf5')
obj.visualize('camtrans')
