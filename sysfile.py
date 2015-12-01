#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid

class SysFile():
    
    def __init__(self, sysfile):
        
        self.sysfile = sysfile
        
        return
        
    def read_magnitudes(self):
        
        return
        
    def read_camtrans(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            idx = f['data/camtrans/idx'].value
            z = f['data/camtrans/z'].value
            nobs = f['data/camtrans/nobs'].value
            
            nx = f['data/camtrans'].attrs['nx']
            ny = f['data/camtrans'].attrs['ny']
            
        pg = PolarGrid(nx, ny)
        z = pg.put_values_on_grid(z, idx, np.nan)
        nobs = pg.put_values_on_grid(nobs, idx)

        return z, nobs
        
    def read_intrapix(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            idx = f['data/intrapix/idx'].value
            a = f['data/intrapix/a'].value
            b = f['data/intrapix/b'].value
            c = f['data/intrapix/c'].value
            d = f['data/intrapix/d'].value
            nobs = f['data/intrapix/nobs'].value
            
            nx = f['data/intrapix'].attrs['nx']
            ny = f['data/intrapix'].attrs['ny']
            
        pg = PolarGrid(nx, ny)
        a = pg.put_values_on_grid(a, idx, np.nan)
        b = pg.put_values_on_grid(b, idx, np.nan)
        c = pg.put_values_on_grid(c, idx, np.nan)
        d = pg.put_values_on_grid(d, idx, np.nan)
        nobs = pg.put_values_on_grid(nobs, idx)
        
        return a, b, c, d, nobs
    
    def read_skytrans(self):
        
        with h5py.File(self.sysfile, 'r') as f:
            
            idx = f['data/skytrans/idx'].value
            lstseq = f['data/skytrans/lstseq'].value
            s = f['data/skytrans/s'].value
            nobs = f['data/skytrans/nobs'].value
            
            nx = f['data/skytrans'].attrs['nx']
            lstmin = f['data/skytrans'].attrs['lstmin']
            lstlen = f['data/skytrans'].attrs['lstlen']
    
        hg = HealpixGrid(nx)
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = s
        s = tmp
        
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = nobs
        nobs = tmp

        return s, nobs, lstmin

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    obj = SysFile('/data2/talens/2015Q2/LPE/sys_201506BLPE.hdf5')
    z, nobs = obj.read_camtrans()
    a, b, c, d, nobs = obj.read_intrapix()
    s, nobs, lstmin = obj.read_skytrans()
    print lstmin
    exit()
    plt.subplot(111, aspect='auto')
    plt.imshow(nobs.T, aspect='auto', origin='lower', interpolation='None')
    plt.colorbar()
    plt.show()
    
    
    
