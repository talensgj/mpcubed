#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package import IO
from filters import masc_harmonic

def trans():

    f = IO.SysFile('/data2/talens/2015Q2/LPW/sys0_201506BLPW.hdf5')
    pg, trans1, nobs1 = f.read_trans()
    
    f = IO.SysFile('/data2/talens/2015Q2/LPW/sys0_vmag_201506BLPW.hdf5')
    pg, trans2, nobs2 = f.read_trans()
    
    plt.imshow((trans1 - trans2).T, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
    plt.colorbar()
    plt.show()
    plt.close()
    
    return

def lightcurve(ascc):
    
    ax = plt.subplot(211)
    
    data = ['/data2/talens/2015Q2/LPE/red0_201506ALPE.hdf5',
            '/data2/talens/2015Q2/LPE/red0_201506BLPE.hdf5',
            '/data2/talens/2015Q2/LPC/red0_201506ALPC.hdf5',
            '/data2/talens/2015Q2/LPC/red0_201506BLPC.hdf5',
            '/data2/talens/2015Q2/LPW/red0_201506ALPW.hdf5',
            '/data2/talens/2015Q2/LPW/red0_201506BLPW.hdf5']
    
    for filename in data:
    
        with h5py.File(filename, 'r') as f:
            try:
                lc = f['data/' + ascc].value
            except:
                continue
            
        plt.plot(lc['jdmid'], lc['mag0'], '.')
    
    ax = plt.subplot(212, sharex=ax, sharey=ax)
    
    data = ['/data2/talens/2015Q2/LPE/red0_vmag_201506ALPE.hdf5',
            '/data2/talens/2015Q2/LPE/red0_vmag_201506BLPE.hdf5',
            '/data2/talens/2015Q2/LPC/red0_vmag_201506ALPC.hdf5',
            '/data2/talens/2015Q2/LPC/red0_vmag_201506BLPC.hdf5',
            '/data2/talens/2015Q2/LPW/red0_vmag_201506ALPW.hdf5',
            '/data2/talens/2015Q2/LPW/red0_vmag_201506BLPW.hdf5']
    
    for filename in data:
    
        with h5py.File(filename, 'r') as f:
            try:
                lc = f['data/' + ascc].value
            except:
                continue

        plt.plot(lc['jdmid'], lc['mag0'], '.')

    plt.ylim(.2, -.2)
    plt.show()
    
    return

def main():
    
    lightcurve('1480383')
    
    return

if __name__ == '__main__':
    main()
