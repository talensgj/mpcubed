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
    
    data = ['/data2/talens/2015Q2/LPE/tmp/red0_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5',
            '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5']
    
    for filename in data:
    
        with h5py.File(filename, 'r') as f:
            try:
                grp = f['data/' + ascc]
            except:
                continue
            else:
                jdmid = grp['jdmid'].value
                mag = grp['mag0'].value
            
        plt.plot(jdmid, mag, '.')
    
    ax = plt.subplot(212, sharex=ax, sharey=ax)
    
    data = ['/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPC/red0_vmag_2015Q2LPC.hdf5',
            '/data2/talens/2015Q2/LPW/red0_vmag_2015Q2LPW.hdf5']
    
    for filename in data:
    
        with h5py.File(filename, 'r') as f:
            try:
                grp = f['data/' + ascc]
            except:
                continue
            else:
                jdmid = grp['jdmid'].value
                mag = grp['mag0'].value

        plt.plot(jdmid, mag, '.')

    plt.ylim(1.5, -1.5)
    plt.show()
    
    return

def lc_diff(file1, file2, file3, ascc):
    
    with h5py.File(file1, 'r') as f:
        grp = f['data/' + ascc]
        jdmid1 = grp['jdmid'].value
        mag1 = grp['mag0'].value
        
    with h5py.File(file2, 'r') as f:
        grp = f['data/' + ascc]
        jdmid2 = grp['jdmid'].value
        mag2 = grp['mag0'].value
        
    with h5py.File(file3, 'r') as f:
        grp = f['data/' + ascc]
        jdmid3 = grp['jdmid'].value
        mag3 = grp['mag0'].value
    
    ax = plt.subplot(311)
    plt.title('ASCC {}'.format(ascc))
    plt.plot(jdmid1, mag1 - np.nanmedian(mag1), '.')
    
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jdmid2, mag2 - np.nanmedian(mag2), '.')
    
    plt.subplot(313, sharex=ax, sharey=ax)
    plt.plot(jdmid3, mag3 - np.nanmedian(mag3), '.')
    
    plt.show()
    
    return
    

def main():
    
    file1 = '/data2/talens/2015Q2/LPE/tmp/red0_2015Q2LPE.hdf5'
    file2 = '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5'
    file3 = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    
    
    #lc_diff(file1, file2, file3, '803743')
    #lc_diff(file1, file2, file3, '803977')
    #lc_diff(file1, file2, file3, '804411')
    #lc_diff(file1, file2, file3, '805622')
    #lc_diff(file1, file2, file3, '805861')
    #lc_diff(file1, file2, file3, '806455')
    
    lc_diff(file1, file2, file3, '894573')
    lc_diff(file1, file2, file3, '891488')
    lc_diff(file1, file2, file3, '714995')
    lc_diff(file1, file2, file3, '344250')
    
    #lightcurve('714995')
    
    return

if __name__ == '__main__':
    main()
