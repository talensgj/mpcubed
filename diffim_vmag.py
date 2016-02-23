#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'


from package import IO
import filters

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

def lc_diff_mc(file1, file2, file3, ascc):
    
    with h5py.File(file1, 'r') as f:
        grp = f['data/' + ascc]
        jdmid1 = grp['jdmid'].value
        lst1 = grp['lst'].value
        mag1 = grp['mag0'].value
        emag1 = grp['emag0'].value
        nobs1 = grp['nobs'].value
        
    select = (nobs1 == 50)
    jdmid1 = jdmid1[select]
    lst1 = lst1[select]
    mag1 = mag1[select]
    emag1 = emag1[select]
        
    with h5py.File(file2, 'r') as f:
        grp = f['data/' + ascc]
        jdmid2 = grp['jdmid'].value
        lst2 = grp['lst'].value
        mag2 = grp['mag0'].value
        emag2 = grp['emag0'].value
        nobs2 = grp['nobs'].value
        
    select = (nobs2 == 50)
    jdmid2 = jdmid2[select]
    lst2 = lst2[select]
    mag2 = mag2[select]
    emag2 = emag2[select]
        
    with h5py.File(file3, 'r') as f:
        grp = f['data/' + ascc]
        jdmid3 = grp['jdmid'].value
        lst3 = grp['lst'].value
        mag3 = grp['mag0'].value
        emag3 = grp['emag0'].value
        nobs3 = grp['nobs'].value
        
    select = (nobs3 == 50)
    jdmid3 = jdmid3[select]
    lst3 = lst3[select]
    mag3 = mag3[select]
    emag3 = emag3[select]
    
    plt.figure(figsize=(16,4))
    
    ax = plt.subplot(111)
    plt.title('ASCC {}'.format(ascc))
    plt.plot(jdmid1, mag1, '.')
    plt.plot(jdmid2, mag2, '.')
    plt.plot(jdmid3, mag3, '.')
    
    plt.xlabel('Time [JD]')
    
    plt.tight_layout()
    plt.show()
    
    return

def lc_diff(file1, file2, file3, ascc):
    
    with h5py.File(file1, 'r') as f:
        grp = f['data/' + ascc]
        jdmid1 = grp['jdmid'].value
        lst1 = grp['lst'].value
        mag1 = grp['mag0'].value
        emag1 = grp['emag0'].value
        nobs1 = grp['nobs'].value
        
    select = (nobs1 == 50)
    jdmid1 = jdmid1[select]
    lst1 = lst1[select]
    mag1 = mag1[select]
    emag1 = emag1[select]
        
    with h5py.File(file2, 'r') as f:
        grp = f['data/' + ascc]
        jdmid2 = grp['jdmid'].value
        lst2 = grp['lst'].value
        mag2 = grp['mag0'].value
        emag2 = grp['emag0'].value
        nobs2 = grp['nobs'].value
        
    select = (nobs2 == 50)
    jdmid2 = jdmid2[select]
    lst2 = lst2[select]
    mag2 = mag2[select]
    emag2 = emag2[select]
        
    with h5py.File(file3, 'r') as f:
        grp = f['data/' + ascc]
        jdmid3 = grp['jdmid'].value
        lst3 = grp['lst'].value
        mag3 = grp['mag0'].value
        emag3 = grp['emag0'].value
        nobs3 = grp['nobs'].value
        
    select = (nobs3 == 50)
    jdmid3 = jdmid3[select]
    lst3 = lst3[select]
    mag3 = mag3[select]
    emag3 = emag3[select]
    
    plt.figure(figsize=(16,9))
    
    ax = plt.subplot(311)
    plt.title('ASCC {}'.format(ascc))
    plt.plot(jdmid1, mag1 - np.nanmedian(mag1), '.')
    
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jdmid2, mag2 - np.nanmedian(mag2), '.')
    
    plt.subplot(313, sharex=ax, sharey=ax)
    plt.plot(jdmid3, mag3 - np.nanmedian(mag3), '.')
    
    plt.xlabel('Time [JD]')
    
    plt.tight_layout()
    plt.show()
    
    return
  
def as_array(filename, ascc, day0, ndays):
    
    with h5py.File(filename, 'r') as f:
        grp = f['data/' + ascc]
        lstseq = grp['lstseq'].value
        mag = grp['mag0'].value
        emag = grp['emag0'].value
        nobs = grp['nobs'].value
        
    select = (nobs == 50)
    lstseq = lstseq[select]
    mag = mag[select]
    emag = emag[select]
    
    lstday = lstseq//270
    lstidx = lstseq%270
    
    lstday = lstday - day0
    
    im = np.full((ndays, 270), fill_value=np.nan)
    im[lstday, lstidx] = mag
    
    return im, np.amin(lstidx), np.amax(lstidx)
    
def plot_array(file1, file2, file3, ascc):
    
    im1, xmin, xmax = as_array(file1, ascc, 827, 93)
    im2, xmin, xmax = as_array(file2, ascc, 827, 93)
    im3, xmin, xmax = as_array(file3, ascc, 827, 93)
    
    ax = plt.subplot(311)
    plt.imshow(im1, aspect='auto')
    
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.imshow(im2, aspect='auto')
    
    plt.subplot(313, sharex=ax, sharey=ax)
    plt.imshow(im3, aspect='auto')
    
    plt.xlim(xmin-.5, xmax+.5)
    
    plt.show()
    
    ax = plt.subplot(211)
    plt.imshow(im1-im2, aspect='auto')
    
    plt.subplot(212, sharex=ax, sharey=ax)
    plt.imshow(im1-im3, aspect='auto')
    
    plt.xlim(xmin-.5, xmax+.5)
    
    plt.show()
    
    return
    
def main():
    
    file1 = '/data2/talens/2015Q2/LPE/tmp/red0_2015Q2LPE.hdf5'
    file2 = '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5'
    file3 = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    
    #plot_array(file1, file2, file3, '803743')
    #plot_array(file1, file2, file3, '804411')
    #plot_array(file1, file2, file3, '806455')
    
    #plot_array(file1, file2, file3, '894573')
    #plot_array(file1, file2, file3, '714995')
    #plot_array(file1, file2, file3, '344250')
    
    lc_diff(file1, file2, file3, '803743')
    lc_diff(file1, file2, file3, '804411')
    lc_diff(file1, file2, file3, '806455')
    
    lc_diff(file1, file2, file3, '894573')
    lc_diff(file1, file2, file3, '714995')
    lc_diff(file1, file2, file3, '344250')
    
    file1 = '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5'
    file2 = '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5'
    file3 = '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5'
    
    lc_diff_mc(file1, file2, file3, '803743')
    lc_diff_mc(file1, file2, file3, '804411')
    lc_diff_mc(file1, file2, file3, '806455')
    
    lc_diff_mc(file1, file2, file3, '894573')
    lc_diff_mc(file1, file2, file3, '714995')
    lc_diff_mc(file1, file2, file3, '344250')
    
    file1 = '/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5'
    file2 = '/data2/talens/2015Q2/LPC/red0_vmag_2015Q2LPC.hdf5'
    file3 = '/data2/talens/2015Q2/LPW/red0_vmag_2015Q2LPW.hdf5'
    
    lc_diff_mc(file1, file2, file3, '803743')
    lc_diff_mc(file1, file2, file3, '804411')
    lc_diff_mc(file1, file2, file3, '806455')
    
    lc_diff_mc(file1, file2, file3, '894573')
    lc_diff_mc(file1, file2, file3, '714995')
    lc_diff_mc(file1, file2, file3, '344250')
    
    return

if __name__ == '__main__':
    main()
