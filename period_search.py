#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import h5py
import numpy as np

import time

import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package import IO
from package.coordinates import grids

from boxlstsq import boxlstsq
from boxlstsq_ms import boxlstsq_ms

def period_search():

    filename = '/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5'
    
    with h5py.File(filename, 'r') as f:
        ascc = f['header/ascc'].value
        ra = f['header/ra'].value
        dec = f['header/dec'].value

    hg = grids.HealpixGrid(8)
    idx = hg.radec2idx(ra, dec)
    
    arg, = np.where(ascc == '807144')
    select = (idx == idx[arg])
    ascc = ascc[select]
    
    lstseq = np.array([], dtype='int')
    lst = np.array([])
    jdmid = np.array([])
    mag0 = np.array([])
    emag0 = np.array([])
    nobs = np.array([])
    staridx = np.array([], dtype='int')

    with h5py.File(filename, 'r') as f:
        for i in range(len(ascc)):
            grp = f['data/' + ascc[i]]
            lstseq = np.append(lstseq, grp['lstseq'])
            staridx = np.append(staridx, [i]*len(grp['lstseq']))
            lst = np.append(lst, grp['lst'])
            jdmid = np.append(jdmid, grp['jdmid'])
            mag0 = np.append(mag0, grp['mag0'])
            emag0 = np.append(emag0, grp['emag0'])
            nobs = np.append(nobs, grp['nobs'])
            
    emag0 = emag0/np.sqrt(nobs)
    
    select = (nobs == 50)
    lstseq = lstseq[select]
    staridx = staridx[select]
    lst = lst[select]
    jdmid = jdmid[select]
    mag0 = mag0[select]
    emag0 = emag0[select]
    nobs = nobs[select]
    
    lstseq, args, idx = np.unique(lstseq, return_index=True, return_inverse=True)
    
    npoints = len(lstseq)
    nstars = len(ascc)
    
    jdmid = jdmid[args]
    
    tmp = np.full((nstars, npoints), fill_value=0)
    tmp[staridx,idx] = mag0
    mag0 = tmp
    
    tmp = np.full((nstars, npoints), fill_value=0)
    tmp[staridx,idx] = 1/emag0**2
    weights = tmp
    
    jdmid = np.append(jdmid, jdmid)
    mag0 = np.append(mag0, mag0, axis=1)
    weights = np.append(weights, weights, axis=1)
    
    plt.imshow(mag0, aspect='auto')
    plt.xlabel('Time')
    cb = plt.colorbar()
    cb.set_label('Magnitude')
    plt.show()
    
    #start = time.time()
    #freq, dchisq, depth, hchisq = boxlstsq_ms(jdmid, mag0.T, weights.T)
    #print time.time() - start
    
    start = time.time()
    for i in range(len(ascc)):
        freq_0, dchisq_0, depth_0, hchisq_0 = boxlstsq(jdmid, mag0[i], weights[i])
        
        #crit = np.allclose(dchisq[:,i], dchisq_0) & np.allclose(depth[:,i], depth_0)
        #crit = crit & np.allclose(hchisq[:,i], hchisq_0)
        
        #if not crit:
            #print 'Difference found', i
        
        #plt.plot(freq, dchisq[:,i])
        #plt.plot(freq_0, dchisq_0)
        #plt.xlabel(r'Frequency [day$^{-1}$]')
        #plt.ylabel(r'$\Delta\chi^2$')
        #plt.show()
        #plt.close()
    print time.time() - start
    
    return
    
if __name__ == '__main__':
    period_search()
