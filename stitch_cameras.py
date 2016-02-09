#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package.coordinates import grids

import filters


def read_skybin(filename):
    
    with h5py.File(filename, 'r') as f:
        
        # Read the header.
        grp = f['header']
        ascc = grp['ascc'].value
        ra = grp['ra'].value
        dec = grp['dec'].value
        
        # Select the stars in the skybin.
        hg = grids.HealpixGrid(8)
        skyidx = hg.radec2idx(ra, dec)
        select = (skyidx == 266)
        ascc = ascc[select]
        
        lstseq = np.array([])
        staridx = np.array([])
        jdmid = np.array([])
        lst = np.array([])
        mag0 = np.array([])
        emag0 = np.array([])
        
        for i in range(len(ascc)):
            
            grp = f['data/' + ascc[i]]
            lstseq_ = grp['lstseq'].value
            jdmid_ = grp['jdmid'].value
            lst_ = grp['lst'].value
            mag0_ = grp['mag0'].value
            emag0_ = grp['emag0'].value
            nobs_ = grp['nobs'].value
            
            emag0_ = emag0_/np.sqrt(nobs_)
            
            select = (nobs_ == 50)
            lstseq_ = lstseq_[select]
            jdmid_ = jdmid_[select]
            lst_ = lst_[select]
            mag0_ = mag0_[select]
            emag0_ = emag0_[select]
            
            lstseq = np.append(lstseq, lstseq_)
            staridx = np.append(staridx, [i]*len(lstseq_))
            jdmid = np.append(jdmid, jdmid_)
            lst = np.append(lst, lst_)
            mag0 = np.append(mag0, mag0_)
            emag0 = np.append(emag0, emag0_)
          
    staridx = staridx.astype('int')
            
    lstseq, args, idx = np.unique(lstseq, return_index=True, return_inverse=True)
    
    nstars = len(ascc)
    npoints = len(lstseq)
    
    jdmid = jdmid[args]
    lst = lst[args]
    
    tmp = np.zeros((nstars, npoints))
    tmp[staridx, idx] = mag0
    mag0 = tmp
    
    tmp = np.zeros((nstars, npoints))
    tmp[staridx, idx] = 1/emag0**2
    weights = tmp
        
    return ascc, jdmid, lst, mag0, weights
        
def main():
    
    data = ['/data2/talens/2015Q2/LPN/red0_2015Q2LPN.hdf5',
            '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPS/red0_2015Q2LPS.hdf5',
            '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5',
            '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5']
    
    ascc = '807144'
    
    for filename in data:
    
        with h5py.File(filename, 'r') as f:
            grp = f['data/' + ascc[i]]
            lst = grp['lst'].value
            jdmid = grp['jdmid'].value
            mag0 = grp['mag0'].value
            emag0 = grp['emag0'].value
            nobs = grp['nobs'].value
            
        emag0 = emag0/np.sqrt(nobs)
        
        select = (nobs == 50)
        lst = lst[select]
        jdmid = jdmid[select]
        mag0 = mag0[select]
        emag0 = emag0[select]
        
        weights = 1/emag0**2
        base = np.ptpt(jdmid)
        chisq, pars, fit = filters.harmonic(jdmid, lst, mag0, weights, 2*base, 5) 
        
        plt.errorbar(jdmid, mag0 - fit, yerr=emag0)
    plt.show()
    
    return
        
if __name__ == '__main__':
    main()
        
        
            
            
            
            
            
