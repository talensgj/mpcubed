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

from package.coordinates import grids

import filters

def read_stars(filename, ascc):
    
    with h5py.File(filename, 'r') as f:
    
        lstseq = np.array([])
        staridx = np.array([])
        jdmid = np.array([])
        lst = np.array([])
        mag0 = np.array([])
        emag0 = np.array([])
        
        for i in range(len(ascc)):
            
            try:
                grp = f['data/' + ascc[i]]
            except:
                continue
                
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
            
            try:
                grp = f['data/' + ascc[i]]
            except:
                continue
                
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
     
def data_as_arrays():
    
    data = ['/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5',
            '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5']
    
    with h5py.File(data[0], 'r') as f:
        grp = f['header']
        ascc = grp['ascc'].value
        ra = grp['ra'].value
        dec = grp['dec'].value
        
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)
    select = (skyidx == 266)
    ascc = ascc[select]
    
    for filename in data:
        ascc, jdmid, lst, mag0, weights = read_stars(filename, ascc)
        
        fit = np.zeros(mag0.shape)
        for i in range(len(ascc)):
            chisq, pars, fit[i] = filters.harmonic(lst, mag0[i], weights[i], 24., 8)
            ##chisq, pars, fit = filters.harmonic(jdmid, mag0, weights, 180., 20) 
            ##chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag0, weights, 180., 5)
        
        plt.subplot(311)
        plt.imshow(mag0, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
        plt.colorbar()
        
        plt.subplot(312)
        plt.imshow(fit, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
        plt.colorbar()
        
        plt.subplot(313)
        plt.imshow(mag0 - fit, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
        plt.colorbar()
        
        plt.show()
        
    return
    
def data_per_star(ascc):
    
    data = ['/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPC/red0_vmag_2015Q2LPC.hdf5',
            '/data2/talens/2015Q2/LPW/red0_vmag_2015Q2LPW.hdf5']      

    lstseq = np.array([])
    camidx = np.array([])
    mag0 = np.array([])
    emag0 = np.array([])
    fit = np.array([])

    for i in range(3):
        
        with h5py.File(data[i], 'r') as f:
            
            try:
                grp = f['data/' + ascc]
            except:
                continue
                
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
        
        if len(jdmid_) == 0: continue
        
        weights = 1/emag0_**2
        
        #chisq, pars, fit_ = filters.harmonic(lst_, mag0_, weights, 24., 8)
        #chisq, pars, fit = filters.harmonic(jdmid, mag0, weights, 180., 20) 
        chisq, pars, fit_ = filters.masc_harmonic(jdmid_, lst_, mag0_, weights, 180., 20)
        
        lstseq = np.append(lstseq, lstseq_)
        camidx = np.append(camidx, [i]*len(lstseq_))
        mag0 = np.append(mag0, mag0_)
        emag0 = np.append(emag0, emag0_)
        fit = np.append(fit, fit_)

    lstseq, idx = np.unique(lstseq, return_inverse=True)
    camidx = camidx.astype('int')
        
    tmp = np.full((5, len(lstseq)), fill_value=np.nan)
    tmp[camidx, idx] = mag0
    mag0 = tmp
       
    tmp = np.full((5, len(lstseq)), fill_value=np.nan)
    tmp[camidx, idx] = emag0
    emag0 = tmp
    
    tmp = np.full((5, len(lstseq)), fill_value=np.nan)
    tmp[camidx, idx] = fit
    fit = tmp
       
    plt.figure(figsize=(16,9))
        
    ax = plt.subplot(311)
    plt.title('ASCC {}'.format(ascc))
    plt.plot(np.arange(len(lstseq)), mag0.T, '.')
    plt.ylabel('Magnitude')

    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(np.arange(len(lstseq)), fit.T, '.')
    plt.ylabel('Magnitude')
    plt.ylim(.2, -.2)

    plt.subplot(313, sharex=ax)
    plt.plot(np.arange(len(lstseq)), mag0.T - fit.T, '.')
    plt.ylim(.1, -.1)
    plt.xlabel('Time')
    plt.ylabel('Magnitude')
        
    plt.tight_layout()
    plt.savefig('ASCC{}_vmag.png'.format(ascc))
    #plt.show()
    plt.close()
    
    return
    
def on_N_cameras():
    
    data = ['/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5',
            '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5']    

    ascc = np.array([])
    ra = np.array([])
    dec = np.array([])
    
    for filename in data:
        
        with h5py.File(filename, 'r') as f:
            grp = f['header']
            ascc_ = grp['ascc'].value
            ra_ = grp['ra'].value
            dec_ = grp['dec'].value
            
            ascc = np.append(ascc, ascc_)
            ra = np.append(ra, ra_)
            dec = np.append(dec, dec_)
            
    ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
    ra = ra[args]
    dec = dec[args]
    ncams = np.bincount(idx)    
    
    return ascc, ra, dec, ncams
        
def main():
    
    ascc, ra, dec, ncams = on_N_cameras()
    
    plt.subplot(111, projection='mollweide')
    plt.scatter(ra*np.pi/180-np.pi, dec*np.pi/180, c=ncams, vmin=1, vmax=5)
    plt.colorbar()
    plt.show()
    
    select = (ncams > 1)
    ascc = ascc[select]
    np.random.seed(19910909)
    ascc = np.random.choice(ascc, 100, replace=False)
    
    for i in range(len(ascc)):
        data_per_star(ascc[i])
    
    return
        
if __name__ == '__main__':
    main()
        
        
            
            
            
            
            
