#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
import time

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

from package.models import transit
from package.coordinates import grids

import filters
from boxlstsq import boxlstsq
from boxlstsq_ms import boxlstsq_ms

# Set random seed.
np.random.seed(19910909)

def data_as_array(filename, ascc):
    
    lstseq = np.array([])
    staridx = np.array([])
    jdmid = np.array([])
    lst = np.array([])
    mag0 = np.array([])
    emag0 = np.array([])
    
    with h5py.File(filename, 'r') as f:
        
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
    jdmid = jdmid[args]
    lst = lst[args]
    
    tmp = np.zeros((len(ascc), len(lstseq)))
    tmp[staridx, idx] = mag0
    mag0 = tmp
    
    tmp = np.zeros((len(ascc), len(lstseq)))
    tmp[staridx, idx] = emag0
    emag0 = tmp
    
    tmp = np.zeros((len(ascc), len(lstseq)), dtype='bool')
    tmp[staridx, idx] = True
    mask = ~tmp
    
    return ascc, jdmid, lst, mag0, emag0, mask

def inject_signals(jdmid, mag0):
    
    nstars = mag0.shape[0]
    
    P = 1 + (np.ptp(jdmid)/9 - 1)*np.random.rand(nstars)
    Tp = jdmid[0] + P*np.random.rand(nstars)
    delta = .005 + (.03 - .005)*np.random.rand(nstars)
    eta = (1.8/24)*P**(1./3)
    
    for i in range(nstars):
        
        signal = transit.softened_box_model(jdmid, P[i], Tp[i], delta[i], eta[i])
        mag0[i] = mag0[i] - signal # add or subtract?
        
    return P, mag0

def trend_filter(jdmid, lst, mag0, emag0):
    
    nstars = mag0.shape[0]
    weights = np.where(emag0 > 1e-3, 1./emag0**2, 0.)
    
    for i in range(nstars):
        #chisq, pars, fit = filters.harmonic(lst, mag0[i], weights[i], 24., 8)
        chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag0[i], weights[i], 180., 20)
        mag0[i] = mag0[i] - fit
    
    return mag0
    
def test_on_skybin():
    
    filename = '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5'
    
    with h5py.File(filename, 'r') as f:
        grp = f['header']
        ascc = grp['ascc'].value
        ra = grp['ra'].value
        dec = grp['dec'].value
        
    hg = grids.HealpixGrid(8)
    skyidx = hg.radec2idx(ra, dec)
    select = (skyidx == 266)
    ascc = ascc[select]
    
    ascc, jdmid, lst, mag0, emag0, mask = data_as_array(filename, ascc)
    
    plt.figure(figsize=(16,9))
    plt.imshow(mag0, aspect='auto', interpolation='None', vmin=-.1, vmax=.1, cmap=cm.Greys)
    cb = plt.colorbar()
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    plt.xlabel('Time')
    plt.tight_layout()
    plt.show()
    
    P, mag0 = inject_signals(jdmid, mag0)
    
    plt.figure(figsize=(16,9))
    plt.imshow(mag0, aspect='auto', interpolation='None', vmin=-.1, vmax=.1, cmap=cm.Greys)
    cb = plt.colorbar()
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    plt.xlabel('Time')
    plt.tight_layout()
    plt.show()
    
    mag0 = trend_filter(jdmid, lst, mag0, emag0)
    
    plt.figure(figsize=(16,9))
    plt.imshow(mag0, aspect='auto', interpolation='None', vmin=-.1, vmax=.1, cmap=cm.Greys)
    cb = plt.colorbar()
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    plt.xlabel('Time')
    plt.tight_layout()
    plt.show()
    
    for i in range(len(ascc)):
        if ascc[i] != '891113': continue
        
        phase = np.mod(jdmid/P[i], 1)
        phase = np.mod(phase + .5, 1) - .5 
        
        plt.figure(figsize=(16,5))
        plt.title('ASCC {}'.format(ascc[i]))
        plt.errorbar(phase, mag0[i], yerr=emag0[i], fmt='.')
        plt.xlim(-.5,.5)
        plt.ylim(.1,-.1)
        plt.xlabel('Phase')
        plt.ylabel('Magnitude')
        plt.tight_layout()
        plt.show()
        plt.close()
        
    return
    
def main():
    
    filename = '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5'
    
    with h5py.File(filename, 'r') as f:
        grp = f['header']
        ascc = grp['ascc'].value
    
    ascc = np.random.choice(ascc, 500, replace=False)
    
    ascc, jdmid, lst, mag0, emag0, mask = data_as_array(filename, ascc)
    weights = np.where(mask, 0, 1/emag0**2)
        
    plt.figure(figsize=(16,9))
    plt.imshow(mag0, aspect='auto', interpolation='None', vmin=-.1, vmax=.1, cmap=cm.Greys)
    cb = plt.colorbar()
    cb.ax.invert_yaxis()
    cb.set_label('Magnitude')
    plt.xlabel('Time')
    plt.tight_layout()
    plt.show()
    
    start = time.time()
    #boxlstsq_ms(jdmid, mag0.T, weights.T)
    for i in range(500):
        boxlstsq(jdmid, mag0[i], weights[i])
    print time.time() - start
        
    return 0

if __name__ == '__main__':
    main()
