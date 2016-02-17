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
    
    ascc, jdmid, lst, mag0, weights = read_skybin('/data2/talens/2015Q2/LPE/red0_vmag_2015Q2LPE.hdf5')
    
    base = np.ptp(jdmid)
    for i in range(mag0.shape[0]):
        #if ascc[i] != '807144': continue
        
        chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag0[i], weights[i], 180., 20) 
        
        plt.figure(figsize=(16,9))
        
        ax = plt.subplot(311)
        ax.invert_yaxis()
        plt.title(r'ASCC {}, $\chi^2 = {:.2f}$'.format(ascc[i], chisq))
        plt.plot(mag0[i], '.')
        plt.ylim(-.1, .1)
        
        plt.subplot(312, sharex=ax, sharey=ax)
        plt.plot(fit, '.')
        plt.ylim(-.1, .1)
    
        plt.subplot(313, sharex=ax, sharey=ax)
        plt.plot(mag0[i] - fit, '.')
        plt.ylim(.1, -.1)
    
        plt.savefig('ASCC{}_vmag.png'.format(ascc[i]))
        #plt.show()
        plt.close()
    
    
    #avg = np.nansum(weights*mag0, axis=1, keepdims=True)/np.nansum(weights, axis=1, keepdims=True)
    #mag0 = mag0 - avg
    
    #plt.imshow(mag0, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
    #plt.show()
    
    #fit = np.zeros(mag0.shape)
    #for i in range(5):
        #tmp = filters.sysrem(mag0 - fit, weights)
        #fit = fit + tmp
    
    #plt.imshow(fit, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
    #plt.show()
    
    #for i in range(mag0.shape[0]):
        ##if ascc[i] != '807144': continue
        
        #plt.figure(figsize=(16,9))
        
        #ax = plt.subplot(311)
        #ax.invert_yaxis()
        #plt.title('ASCC {}'.format(ascc[i]))
        #plt.plot(mag0[i], '.')
        #plt.ylim(-.1, .1)
        
        #plt.subplot(312, sharex=ax, sharey=ax)
        #plt.plot(fit[i], '.')
        #plt.ylim(-.1, .1)
    
        #plt.subplot(313, sharex=ax, sharey=ax)
        #plt.plot(mag0[i] - fit[i], '.')
        #plt.ylim(.1, -.1)
    
        #plt.savefig('ASCC{}.png'.format(ascc[i]))
        ##plt.show()
        #plt.close()
    
    return
        
if __name__ == '__main__':
    main()
        
        
            
            
            
            
            
