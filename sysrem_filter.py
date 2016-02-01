#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import h5py
import numpy as np

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

import fourierfuncs

def sysrem(mag0, emag0):
    
    weights = 1/emag0**2
    par2 = np.ones((mag0.shape[0], 1))
    for i in range(5):
        
        par1 = np.nansum(weights*mag0*par2, axis=1, keepdims=True)/np.nansum(weights*par2**2, axis=1, keepdims=True)
        par2 = np.nansum(weights*mag0*par1, axis=0, keepdims=True)/np.nansum(weights*par1**2, axis=0, keepdims=True)
    
    fit = np.outer(par1, par2)
    
    return fit
    
def harmonics_filter(jdmid, lst, mag0, emag0, Pjd, njd, Plst=24., nlst=5):
    
    weights = 1/emag0**2
    
    fmat = np.ones((len(mag0),1))
    
    if (njd > 0):
        mat_jd = fourierfuncs.fourier_mat(jdmid, 1/Pjd, njd)
        fmat = np.hstack([fmat, mat_jd])
        
    if (nlst > 0):
        mat_lst = fourierfuncs.fourier_mat(lst, 1/Plst, nlst)
        fmat = np.hstack([fmat, mat_lst])
    
    pars = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), mag0*np.sqrt(weights))[0]
    fit = np.dot(fmat, pars)
    
    chisq = np.sum(weights*(mag0 - fit)**2)/(len(mag0) - len(pars))
    
    return chisq, pars, fit

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
    
    for i in range(len(ascc)):
        if ascc[i] != '891488': continue
        
        select = (staridx == i)
        
        BL = 2*np.ptp(jdmid[select])
        chisq, pars, fit = harmonics_filter(jdmid[select], lst[select], mag0[select], emag0[select], 2*BL, 5)
    
        plt.figure(figsize=(16,9))
        
        ax = plt.subplot(311)
        ax.invert_yaxis()
        plt.title('ASCC {}'.format(ascc[i]))
        plt.plot(mag0[select], '.')
        plt.ylabel('Magnitude')
        
        plt.subplot(312, sharex=ax, sharey=ax)
        plt.plot(fit, '.')
        plt.ylabel('Magnitude')
        
        plt.subplot(313, sharex=ax)
        plt.plot(mag0[select] - fit, '.')
        plt.ylim(.1, -.1)
        plt.xlabel('Time')
        plt.ylabel('Magnitude')
        
        plt.show()
        plt.close()
    
    #lstseq, idx = np.unique(lstseq, return_inverse=True)
    
    #npoints = len(lstseq)
    #nstars = len(ascc)
    
    #tmp = np.full((nstars, npoints), fill_value=np.nan)
    #tmp[staridx,idx] = mag0
    #mag0 = tmp
    
    #tmp = np.full((nstars, npoints), fill_value=np.nan)
    #tmp[staridx,idx] = emag0
    #emag0 = tmp
    
    #m = np.nanmean(mag0/emag0**2, axis=1, keepdims=True)/np.nanmean(1/emag0**2, axis=1, keepdims=True)
    #mag0 = mag0 - m
        
    #fit = 0
    #for i in range(5):
        
        #comp = sysrem(mag0 - fit, emag0)
        
        #plt.imshow(comp, aspect='auto', interpolation='None')
        #plt.show()
        
        #fit = fit + comp
    
    #plt.subplot(311)
    #plt.imshow(mag0, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
    
    #plt.subplot(312)    
    #plt.imshow(fit, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
    
    #plt.subplot(313)
    #plt.imshow(mag0 - fit, aspect='auto', interpolation='None', vmin=-.1, vmax=.1)
    
    #plt.show()
    
    #i, = np.where(ascc == '807144')
    
    #for i in range(len(ascc)):
        #plt.figure(figsize=(16,9))
        
        #ax = plt.subplot(311)
        #ax.invert_yaxis()
        #plt.title('ASCC {}'.format(ascc[i]))
        #plt.plot(mag0[i].T, '.')
        #plt.ylabel('Magnitude')
        
        #plt.subplot(312, sharex=ax, sharey=ax)
        #plt.plot(fit[i].T, '.')
        #plt.ylabel('Magnitude')
        
        #plt.subplot(313, sharex=ax, sharey=ax)
        #plt.plot((mag0 - fit)[i].T, '.')
        #plt.ylim(.1, -.1)
        #plt.xlabel('Time')
        #plt.ylabel('Magnitude')
        
        #plt.show()
        #plt.close()
    
    return
    
if __name__ == '__main__':
    period_search()
