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

import filters

def trend_filter(jdmid, lst, mag, emag):
    
    weights = 1./emag**2
    
    #chisq, pars, fit = filters.harmonic(lst, mag, weights, 24., 8)
    chisq, pars, fit = filters.masc_harmonic(jdmid, lst, mag, weights, 180., 20)
    mag = mag - fit
    
    return mag

def main():
    
    data = ['/data2/talens/2015Q2/LPN/red0_2015Q2LPN.hdf5',
            '/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2/LPS/red0_2015Q2LPS.hdf5',
            '/data2/talens/2015Q2/LPW/red0_2015Q2LPW.hdf5',
            '/data2/talens/2015Q2/LPC/red0_2015Q2LPC.hdf5']
          
    boxlstsq = '/data2/talens/2015Q2/bls0_2015Q2_patch266.hdf5'
    with h5py.File(boxlstsq, 'r') as f:
        grp = f['header']
        ascc = grp['ascc'].value
        flag = grp['flag'].value
        period = grp['period'].value
        
        grp = f['data']
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
            
    select = (flag == 0)
    ascc = ascc[select]
    period = period[select]
    dchisq = dchisq[:,select]
    
    for i in range(len(ascc)):
        
        fig = plt.figure(figsize=(16,9))
        
        plt.subplot(211)
        plt.title(r'ASCC {}, $P={:.2f}$ days'.format(ascc[i], period[i]))
        
        for filename in data:
        
            with h5py.File(filename, 'r') as f:
                try:
                    grp = f['data/' + ascc[i]]
                except:
                    continue
                    
                jdmid = grp['jdmid'].value
                lst = grp['lst'].value
                mag = grp['mag0'].value
                emag = grp['emag0'].value
                nobs = grp['nobs'].value
                
                emag = emag/np.sqrt(nobs)
                
                select = (nobs == 50)
                jdmid = jdmid[select]
                lst = lst[select]
                mag = mag[select]
                emag = emag[select]
                
                mag = trend_filter(jdmid, lst, mag, emag)
        
                phase = np.mod(jdmid/period[i], 1.)
        
                #q = (1.8/24)*period[i]**(-2./3)
                #nbins = np.ceil(9/q)
                #edges = np.linspace(0, 1, nbins+1)
        
                #npbin, edges = np.histogram(phase, bins=edges)
                #x1, edges = np.histogram(phase, bins=edges, weights=mag)
                #x2, edges = np.histogram(phase, bins=edges, weights=mag**2)
                #mag = x1/npbin
                #emag = np.sqrt(x2/npbin - (x1/npbin)**2)/np.sqrt(npbin)
                #phase = (edges[1:] + edges[:-1])/2
        
                #x0, edges = np.histogram(phase, bins=np.linspace(0,1,251), weights=1/emag**2)
                #x1, edges = np.histogram(phase, bins=np.linspace(0,1,251), weights=mag/emag**2)
                
                #mag = x1/x0
                #emag = np.sqrt(1./x0)
                #phase = (edges[1:] + edges[:-1])/2
        
                plt.errorbar(phase, mag, yerr=emag, fmt='.')
                
        plt.xlim(0, 1)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel('Magnitude')
        
        plt.subplot(212)
        plt.plot(freq, dchisq[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.xlabel(r'Frequency [day$^{-1}$]')
        plt.ylabel(r'$\Delta\chi^2$')
        
        plt.savefig('candidate_patch266_ASCC{}.png'.format(ascc[i]))
        plt.show()
        plt.close()
        
    return

if __name__ == '__main__':
    main()
