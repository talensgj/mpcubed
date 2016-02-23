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

from scipy.optimize import curve_fit

from package.models import transit

import filters
from boxlstsq import boxlstsq

def read_data(filename, ascc):
    
    with h5py.File(filename, 'r') as f:
        
        try:
            grp = f['data/' + ascc]
        except:
            jdmid = None
            lst = None
            mag = None
            emag = None
        else:
        
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
        
    return jdmid, lst, mag, emag

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
          
    blsfile = '/data2/talens/2015Q2/bls0_2015Q2_patch266.hdf5'
    with h5py.File(blsfile, 'r') as f:
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
        
        jdmid = np.array([])
        lst = np.array([])
        mag = np.array([])
        emag = np.array([])
        
        for filename in data:
        
            jdmid_, lst_, mag_, emag_ = read_data(filename, ascc[i])
            
            if jdmid_ is None: continue
            
            mag_ = trend_filter(jdmid_, lst_, mag_, emag_)
        
            jdmid = np.append(jdmid, jdmid_)
            lst = np.append(lst, lst_)
            mag = np.append(mag, mag_)
            emag = np.append(emag, emag_)

        freq1, dchisq1, depth, hchisq, chisq0, epoch, duration = boxlstsq(jdmid, mag, 1/emag**2)
        
        epoch = epoch + duration/2
        depth = -depth
        
        arg = np.argmax(dchisq1)
        pin = [1/freq1[arg], jdmid[0] + epoch[arg], depth[arg], duration[arg]]
        
        model, mask = transit.box_model(jdmid, *pin)
        nt = np.sum(mask)*(320./(24*3600))/pin[3]
        
        phase = np.mod((jdmid - pin[1])/pin[0], 1)
        
        time = np.linspace(pin[1], pin[1] + pin[0], 500)
        model, mask = transit.box_model(time, *pin)
        mphase = np.mod((time - pin[1])/pin[0], 1)

        phase = np.mod(phase+.5, 1)-.5
        mphase = np.mod(mphase+.5, 1)-.5
        
        sort = np.argsort(mphase)
        mphase = mphase[sort]
        model = model[sort]
        
        fig = plt.figure(figsize=(16,9))
        
        plt.subplot(411)
        plt.title(r'ASCC {}, $P={:.2f}$ days, $N_t={:.1f}$'.format(ascc[i], pin[0], nt))
        plt.errorbar(phase, mag, yerr=emag, fmt='.', c='k')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        
        
        pin.append(10.)
        
        pin[1] = pin[1] - np.amin(jdmid)
        jdmid = jdmid - np.amin(jdmid)
        popt, pcov = curve_fit(transit.softened_box_model, jdmid, mag, pin, emag, absolute_sigma=True)
        
        print pin
        print popt
        
        phase = np.mod((jdmid - popt[1])/popt[0], 1)
        
        time = np.linspace(popt[1], popt[1] + popt[0], 500)
        model = transit.softened_box_model(time, *popt)
        mphase = np.mod((time - popt[1])/popt[0], 1)

        phase = np.mod(phase+.5, 1)-.5
        mphase = np.mod(mphase+.5, 1)-.5
        
        sort = np.argsort(mphase)
        mphase = mphase[sort]
        model = model[sort]
        
        plt.subplot(412)
        plt.errorbar(phase, mag, yerr=emag, fmt='.', c='k')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        
        nbins = 9/(popt[3]/popt[0])
        edges = np.linspace(-.5, .5, nbins+1)
        x0, edges = np.histogram(phase, bins=edges, weights=1/emag**2)
        x1, edges = np.histogram(phase, bins=edges, weights=mag/emag**2)
        
        bmag = x1/x0
        bemag = np.sqrt(1./x0)
        bphase = (edges[1:] + edges[:-1])/2
        
        plt.subplot(413)
        plt.errorbar(bphase, bmag, yerr=bemag, fmt='.', c='k')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        
        plt.subplot(414)
        plt.plot(freq, dchisq[:,i], c='k')
        plt.axvline(1/period[i], c='r', ls='--')
        plt.xlabel(r'Frequency [day$^{-1}$]')
        plt.ylabel(r'$\Delta\chi^2$')
        
        plt.savefig('candidate_patch266_ASCC{}.png'.format(ascc[i]))
        plt.show()

    return

if __name__ == '__main__':
    main()
