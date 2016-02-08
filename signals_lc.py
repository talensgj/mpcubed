#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

import filters

## Constants ##
G = 6.67384e-11 # [m^3 kg^-1 s^-2]
SecInDay = 60*60*24 # number of seconds in one day
Msun = 1.9891e30 # [kg]
Rsun = 696342e3 # [m] 

index = '/data2/talens/inj_signals/signals/signals_index.hdf5'
boxlstsq = '/data2/talens/inj_signals/signals/signals_boxlstsq_minusm.hdf5'

with h5py.File(index, 'r') as f:
    ascc = f['ascc'].value
    P = f['P'].value
    
with h5py.File(boxlstsq, 'r') as f:
    
    flag = f['flag'].value
    flag = flag.astype('int')
    
    for i in range(len(ascc)):
        
        try:
            grp = f[ascc[i]]
        except:
            continue
        else:
            freq = grp['freq'].value
            dchisq = grp['dchisq'].value
            depth = grp['depth'].value
            nobs = grp.attrs['nobs']
            chisq0 = grp.attrs['chisq0']
            chisq = grp.attrs['chisq']
        
        arg = np.argmax(dchisq)
        best_chisq = chisq0 - dchisq[arg]
        best_freq = freq[arg]
        best_depth = depth[arg]
        
        print chisq0, dchisq[arg], nobs
        continue
        with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as g:
            grp = g['data/' + ascc[i]]
            jdmid = grp['jdmid'].value
            lst = grp['lst'].value
            mag0 = grp['mag0'].value
            emag0 = grp['emag0'].value
            nobs = grp['nobs'].value
        
        emag0 = emag0/np.sqrt(nobs)
            
        select = (nobs == 50)
        jdmid = jdmid[select]
        lst = lst[select]
        mag0 = mag0[select]
        emag0 = emag0[select]
        
        weights = 1/emag0**2
        
        base = np.ptp(jdmid)
        chisq, pars, fit = filters.harmonic(jdmid, lst, mag0, weights, 2*base, 5)
        mag0 = mag0 - fit
        
        phase = jdmid*best_freq
        phase = np.mod(phase, 1)
        
        plt.figure(figsize=(16,5))
        ax = plt.subplot(111)
        ax.invert_yaxis()
        plt.title('ASCC {}, $P={:.2f}$, $P_r={:.2f}$, flag = {}'.format(ascc[i], P[i], 1/best_freq, flag[i]))
        plt.errorbar(phase, mag0, yerr=emag0, fmt='.')
        plt.ylim(.1, -.1)
        plt.xlim(0, 1)
        plt.xlabel('Phase')
        plt.ylabel('Magnitude')
        plt.tight_layout()
        #plt.show()
        plt.savefig('/data2/talens/inj_signals/signals/lightcurves/flag{}_ASCC{}.png'.format(flag[i], ascc[i]))
        plt.close()
    
        
