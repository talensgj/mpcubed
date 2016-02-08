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
    
data = np.zeros((len(ascc), 4))
n1 = 0
n2 = 0
with h5py.File(boxlstsq, 'r') as f:
    
    for i in range(len(ascc)):
        
        flag = 0
        
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
        
        with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as g:
            grp = g['data/' + ascc[i]]
            jdmid = grp['jdmid'].value
    
        # Best fit.
        arg = np.argmax(dchisq)
        best_chisq = chisq0 - dchisq[arg]
        best_freq = freq[arg]
        best_depth = depth[arg]
            
        # Good fit?
        quality = best_chisq/nobs
        data[i,0] = quality
        
        if (quality > 3.5):
            flag += 1
            
        # Anti-transit ratio.
        args1, = np.where(depth > 0)
        args2, = np.where(depth < 0)
        ratio = np.max(dchisq[args1])/np.max(dchisq[args2])
        data[i,1] = ratio
        
        if (ratio < 1.5):
            flag += 2
        
        # Data gaps.
        q = (1.8/24)*best_freq**(2./3) # Fractional transit duration.
        phase = np.mod(jdmid*best_freq, 1)
        phase = np.sort(phase)
        gapsizes = np.diff(phase)
        gapsizes = np.append(gapsizes, 1-np.ptp(phase))
        gapsizes = np.amax(gapsizes)/q
        data[i,2] = gapsizes
        
        if (gapsizes > 2.5):
            flag += 4
        
        #
        if (np.abs(P[i]*best_freq - 1) < .05):
            n1 += 1
            
            if (flag == 0):
                n2 +=1
    
        data[i,3] = flag
    
        # Plot the point.
        plt.scatter(P[i], 1/best_freq, c=flag, vmin=0, vmax=1)

print n1, n2

plt.xlabel(r'$P$ [days]')
plt.ylabel(r'$P_{rec}$ [days]')
plt.show()
    
#with h5py.File(boxlstsq) as f:
    #f.create_dataset('flag', data=data[:,3])
    
plt.hist(data[:,0], bins=np.linspace(0,20,41))
plt.axvline(3.5, c='k')
plt.xlabel(r'$\chi^2/N$')
plt.show()
    
plt.hist(data[:,1], bins = np.linspace(0, 10, 21))
plt.axvline(1.5, c='k')
plt.xlabel(r'$\Delta\chi^2/\Delta\chi^2_{-}$')
plt.show()
    
plt.hist(data[:,2], bins = np.linspace(0,10,21))
plt.axvline(2.5, c='k')
plt.xlabel('$\Delta \phi/q$')
plt.show()
    
    
    
