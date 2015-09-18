#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from coordinate_grids import HealpixGrid
from index_functions import index_statistics
from BLS_ms import BLS

period = 2.21857312      						#planet orbital period (days)
t0 = 2454037.612 

filelist = glob.glob('/data2/talens/Jul2015/fLC_201507??LPE.hdf5')

hg = HealpixGrid(16)

with h5py.File(filelist[0]) as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value

hgidx = hg.find_gridpoint(ra, dec)
print hgidx[ascc=='807144']
here = hgidx == 981
ascc = ascc[here]
print ascc
arg, = np.where(ascc=='807144')
print arg, ascc[arg]

jdmid = np.array([])
lstidx = np.array([])
flux = np.array([])
eflux = np.array([])
flags = np.array([])
staridx = np.array([])

for i in range(len(filelist)):

    for j in range(len(ascc)):

        with h5py.File(filelist[i]) as f:
            
            try:
                lc = f['data/'+ascc[j]]
            except:
                continue
            
            jdmid = np.append(jdmid, lc['jdmid'])
            lstidx = np.append(lstidx, lc['lstidx'])
            flux = np.append(flux, lc['flux0'])
            eflux = np.append(eflux, lc['eflux0'])
            flags = np.append(flags, lc['flag'])
            staridx = np.append(staridx, [j]*len(lc['jdmid']))
            
weights = 1./eflux**2
weights[flags>0] = 0.
weights[flux<0] = 0.
weights[eflux<0] = 0.

dayidx = np.floor(jdmid).astype('int')
dayidx = dayidx - np.amin(dayidx)
lstidx = lstidx.astype('int')
staridx = staridx.astype('int')

Fidx = np.ravel_multi_index((staridx, dayidx), (np.amax(staridx)+1, np.amax(dayidx)+1))
Tidx = np.ravel_multi_index((staridx, lstidx), (np.amax(staridx)+1, 13500))
Sidx = np.ravel_multi_index((dayidx, lstidx), (np.amax(dayidx)+1, 13500))

Flen = np.amax(Fidx)+1
Tlen = np.amax(Tidx)+1
Slen = np.amax(Sidx)+1

#T = np.ones(Tlen)
#S = np.ones(Slen)

#for i in range(10):
    
    #F = np.bincount(Fidx, weights*flux*(T[Tidx]*S[Sidx]), minlength=Flen)/np.bincount(Fidx, weights*(T[Tidx]*S[Sidx])**2, minlength=Flen) # F    
    #T = np.bincount(Tidx, weights*flux*(F[Fidx]*S[Sidx]), minlength=Tlen)/np.bincount(Tidx, weights*(F[Fidx]*S[Sidx])**2, minlength=Tlen) # F  
    #S = np.bincount(Sidx, weights*flux*(F[Fidx]*T[Tidx]), minlength=Slen)/np.bincount(Sidx, weights*(F[Fidx]*T[Tidx])**2, minlength=Slen) # F  
    
    #F[np.isnan(F)] = 0.
    #T[np.isnan(T)] = 0.
    #S[np.isnan(S)] = 0.
    
    #if i > 0:
        #print np.nanmax(np.abs((F-F_old)/F_old)), np.nanmax(np.abs((T-T_old)/T_old)), np.nanmax(np.abs((S-S_old)/S_old))

    #F_old = np.copy(F)
    #T_old = np.copy(T)
    #S_old = np.copy(S)

    
    ##plt.plot(F, '.')
    ##plt.show()
    
    ##plt.plot(T, '.')
    ##plt.show()
    
    ##plt.plot(S, '.')
    ##plt.show()
    
T = np.ones(Tlen)
S = np.zeros(Slen)
    
for i in range(10):
    
    F = np.bincount(Fidx, weights*flux*(T[Tidx]-S[Sidx]), minlength=Flen)/np.bincount(Fidx, weights*(T[Tidx]-S[Sidx])**2, minlength=Flen) # F    
    T = np.bincount(Tidx, weights*(flux+F[Fidx]*S[Sidx])*F[Fidx], minlength=Tlen)/np.bincount(Tidx, weights*(F[Fidx])**2, minlength=Tlen) # F  
    S = -np.bincount(Sidx, weights*(flux-F[Fidx]*T[Tidx])*F[Fidx], minlength=Slen)/np.bincount(Sidx, weights*(F[Fidx])**2, minlength=Slen) # F  
    
    F[np.isnan(F)] = 0.
    T[np.isnan(T)] = 0.
    S[np.isnan(S)] = 0.
    
    if i > 0:
        print np.nanmax(np.abs((F-F_old)/F_old)), np.nanmax(np.abs((T-T_old)/T_old)), np.nanmax(np.abs((S-S_old)/S_old))

    F_old = np.copy(F)
    T_old = np.copy(T)
    S_old = np.copy(S)

    
    plt.plot(F, '.')
    plt.show()
    
    plt.plot(T, '.')
    plt.show()
    
    plt.plot(S, '.')
    plt.show()
    
plt.plot(F, '.')
plt.show()

plt.plot(T, '.')
plt.show()

plt.plot(S, '.')
plt.show()
    
fit = F[Fidx]*(T[Tidx]-S[Sidx])

here = staridx == 0
plt.plot(jdmid[here], flux[here], '.')
plt.plot(jdmid[here], fit[here], '.')
plt.show()

flux = flux/fit
eflux = eflux/fit

phase = np.mod(jdmid-t0, period)/period

here = here&np.isfinite(flux)


freq, chisq = BLS(jdmid[here], flux[here], eflux[here])

plt.subplot(211)
plt.plot(phase[here], flux[here], '.')
plt.subplot(212)
plt.plot(freq, chisq)
plt.axvline(1/period)
plt.show()

binidx = np.ravel_multi_index((dayidx, lstidx//50), (np.amax(dayidx)+1, 270))

count = index_statistics(binidx[here], jdmid[here], statistic='count')
bin_jd = index_statistics(binidx[here], jdmid[here], statistic='mean')
bin_flux = index_statistics(binidx[here], flux[here], statistic='mean')
bin_eflux = index_statistics(binidx[here], flux[here], statistic='std')/np.sqrt(count)

here = count ==50
freq, chisq = BLS(bin_jd[here], bin_flux[here], bin_eflux[here])
phase = np.mod(bin_jd-t0, period)/period
plt.subplot(211)
plt.plot(bin_jd[here], flux[here], '.')
plt.subplot(212)
plt.plot(freq, chisq)
plt.axvline(1/period)
plt.show()





