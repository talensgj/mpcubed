#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

from fourierfuncs import lst_trend, fourier_fit
from BLS_ms import BLS
from package.statistics import statistics
from scipy.signal import lombscargle

P = 2.21857312
Tp = 2454037.612

with h5py.File('/data2/talens/2015Q2/LPE/red0_2015Q2LPE.hdf5') as f:
    grp = f['data/807144']
    jdmid = grp['jdmid'].value
    mag0 = grp['mag0'].value
    emag0 = grp['emag0'].value
    lst = grp['lst'].value
    lstseq = grp['lstseq'].value
    nobs = grp['nobs'].value
    
#with h5py.File('/data2/talens/2015Q2/LPE/new_solver/red0_2015Q2LPE.hdf5') as f:
    #grp = f['data/807144']
    #jdmid = grp['jdmid'].value
    #mag0 = grp['mag0'].value
    #emag0 = grp['emag0'].value
    #lst = grp['lst'].value
    #lstseq = grp['lstseq'].value
    #nobs = grp['nobs'].value

# Select full bins.
select = (nobs == 50)
jdmid = jdmid[select]
mag0 = mag0[select]
emag0 = emag0[select]
lst = lst[select]
lstseq = lstseq[select]

emag0 = emag0/np.sqrt(50)

jdmid = jdmid.astype('float64')
mag0 = mag0.astype('float64')
freqs = np.linspace(1/np.ptp(jdmid), 12., 1000)

mag0 = mag0 - np.nanmean(mag0)

for i in range(10):

    pgram = lombscargle(jdmid, mag0, 2*np.pi*freqs)

    plt.plot(freqs, pgram)
    plt.show()

    arg = np.argmax(pgram)
    pars, fit = fourier_fit(jdmid, mag0, freqs[arg], 1, 1/emag0**2)

    plt.plot(mag0, '.')
    plt.plot(fit, '.')
    plt.show()

    mag0 = mag0 - fit

exit()
# Reduction method onvolving 1 fourier fit.
dayidx = lstseq//270
pars1, pars2, fit = lst_trend(dayidx, lst, mag0, 1/24., 5, 1/emag0**2)

chisq = (mag0 - fit)**2/emag0**2
print np.sum(chisq)

mag_1 = mag0 - fit

# Reduction method involving 2 fourier fits.
fit2 = np.zeros(len(mag0))
for i in range(5):
    pars1, fit1 = fourier_fit(jdmid, mag0 - fit2, 1/np.ptp(jdmid), 35, 1/emag0**2)
    pars2, fit2 = fourier_fit(lst, mag0 - fit1, 1/24., 5, 1/emag0**2)
    print pars2

chisq = (mag0 - fit1 - fit2)**2/emag0**2
print np.sum(chisq)

plt.plot(jdmid, mag0, '.')
plt.plot(jdmid, fit, '.')
plt.plot(jdmid, fit1, '.')
plt.show()

mag_2 = mag0 - fit1 - fit2

# Plot the resulting lightcurves.
phase = np.mod((jdmid - Tp)/P, 1.)
phase = np.mod(phase + .5, 1) - .5

plt.figure(figsize=(18, 5))
plt.plot(phase, mag_1, '.')
plt.plot(phase, mag_2, '.')
plt.xlim(-.5, .5)
plt.ylim(.1, -.1)
plt.xlabel('Phase')
plt.ylabel(r'$\Delta m$')
plt.tight_layout()
plt.show()

freq_1, dchisq_1, depth, hchisq = BLS(jdmid, mag_1, emag0)
freq_2, dchisq_2, depth, hchisq = BLS(jdmid, mag_2, emag0)

plt.figure(figsize=(18, 10))
plt.plot(freq_1, dchisq_1)
plt.plot(freq_2, dchisq_2)
plt.axvline(1/P, c='k')
plt.ylim(0, 1500)
plt.xlabel('Frequency [1/day]')
plt.ylabel(r'$\Delta\chi^2$')
plt.tight_layout()
plt.show()
