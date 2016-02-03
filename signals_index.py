#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package import IO
from package.statistics import statistics

def create_sample(ascc, ra, dec, vmag, nobs, jdstart):
    
    np.random.seed(19910909)
    
    # Select stars with some minimum number of datapoints.
    select = (nobs > 6*13500)
    ascc = ascc[select]
    ra = ra[select]
    dec = dec[select]
    vmag = vmag[select]
    nobs = nobs[select]
    jdstart = jdstart[select]

    # Select a sample of N stars.
    p = 1/vmag**2
    p = p/np.sum(p)
    sample = np.random.choice(np.arange(len(ascc)), 500, replace=False, p = p)

    ascc = ascc[sample]
    ra = ra[sample]
    dec = dec[sample]
    vmag = vmag[sample]
    nobs = nobs[sample]
    jdstart = jdstart[sample]

    # Generate a set of plausible transit parameters.
    delta = .005 + .025*np.random.rand(500)
    P = 1 + 14*np.random.rand(500)
    Tp = jdstart + P*np.random.rand(500)
    mu = (1.8/24)*P**(1./3) # Assuming solar mass and radius.

    return ascc, ra ,dec, vmag, nobs, jdstart, delta, P, Tp, mu

filelist = glob.glob('/data2/talens/inj_signals/reference/fLC_*.hdf5')

ascc = np.array([])
ra = np.array([])
dec = np.array([])
vmag = np.array([])
nobs = np.array([])
jdstart = np.array([])

for filename in filelist:
    
    f = IO.fLCfile(filename)
    fields = ['ascc', 'ra', 'dec', 'vmag', 'nobs', 'jdstart']
    data = f.read_header(fields)
    
    ascc = np.append(ascc, data[0])
    ra = np.append(ra, data[1])
    dec = np.append(dec, data[2])
    vmag = np.append(vmag, data[3])
    nobs = np.append(nobs, data[4])
    jdstart = np.append(jdstart, data[5])

ascc, args, idx = np.unique(ascc, return_index=True, return_inverse=True)
ra = ra[args]
dec = dec[args]
vmag = vmag[args]
nobs = np.bincount(idx, nobs)
jdstart = statistics.idxstats(idx, jdstart, statistic=np.amin)

plt.scatter(ra, dec, c=jdstart)
plt.show()

sample = create_sample(ascc, ra, dec, vmag, nobs, jdstart)

plt.subplot(221)
plt.hist(sample[6])
plt.subplot(222)
plt.hist(sample[7])
plt.subplot(223)
plt.hist(sample[8])
plt.subplot(224)
plt.hist(sample[9])
plt.show()


with h5py.File('/data2/talens/inj_signals/signals/signals_index.hdf5') as f:
    f.create_dataset('ascc', data = sample[0])
    f.create_dataset('ra', data = sample[1])
    f.create_dataset('dec', data = sample[2])
    f.create_dataset('vmag', data = sample[3])
    f.create_dataset('nobs', data = sample[4])
    f.create_dataset('jdstart', data = sample[5])
    f.create_dataset('delta', data = sample[6])
    f.create_dataset('P', data = sample[7])
    f.create_dataset('Tp', data = sample[8])
    f.create_dataset('mu', data = sample[9])
