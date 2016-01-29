#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('/data2/talens/inj_signals/signals/red0_2015Q2LPE.hdf5', 'r') as f:
    jdmid = f['data/807144/jdmid'].value
    mag0 = f['data/807144/mag0'].value
    emag0 = f['data/807144/emag0'].value
    nobs = f['data/807144/nobs'].value
    
select = (nobs == 50)
jdmid = jdmid[select]
mag0 = mag0[select]
emag0 = emag0[select]

from BLS_ms import BLS
from boxlstsq import boxlstsq
from boxlstsq_dev import boxlstsq_dev

freq1, dchisq1, depth1, hchisq1 = BLS(jdmid, mag0, emag0)
freq2, dchisq2, depth2, hchisq2 = boxlstsq(jdmid, mag0, emag0)
freq3, dchisq3, depth3, hchisq3 = boxlstsq_dev(jdmid, mag0, emag0)

plt.subplot(311)
plt.plot(freq1, dchisq1)
plt.plot(freq2, dchisq3)
plt.plot(freq3, dchisq3)

plt.subplot(312)
plt.plot(freq1, depth1)
plt.plot(freq2, depth3)
plt.plot(freq3, depth3)

plt.subplot(313)
plt.plot(freq1, hchisq1)
plt.plot(freq2, hchisq3)
plt.plot(freq3, hchisq3)

plt.show()
