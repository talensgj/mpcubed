#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

path = '/data2/talens/inj_signals/signals/'

with h5py.File(path+'HF_njd5_2BL_nlst5_24h.hdf5', 'r') as f:
    
    pars1 = f['pars'].value
    chisq1 = f['chisq'].value
 
with h5py.File(path+'AF_njd5_2BL.hdf5', 'r') as f:
    
    #pars2 = f['pars'].value
    chisq2 = f['chisq'].value

plt.plot(chisq1, chisq1 - chisq2, '.')
plt.show()
