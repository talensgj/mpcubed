#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

with h5py.File('HF_params.hdf5', 'r') as f:
    
    grp = f['nlst5']
    chisq0 = grp['chisq'].value
    
    grp = f['nlst5min2']
    chisq = grp['chisq'].value
    pars = grp['pars'].value
    
print pars[:,0]
    
plt.hist(chisq0/chisq, bins=np.linspace(0,1,21))
plt.show()
