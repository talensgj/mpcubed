#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from core import index_functions

def phase(time, period, time_ref=0., fold=True):
    
    phase = (time - time_ref)/period
    
    if fold:
        phase = np.mod(phase, 1.)
        
    return phase
    
def bindata(binidx, time, data, minpoints=2):
    
    npoints = index_functions.index_statistics(binidx, None, statistic='count')
    bin_time = index_functions.index_statistics(binidx, time, statistic='mean')
    bin_data = index_functions.index_statistics(binidx, data, statistic='mean')
    bin_edata = index_functions.index_statistics(binidx, data, statistic='std')/np.sqrt(npoints)
    
    here = (npoints > minpoints)
    
    bin_time = bin_time[here]
    bin_data = bin_data[here]
    bin_edata = bin_edata[here]
    
    return bin_time, bin_data, bin_edata
    
def mad(data, K=1.4826, axis=None):
    
    med = np.nanmedian(data, axis=axis, keepdims=True)
    mad = np.nanmedian(np.abs(data-med), axis=axis)
    
    return K*mad
