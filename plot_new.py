#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

from index_functions import index_statistics

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

plt.figure(figsize=(16,8))

for camera in ['LPN', 'LPE', 'LPS', 'LPW', 'LPC']:

    with h5py.File('/data2/talens/Jul2015/fLC_20150714%s.hdf5'%camera) as f:
        
        try:
            lc = f['data/807144'].value
        except:
            continue

    with h5py.File('/data2/talens/Jul2015/red_20150714%s.hdf5'%camera) as f:

        rc = f['data/807144'].value
        rc2 = f['data2/807144'].value

    norm = np.mean(lc['flux0']/(rc['camtrans0']))
    jdmid = lc['jdmid'] - np.floor(lc['jdmid'])

    ax = plt.subplot(311)
    plt.plot(jdmid, lc['flux0']/norm, '.')
    plt.plot(jdmid, rc['camtrans0'], '.')

    plt.ylabel('Flux [ADU]')

    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jdmid, rc['cflux0']/norm, '.')
    plt.plot(jdmid, rc2['skytrans0'], '.')

    plt.ylabel('cFlux [ADU]')

    plt.subplot(313, sharex=ax, sharey=ax)
    plt.plot(jdmid, rc2['scflux0']/norm, '.')
        
    #x = index_statistics(lc['lstidx']//50, jdmid, statistic='mean')
    #y = index_statistics(lc['lstidx']//50, rc2['scflux0']/norm, statistic='mean')
    #yerr = index_statistics(lc['lstidx']//50, rc2['scflux0']/norm, statistic='std')/np.sqrt(index_statistics(lc['lstidx']//50, rc2['scflux0'], statistic='count'))

    #plt.errorbar(x, y, yerr=yerr, fmt='.', ecolor='k')

    plt.xlim(0.25, 0.75)
    plt.xlabel('Time [JD]')
    plt.ylabel('scFlux [ADU]')

plt.tight_layout()
plt.show()


