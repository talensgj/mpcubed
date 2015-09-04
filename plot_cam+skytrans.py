#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/Jul2015/red_20150716LPC.hdf5') as f:
    
    ascc = f['data'].keys()
    Ntot = 0

    Nbad1 = 0
    Nbad2 = 0
    
    for sid in ascc:
        #Ntot += len(f['data/%s'%sid]['flags'])
        #Nbad1 += np.sum(f['data/%s'%sid]['flags']>0)
        #Nbad2 += np.sum(f['data2/%s'%sid]['flags']>0)

        print f['data/%s'%sid]['flags'][f['data/%s'%sid]['flags']>0]

    print Nbad1, Ntot, float(Nbad1)/Ntot
    print Nbad2, Ntot, float(Nbad2)/Ntot
exit()
plt.figure(figsize=(16,10))

for camera, color in zip(['LPC', 'LPE', 'LPW'], ['c', 'r', 'y']):
    with h5py.File('/data2/talens/Jul2015/red_20150716%s.hdf5'%camera, 'r') as f, h5py.File('/data2/talens/Jul2015/fLC_20150716%s.hdf5'%camera, 'r') as g:
        
        raw = g['data/807144']
        cam = f['data/807144']
        sky = f['data2/807144']

        jd = raw['jdmid']-np.floor(raw['jdmid'])
        corr = sky['skytrans0']
        corr /= np.median(corr)

        print np.sum(cam['flags'])
        print np.sum(sky['flags'])

        ax = plt.subplot(511)
        plt.plot(jd, raw['flux0']/np.median(raw['flux0']), '.', c=color, label=camera)
        plt.ylim(0,1.5)
        plt.ylabel('Flux')
        plt.subplot(512, sharex=ax)
        plt.plot(jd, cam['camtrans0']/np.median(cam['camtrans0']), '.', c=color)
        plt.ylim(0,1.5)
        plt.ylabel('Trans')
        plt.subplot(513, sharex=ax)
        plt.plot(jd, cam['cflux0']/np.median(cam['cflux0']), '.', c=color)
        plt.ylim(0,1.5)
        plt.ylabel('cFlux')
        plt.subplot(514, sharex=ax)
        plt.plot(jd, corr, '.', c=color)
        plt.ylim(0,1.5)
        plt.ylabel('Sky')
        plt.subplot(515, sharex=ax)
        plt.plot(jd, sky['scflux0']/np.median(sky['scflux0']), '.', c=color)
        plt.ylim(0,1.5)
        plt.xlim(0.25, .75)
        plt.ylabel('scFlux')
        plt.xlabel('Time [JD-2457218]')
    
plt.sca(ax)
plt.legend(loc=2)
plt.tight_layout()
plt.show()
