#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import matplotlib.pyplot as plt

from core import intrapix_weights

#ip = intrapix_weights.IntraPixel()
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix1.hdf5')

cf = intrapix_weights.CameraFile('/data2/talens/Jul2015/verify_intrapix1.hdf5')
cf.visualize(wrap=True)
cf.correct('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix2.hdf5')

#with h5py.File('/data2/talens/Jul2015/verify_intrapix2.hdf5', 'r') as f:
    #ascc = f['data'].keys()
    
    #for i in range(0, len(ascc), 50):
        #lc = f['data/'+ascc[i]]
        
        #ax = plt.subplot(311)
        #plt.title(ascc[i])
        #plt.plot(lc['mag0'], '.')
        
        #plt.subplot(312, sharex=ax)
        #plt.plot(lc['camtrans0']+lc['intrapix0'], '.')

        #plt.subplot(313, sharex=ax)
        #plt.plot(lc['ipc_mag0'], '.')
        
        #plt.show()
        #plt.close()
