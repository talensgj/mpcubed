#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import matplotlib.pyplot as plt

from core import intrapix

#ip = intrapix.IntraPixel()
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5')

#cf = intrapix.CameraFile('/data2/talens/Jul2015/camip_20150716LPC.hdf5')
#cf.visualize(wrap=True)
#cf.correct('/data2/talens/Jul2015/fLC_20150714LPC.hdf5')

#with h5py.File('/data2/talens/Jul2015/fLC_20150714LPC.hdf5', 'r') as f, h5py.File('/data2/talens/Jul2015/red_20150714LPC.hdf5', 'r') as g:
    #ascc = f['data'].keys()
    
    #for i in range(0, len(ascc), 50):
        #lc = f['data/'+ascc[i]]
        #rc = g['data/'+ascc[i]]
        
        #x = range(len(lc['flux0']))
        
        #ax = plt.subplot(311)
        #plt.title(ascc[i])
        #plt.plot(x, lc['flux0'], '.')
        
        #plt.subplot(312, sharex=ax)
        #plt.plot(x, rc['camtrans0']*rc['intrapix0'], '.')

        #plt.subplot(313, sharex=ax)
        #plt.plot(x, rc['ipc_flux0'], '.')
        
        #plt.show()
        #plt.close()
