#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import matplotlib.pyplot as plt

from core import skytrans

#st = skytrans.SkyTransmission()
#st.calculate('/data2/talens/Jul2015/fLC_20150714LPC.hdf5', '/data2/talens/Jul2015/red_20150714LPC.hdf5')

#sf = skytrans.SkyFile('/data2/talens/Jul2015/sky_20150714LPC.hdf5')
#sf.visualize(wrap=True)
#sf.correct()

with h5py.File('/data2/talens/Jul2015/red_20150714LPC.hdf5', 'r') as f:
    ascc = f['data'].keys()
    
    for i in range(0, len(ascc), 50):
        lc = f['data/'+ascc[i]]
        rc = f['data2/'+ascc[i]]
        
        x = range(len(lc['ipc_flux0']))
        
        ax = plt.subplot(311)
        plt.title(ascc[i])
        plt.plot(x, lc['ipc_flux0'], '.')
        
        plt.subplot(312, sharex=ax)
        plt.plot(x, rc['skytrans0'], '.')

        plt.subplot(313, sharex=ax)
        plt.plot(x, rc['sipc_flux0'], '.')
        
        plt.show()
        plt.close()
