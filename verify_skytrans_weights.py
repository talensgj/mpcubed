#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import matplotlib.pyplot as plt

from core import skytrans_weights

#st = skytrans_weights.SkyTransmission()
#st.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix2.hdf5', '/data2/talens/Jul2015/verify_skytrans1.hdf5')

#sf = skytrans_weights.SkyFile('/data2/talens/Jul2015/verify_skytrans1.hdf5')
#sf.visualize(wrap=True)
#sf.correct('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix2.hdf5')

with h5py.File('/data2/talens/Jul2015/verify_intrapix2.hdf5', 'r') as f:
    ascc = f['data'].keys()
    
    for i in range(0, len(ascc), 50):
        lc = f['data/'+ascc[i]]
        rc = f['data2/'+ascc[i]]
        
        plt.subplot(311)
        plt.title(ascc[i])
        plt.plot(lc['ipcmag0'], '.')
        
        plt.subplot(312)
        plt.plot(rc['skytrans0'], '.')

        plt.subplot(313)
        plt.plot(rc['scmag0'], '.')
        
        plt.show()
        plt.close()
