#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package.red_decor_vmag import CoarseDecorrelation


filename = '/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5'
sysfile = '/data2/talens/2015Q2/LPE/sys0_201506ALPE_vmag.hdf5'

CoarseDecorrelation(filename, 0, sysfile=sysfile, outer_maxiter=5)

#data = ['/data2/talens/2015Q2/LPN/fLC_201504ALPN.hdf5',
        #'/data2/talens/2015Q2/LPN/fLC_201504BLPN.hdf5',
        #'/data2/talens/2015Q2/LPN/fLC_201505ALPN.hdf5',
        #'/data2/talens/2015Q2/LPN/fLC_201505BLPN.hdf5',
        #'/data2/talens/2015Q2/LPN/fLC_201506ALPN.hdf5',
        #'/data2/talens/2015Q2/LPN/fLC_201506BLPN.hdf5']

#pool = mp.Pool(processes = 6)
#for i in range(6):
    #pool.apply_async(correct, args = (data[i],))
#pool.close()
#pool.join()
