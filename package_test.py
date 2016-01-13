#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package import detrend
from package import red_decor_dev

filepath = '/data2/mascara/LaPalma'
outpath = '/data2/talens/LongBaselines'

detrend.create_baseline('201506', 'LPE', 0, filepath, outpath)

#obj = detrend.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPE/ref_06A_iter1.hdf5', outer_maxiter=1)
#obj.run()

#red_decor_dev.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPE/test_06A_iter1.hdf5', outer_maxiter=1)

#obj = detrend.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPE/ref_06A.hdf5')
#obj.run()

#red_decor_dev.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPE/test_06A.hdf5')
