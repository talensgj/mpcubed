#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package import detrend
from package import red_decor_dev

obj = detrend.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPW/ref_06A_iter1.hdf5', outer_maxiter=1)
obj.run()

red_decor_dev.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPW/test_06A_iter1.hdf5', outer_maxiter=1)

obj = detrend.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPW/ref_06A.hdf5')
obj.run()

red_decor_dev.CoarseDecorrelation('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPW/test_06A.hdf5')
