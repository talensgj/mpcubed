#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package import detrend
from package import new_detrend

#detrend.create_baseline('201506', 'LPE', 0, '/data2/mascara/LaPalma', '/data2/talens/LongBaselines')

#obj = detrend.CoarseDecorrelation('/data2/talens/2015Q2/LPW/fLC_201506ALPW.hdf5', 0, '/data2/talens/2015Q2/LPW/test_master.hdf5', outer_maxiter=1)
#obj.run()

new_detrend.CoarseDecorrelation('/data2/talens/2015Q2/LPW/fLC_201506ALPW.hdf5', 0, '/data2/talens/2015Q2/LPW/test.hdf5', outer_maxiter=1)

#obj = detrend.SysCorr('/data2/talens/2015Q2/LPW/fLC_201506BLPW.hdf5', 0, outfile = '/data2/talens/2015Q2/LPW/tmp_test.hdf5')
#obj.run()

#import glob
#import numpy as np

#filelist = glob.glob('/data2/talens/2015Q2/LPW/tmp0_*')
#filelist = np.sort(filelist)

#detrend.make_quarterfile(filelist, redfile = '/data2/talens/2015Q2/LPW/red_test.hdf5')
