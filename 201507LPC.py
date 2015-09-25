#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import numpy as np

from get_intrapix_weights import IntraPixel, CameraFile
from get_skytrans_weights import SkyTransmission, SkyFile

filelist = glob.glob('/data2/talens/Jul2015/fLC_201507??LPC.hdf5')
filelist = np.sort(filelist)

ip = IntraPixel()
ip.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5')

cf = CameraFile('/data2/talens/Jul2015/coarsecam_20150716LPC.hdf5')
for filename in filelist:
    cf.correct(filename)
    
st = SkyTransmission()
for filename in filelist:
    print filename
    st.calculate(filename)
    sf = SkyFile(st.skyfile)
    sf.correct(filename)
