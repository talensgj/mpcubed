#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import numpy as np

from get_intrapix_1cam_classdev import IntraPixel, CameraFile
from get_skytrans_1cam_classdev import SkyTransmission, SkyFile

filelist = glob.glob('/data2/talens/Jul2015/fLC_201507??LPW.hdf5')
filelist = np.sort(filelist)

ip = IntraPixel()
ip.calculate('/data2/talens/Jul2015/fLC_20150716LPW.hdf5')

cf = CameraFile('/data2/talens/Jul2015/camip_20150716LPW.hdf5')
for filename in filelist:
    cf.correct(filename)
    
st = SkyTransmission()
for filename in filelist:
    st.calculate(filename)
    sf = SkyFile(st.skyfile)
    sf.correct(filename)
