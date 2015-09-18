#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import numpy as np

from get_intrapix_1cam_classdev import IntraPixel, CameraFile
from get_skytrans_1cam_classdev import SkyTransmission, SkyFile

filelist = glob.glob('/data2/talens/Jul2015/fLC_201507??LPC.hdf5')
filelist = np.sort(filelist)
print filelist
print len(filelist)

#ip = IntraPixel()
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5')

#cf = CameraFile('/data2/talens/Jul2015/camip_20150716LPC.hdf5')
#for filename in filelist:
    #cf.correct(filename)
    
#st = SkyTransmission()
for filename in filelist:
    #st.calculate(filename)
    head, tail = os.path.split(filename)
    tail = 'sky_'+tail.rsplit('_')[-1]
    skyfile = os.path.join(head, tail)
    sf = SkyFile(skyfile)
    sf.correct(filename)
