#!/usr/bin/env python
# -*- coding: utf-8 -*-

from core import intrapix

#ip = intrapix.IntraPixel()
#ip.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix1.hdf5')

cf = intrapix.CameraFile('/data2/talens/Jul2015/verify_intrapix1.hdf5')
cf.visualize(wrap=True)
cf.correct('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix2.hdf5')
