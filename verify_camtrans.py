#!/usr/bin/env python
# -*- coding: utf-8 -*-

from core import camtrans

#ct = camtrans.CameraTransmission()
#ct.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_camtrans1.hdf5')

cf = camtrans.CameraFile('/data2/talens/Jul2015/verify_camtrans1.hdf5')
#cf.visualize(wrap=True)
cf.correct('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_camtrans2.hdf5')
