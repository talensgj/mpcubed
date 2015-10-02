#!/usr/bin/env python
# -*- coding: utf-8 -*-

from core import skytrans

#st = skytrans.SkyTransmission()
#st.calculate('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix2.hdf5', '/data2/talens/Jul2015/verify_skytrans1b.hdf5')

sf = skytrans.SkyFile('/data2/talens/Jul2015/verify_skytrans1b.hdf5')
#sf.visualize(wrap=True)
sf.correct('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', '/data2/talens/Jul2015/verify_intrapix2.hdf5')
