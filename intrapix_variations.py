#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from core import intrapix_weights
from core import skytrans_weights

#ip = intrapix_weights.IntraPixel()
#ip.calculate('/data2/talens/3mEast/fLC_20150611LPE.hdf5', camfile='/data2/talens/3mEast/camip_aper_w_xy_20150611LPE.hdf5')

cf = intrapix_weights.CameraFile('/data2/talens/3mEast/camip_aper_w_xy_20150611LPE.hdf5')
#cf.visualize()
cf.correct('/data2/talens/3mEast/fLC_20150611LPE.hdf5', redfile='/data2/talens/3mEast/red_aper_w_xy_20150611LPE.hdf5')

