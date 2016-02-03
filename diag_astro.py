#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from mascara.reduction import dbquery

import matplotlib.pyplot as plt

tablename = 'LPEastrometry'

myq = dbquery.MascDB()
myq.connect()

cur = myq.dbconnect.cursor()
cur.execute("SELECT lstseq, alt, az, tho, xo, yo FROM {}".format(tablename))
data = cur.fetchall()

myq.disconnect()

data = np.array(data)
data = data.T

plt.subplot(511)
plt.plot(data[1], '.')

plt.subplot(512)
plt.plot(data[2], '.')

plt.subplot(513)
plt.plot(data[3], '.')

plt.subplot(514)
plt.plot(data[4], '.')

plt.subplot(515)
plt.plot(data[5], '.')

plt.show()

# Further coefficients.
myq = dbquery.MascDB()
myq.connect()

cur = myq.dbconnect.cursor()
cur.execute("SELECT lstseq, crx, cry, cd11, cd12, cd21, cd22 FROM {}".format(tablename))
data = cur.fetchall()

myq.disconnect()

data = np.array(data)
data = data.T

plt.subplot(211)
plt.plot(data[1], '.')

plt.subplot(212)
plt.plot(data[2], '.')

plt.show()

plt.subplot(221)
plt.plot(data[3], '.')

plt.subplot(222)
plt.plot(data[4], '.')

plt.subplot(223)
plt.plot(data[5], '.')

plt.subplot(224)
plt.plot(data[6], '.')

plt.show()

# Polynomial coeffcients.
myq = dbquery.MascDB()
myq.connect()

cur = myq.dbconnect.cursor()
cur.execute("SELECT lstseq, xcoef0, xcoef1, xcoef2, xcoef3 FROM {}".format(tablename))
data = cur.fetchall()

myq.disconnect()

data = np.array(data)
data = data.T

plt.subplot(221)
plt.plot(data[1], '.')

plt.subplot(222)
plt.plot(data[2], '.')

plt.subplot(223)
plt.plot(data[3], '.')

plt.subplot(224)
plt.plot(data[4], '.')

plt.show()
