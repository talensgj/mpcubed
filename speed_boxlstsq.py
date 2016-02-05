#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from package.models import transit
from boxlstsq_ms import boxlstsq_ms

nstars = 100
ndays = 90

jdmid = np.linspace(0, ndays, ndays*270)

mag0 = np.zeros((nstars, ndays*270))
weights = np.ones((nstars, ndays*270))

P = .6 + (ndays/9. - .6)*np.random.rand(nstars)
Tp = P*np.random.rand(nstars)
delta = .005 + .025*np.random.rand(nstars)
eta = (1.8/24)*P**(1./3)

for i in range(nstars):
    print i
    mag0[i] = transit.softened_box_model(jdmid, P[i], Tp[i], delta[i], eta[i])
    
plt.imshow(mag0, aspect='auto', interpolation='None')
plt.show()

freq, dchisq, depth, hchisq = boxlstsq_ms(jdmid, mag0.T, weights.T)

args = np.argmax(dchisq, axis=0)

plt.plot(P, 1/freq[args] - P, '.')
plt.show()
