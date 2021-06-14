#! /usr/bin/env python

import sys
import pyctf
import nolteFwd
import numpy as np
import matplotlib.pyplot as plt

ds = pyctf.dsopen(sys.argv[1])
nolteFwd.dsopen(sys.argv[1])
print("done")

pos = (0., -5., 7.)
##theta = 45 * 3.141 / 180
##ori = (np.sin(theta), 0., np.cos(theta))

s = pyctf.sensortopo(ds)

C = ds.readcov("beta,2,20-30Hz/Orient.cov")
Cinv = np.linalg.inv(C)

plt.subplot(1, 2, 1)
b = nolteFwd.doFwdMoment(pos, Cinv)
s.plot(b, zrange = [b.min(), b.max()])
#plt.show()

plt.subplot(1, 2, 2)
x = Cinv.dot(b)
h = x / x.dot(b)
s.plot(h, zrange = [h.min(), h.max()])

plt.show()
