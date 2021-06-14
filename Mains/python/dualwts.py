#! /usr/bin/env python

import sys
from pylab import *
import pyctf
from pyctf.samiir import *
from pyctf.st import st
from pyctf.sensortopo import sensortopo
from plotst import plotst
from dual import dual
from math import acos

ds = pyctf.dsopen(sys.argv[1])
srate = ds.getSampleRate()
M = ds.getNumberOfPrimaries()
T = ds.getNumberOfSamples()

lo = 15.
hi = 30.
f = mkiir(lo, hi, srate)

print 'reading and filtering'
d = zeros((M, T), 'd')
for m in range(M):
	ch = ds.getPrimaryIndex(m)
	x = ds.getDsData(0, ch)
	d[m, :] = bdiir(x, f)
print 'done'

w, c = ds.readwts("click0,targ4")
#w, c = ds.readwts("targ4")

x = concatenate((w[0, :], w[1, :]))
mx = max(x)
mn = min(x)
zrange = [mn, mx]
print zrange

topo = sensortopo(ds)

subplot(121)
topo.plot(w[0, :], zrange = zrange)
subplot(122)
topo.plot(w[1, :], zrange = zrange)
figure()

def vangle(v1, v2):
	"""Calculate the angle between v1 and v2 using the relationship:
	    cos(theta) = <v1,v2> / (||v1|| * ||v2||)
	"""
	return acos(dot(v1, v2) / sqrt(dot(v1, v1) * dot(v2, v2)))

w0 = w[0, :]
w1 = w[1, :]
print vangle(w0, w1) * 180. / pi

w = dual(w)
w0 = w[0, :]
w1 = w[1, :]

# Convert frequencies in Hz into rows of the ST, given sampling rate and length.

def freq(f, n, srate = srate):
	return int(f * n / srate + .5)

m = ds.marks['click0']

print 'processing trials'
S = 0.
pre = .5 * srate
post = 1.5 * srate
n = pre + post
flo = freq(lo, n)
fhi = freq(hi, n)
a = b = 0.
for tr, t in m:
	samp = ds.getSampleNo(t)
	s = d[:, samp - pre : samp + post]
	s0 = dot(w0, s)
	s1 = dot(w1, s)
	a += s0
	b += s1

	z0 = st(s0, flo, fhi)
	z1 = st(s1, flo, fhi)

	z0 /= abs(z0)
	z1 /= abs(z1)

	S += conjugate(z1) * z0 # relative phase

a /= len(m)
b /= len(m)
plv = abs(S) / len(m)
relph = angle(S)

S = S.imag
#S = S.real

print plv.min(), plv.max()
print relph.min(), relph.max()

x = linspace(-pre / srate, post / srate, len(a))
plot(x, a, x, b)
plotst(relph, srate, start = -pre, end = post, lo = lo, hi = hi, logscale = False, cmap = cm.hsv, range = [-pi, pi])
plotst(plv, srate, start = -pre, end = post, lo = lo, hi = hi, logscale = False)
plotst(S, srate, start = -pre, end = post, lo = lo, hi = hi, logscale = False, cmap = cm.hsv)
show()
