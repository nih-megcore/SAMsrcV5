#! /usr/bin/env python

import sys
from time import sleep
from pylab import *
import pyctf
from pyctf.sensortopo import sensortopo
from numpy.fft import fft
from gtk import events_pending, main_iteration

def update():
	while events_pending():
		main_iteration()

sum = 0.
N = 0

ds = pyctf.dsopen(sys.argv[1])
srate = ds.getSampleRate()
M = ds.getNumberOfPrimaries()
T = ds.getNumberOfSamples()

#for dsname in sys.argv[1:]:
#        ds = pyctf.dsopen(dsname)
#        srate = ds.getSampleRate()
#        M = ds.getNumberOfPrimaries()
#        T = ds.getNumberOfSamples()

marker = 'click0'
pre = 1.1 * srate
post = 2.1 *srate

print 'reading data'
d = zeros((M, T), 'd')
for i in range(M):
	ch = ds.getPrimaryIndex(i)
	x = ds.getDsData(0, ch)
	d[i, :] = x
print 'done'

print 'processing trials'
m = ds.marks[marker]
for tr, t in m:
	samp = ds.getSampleNo(t)
	s = d[:, samp - pre : samp + post].copy()
	s -= s.mean()
	s /= s.std()
	S = fft(s)
	sum += S.imag * S.imag + S.real * S.real
	N += 1

sum /= N

n = pre + post
fr = srate / n

ion()
topo = sensortopo(ds)

while 1:
	for f in linspace(0, 120. / fr, 25):
		if f == 0: f += 1
		clf()
		a = add.reduce(sum[:,f:f+5], axis=1)
		topo.plot(a)
		title("%g Hz" % (f * fr))
		update()
		sleep(.5)
