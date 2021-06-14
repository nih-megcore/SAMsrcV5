#! /usr/bin/env python

import sys, os
import pyctf
import nolteFwd
import gtk
from pylab import array, linspace, sin, cos, pi, figure, show, hold, clf, ion, savefig, colorbar
from time import sleep, time

def gtk_idle(s):
    t = time() + s
    while time() < t:
        gtk.main_iteration(0)
        sleep(.01)

ds = pyctf.dsopen(sys.argv[1])
nolteFwd.dsopen(sys.argv[1])

pos = (0., -5., 7.)

f = open("/tmp/moo", "w")

s = pyctf.sensortopo(ds)
ion()
ctr = 0
for theta in linspace(0., 2. * pi, 50)[:-1]:
    print theta

    ori = (sin(theta), 0., cos(theta))

    b = nolteFwd.doFwd(pos, ori)
    f.write("%g %g\n" % (b.min(), b.max()))

    clf()
    s.plot(b, zrange = [-1.5e-5, 1.5e-5])
    savefig("r%03d.png" % ctr)
    ctr += 1
    gtk_idle(.1)

f.close()
