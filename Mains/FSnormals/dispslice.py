#! /usr/bin/env python

import sys
import nibabel
from pylab import *
from matplotlib import lines, cm
from matplotlib.collections import LineCollection
from bisect import bisect

mri = sys.argv[1]
surf = sys.argv[2]

im = nibabel.load(mri)
data = im.get_data()
toLPI = im.affine
fromLPI = inv(toLPI)

# Load the surface verts and normals, convert from PRI m to LPI mm

surf = loadtxt(surf)
surf = surf[:,:6]
for i in range(len(surf)):
    v = surf[i][:3] * 1000.
    n = surf[i][3:]
    v = -v[1], v[0], v[2]
    n = -n[1], n[0], n[2]
    surf[i][:3] = v
    surf[i][3:] = n

def find_verts(axname, a, b):
    """Return all the vertices that are between a and b along ax."""

    ax = 'LPI'.index(axname)
    idx = surf[:,ax].argsort()  # create a sort based on the given column
    v = surf[idx]               # apply to entire array
    z = v[:,ax]                 # get just the sorted column, to search
    s = bisect(z, a)            # find the start
    e = bisect(z, b)            # and the end
    return v[s:e+1]             # return them

class view(object):
    """v = view(ax) where ax is in 'ASC'"""
    def __init__(self, ax):
        if ax == 'A':   # axial view
            self.axes = ('LR', 'PA')
        elif ax == 'S': # sagittal view
            self.axes = ('AP', 'IS')
        elif ax == 'C': # coronal view
            self.axes = ('LR', 'IS')

# These are in LPI order
axes_ends = {
    'A': [0, 100, 0], 'P': [0, -100, 0],
    'S': [0, 0, 140], 'I': [0, 0, -20],
    'L': [-80, 0, 0], 'R': [80, 0, 0],
}

def get_axis(LPIs, LPIe):
    LPIs = resize(LPIs, (4, 1))
    LPIs[-1] = 1.
    LPIe = resize(LPIe, (4, 1))
    LPIe[-1] = 1.
    s = dot(fromLPI, LPIs)
    e = dot(fromLPI, LPIe)
    return s[:,0], e[:,0]

a, p = get_axis(axes_ends['A'], axes_ends['P'])
i, s = get_axis(axes_ends['I'], axes_ends['S'])
a2p = range(int(a[0]), int(p[0]))
i2s = range(int(i[1]), int(s[1]), -1)
i2s = list(i2s)
i2s.reverse()

slice = 140

slicea = [0, 0, slice, 1]
sliceb = [0, 0, slice + 1, 1]
a = dot(toLPI, slicea)
b = dot(toLPI, sliceb)
h = find_verts('L', a[0], b[0])

v = h[:,:3] # verts
n = h[:,3:] # normal vectors
#v += n * 2.    # move the vert along the normal
e = v + n * 2.  # the end is 2 normals away
v /= 10.        # now cm
e /= 10.
r = range(len(v))

# data: [A->P, S->I, L->R] a.k.a. ASL

d = data[:, :, slice]
d = d[:, i2s]
d = d[a2p, :]

#from numpy import random
#colors = []
#for i in r:
#    colors.append(random.random(4))

imshow(d.T, extent = (10., -10., -2., 14.), cmap = cm.bone)
##plot(v[:,1], v[:,2], 'y.')
##plot(v[:,1], v[:,2], 'y.', markersize = 10)

ll = []
for i in r:
    l = [(v[i][1], v[i][2]), (e[i][1], e[i][2])]
    ll.append(l)
c = LineCollection(ll)
#c.set_color(colors)
c.set_color('y')
c.set_linewidth(2.)
ax = gca()
ax.add_collection(c)

show()
