#! /usr/bin/env python

import sys
import numpy as np

plyHeader = """ply
format ascii 1.0
comment author: Greg Turk
element vertex {nverts}
property float x
property float y
property float z
property float nx
property float ny
property float nz
element face {ntri}
property list uchar int vertex_indices
end_header"""

if len(sys.argv) != 2:
    print("usage: {} hull.shape".format(sys.argv[0]))
    sys.exit(1)

name = sys.argv[1]

shapeFile = open(name)
nverts = int(shapeFile.readline())
verts = np.zeros((nverts, 6))
for i in range(nverts):
    s = shapeFile.readline()
    l = [float(x) for x in s.split()]
    verts[i,:] = l

ntri = int(shapeFile.readline())

print(plyHeader.format(nverts=nverts, ntri=ntri))
for i in range(nverts):
    print(' '.join(str(x) for x in verts[i,:]))
for i in range(ntri):
    s = shapeFile.readline()
    print("3", s, end='')

