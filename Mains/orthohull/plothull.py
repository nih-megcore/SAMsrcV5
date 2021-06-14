#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("usage: {} hull.shape".format(sys.argv[0]))
    sys.exit(1)

name = sys.argv[1]

def readhull(name):
    with open(name) as f:
        n = int(next(f))
        if n != 9002:
            print("warning, non-standard number of vertices")
        for i in range(n):
            yield next(f)

h = np.loadtxt(readhull(name))

fig, ax = plt.subplots(2, 2)

ax[0,0].plot(h[:,0], h[:,1], ',')
ax[0,0].set_xlabel("P <-> A")
ax[0,0].set_ylabel("R <-> L")

ax[1,0].plot(h[:,0], h[:,2], ',')
ax[1,0].set_xlabel("P <-> A")
ax[1,0].set_ylabel("I <-> S")

ax[0,1].plot(h[:,1], h[:,2], ',')
ax[0,1].set_xlabel("R <-> L")
ax[0,1].set_ylabel("I <-> S")

for a in [ax[0,0], ax[1,0], ax[0,1]]:
    a.set_aspect('equal')
    #a.set_aspect('equal', adjustable = 'box-forced')
    a.grid()

ax[1, 1].remove()

fig.text(.5, .95, name, ha = 'center', size = 'large')

fig.tight_layout(rect = (0, 0, 1, .95))

plt.show()
