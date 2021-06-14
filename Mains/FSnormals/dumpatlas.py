#! /usr/bin/env python

import sys
import nibabel
import numpy as np

V = np.empty((0,3))
N = np.empty((0,3))
g = nibabel.load(sys.argv[1])
for da in g.darrays:
    if da.intent == nibabel.nifti1.intent_codes['NIFTI_INTENT_POINTSET']:
        V = np.vstack((V, da.data))
    if da.intent == nibabel.nifti1.intent_codes['NIFTI_INTENT_VECTOR']:
        N = np.vstack((N, da.data))

for v, n in zip(V, N):
    print(v[0], v[1], v[2], n[0], n[1], n[2])
