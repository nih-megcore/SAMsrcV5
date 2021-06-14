#! /usr/bin/env python

import sys, os
import numpy as np
import nibabel
from nibabel.freesurfer.io import read_geometry
from nibabel.gifti import GiftiDataArray

codes = nibabel.nifti1.intent_codes
INDEX = codes['NIFTI_INTENT_NODE_INDEX']
DATA = codes['NIFTI_INTENT_NONE']

try:
    FShome = os.environ['FREESURFER_HOME']
except KeyError:
    printerror("You must set the FREESURFER_HOME environment variable.")
    sys.exit(1)

try:
    Subjdir = os.environ['SUBJECTS_DIR']
except KeyError:
    Subjdir = os.path.join(FShome, "subjects")
    print("Note: Using the default SUBJECTS_DIR: {}".format(Subjdir))

if len(sys.argv) != 2:
    print("Usage: {} giifile".format(sys.argv[0]))
    sys.exit(1)

giiname = sys.argv[1]

if ".lh." in giiname:
    hemi = "lh"
elif ".rh." in giiname:
    hemi = "rh"
else:
    print("Can't determine hemisphere!")
    sys.exit(1)

g = nibabel.load(giiname)
subjid = g.meta.metadata['SubjectID']
surf = g.meta.metadata['SurfaceID']

prefix = "/tmp/moo"
outfile = "{}.{}.gii".format(prefix, hemi)

dir = os.path.join(Subjdir, subjid, "surf")
geomname = os.path.join(dir, "{}.{}".format(hemi, surf))
geom = read_geometry(geomname)
nv = len(geom[0])
v = np.zeros((nv,), 'f')

for da in g.darrays:
    if da.intent == INDEX:
        idx = da.data
        break

for da in g.darrays:
    if da.intent == DATA:
        for i, f in zip(idx, da.data):
            v[i] = f

dalist = []
dalist.append(GiftiDataArray(v,
                             intent = 'NIFTI_INTENT_NONE',
                             datatype = 'NIFTI_TYPE_FLOAT32'))

gim = nibabel.gifti.gifti.GiftiImage(meta = g.meta, darrays = dalist)
nibabel.save(gim, outfile)

