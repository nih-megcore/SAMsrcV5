#! /usr/bin/env python

from __future__ import print_function

import sys, os, signal
import tempfile
import numpy as np
import nibabel
from nibabel.gifti import GiftiDataArray

sys.path.append("@@libdir@@")
from samutil import *

surf = "smoothwm"                   # default surface to use

valid_s = ["pial", "smoothwm", "white"]
surfmsg = "surface must be one of " + ', '.join(valid_s) + " (default %s)" % surf

usage("""[-q] [-t ortho+orig] [-a annotation] ... [-p parcfile] [-s sfrac]
    [-i fac] [-r thresh] [-L|-R|-B] -o outputname subject [surface]

Output annotated surface vertices and normal vectors for the FreeSurfer
subject 'subject'. The output is in GIFTI format; ".gii" is appended to
the outputname if it doesn't already have an extension.
Coordinates are in meters, in PRI order.

-q means quiet, default verbose.

-t rotates the output vertices and normals into the ortho+orig space;
   ortho+orig must be the output of 3dTagalign.

-a restricts the output to the given annotation region. -a may be repeated;
   the value represents either left or right (or both) hemispheres.

-p specifies the parcellation file to use. Default is aparc.a2009s, which is
   the Destrieux atlas. Other options include aparc (the Desikan-Killiany
   atlas) and aparc.DKTatlas (the Desikan-Killiany-Tourville atlas).

-s sfrac means to randomly sample the surface, rather than using every vertex.
   sfrac must be between 0 and 1; a smaller number means fewer vertices.

-i moves each vertex along its normal vector according to v' = v + fac * n
   (fac in mm, may be negative).

-r thresh filters out vertices which have near radial normals. The
   threshold is the absolute value of the dot product of the vertex's
   position vector (relative to the centroid) with the normal vector, so
   ranges from 0 to 1. Large values allow more radial normals, small values
   (e.g. -r .2) pass only more tangential normals.

-L or -R means just output that hemisphere, -B (the default) means both.

-o specifies the output filename (this is a non-optional option). If the
   name does not end in ".gii", it is added.

%s

You must set the FREESURFER_HOME environment variable, and optionally,
SUBJECTS_DIR (default $FREESURFER_HOME/subjects).
""" % surfmsg)

verbose = True
annots = []
aparcfile = 'aparc.a2009s'
hemi = ['rh', 'lh']     # We always write the RH first!
ortho = None
sfrac = 1
fac = 0
radfilt = 0
outfile = None

optlist, args = parseargs("qLRBa:p:t:s:i:r:o:")
for opt, arg in optlist:
    if opt == '-q':
        verbose = False
    elif opt == '-t':
        ortho = arg
    elif opt == '-a':
        annots.append(int(arg))
    elif opt == '-p':
        aparcfile = arg
    elif opt == '-s':
        sfrac = float(arg)
        if sfrac <= 0. or sfrac > 1.:
            printerror("sfrac must be between 0 and 1.")
            sys.exit(1)
    elif opt == '-i':
        fac = float(arg) / 1000.    # convert mm to m
    elif opt == '-r':
        radfilt = float(arg)
        if radfilt <= 0. or radfilt > 1.:
            printerror("radius threshold must be between 0 and 1.")
            sys.exit(1)
    elif opt == '-o':
        outfile = arg
        if outfile[-4:] != ".gii":
            outfile += ".gii"
    elif opt == '-R':
        hemi = ['rh']
    elif opt == '-L':
        hemi = ['lh']

if len(args) < 1 or len(args) > 2:
    printusage()
    sys.exit(1)

if outfile is None:
    printerror("You must specify an output filename with -o.")
    sys.exit(1)

subjid = args[0]
if len(args) > 1:
    surf = args[1]

if surf not in valid_s:
    printerror(surfmsg)
    sys.exit(1)

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

# Get the transform to ortho space.

if ortho:
    # This is in RAI mm

    s = run("3dAttribute TAGALIGN_MATVEC {}".format(ortho), True)
    m = np.array([float(x) for x in s.split()])
    m.shape = (3, 4)

    rot = m[0:3, 0:3]   # rotation
    trans = m[0:3, 3]   # translation (mm)

# Create a temp directory to work in.

TempDir = tempfile.mkdtemp()
assert len(TempDir) != 0

def cleanup(normal = False):
    os.system("rm -rf {}".format(TempDir))
    if not normal:
        sys.exit(1)
    sys.exit(0)

def handleINT(sig, fr):
    cleanup()

signal.signal(signal.SIGINT, handleINT)

# We need the origin offset of the original MRI ...

try:
    cmd = os.path.join(FShome, "bin", "mri_info")
    name = os.path.join(Subjdir, subjid, "mri", "orig", "001.mgz")
    s = run("{} --cras {}".format(cmd, name), True)
    origin = np.array([float(x) for x in s.split()])
except:
    printerror("Can't find origin offset!")
    cleanup()

# Convert strings of floats to an array.

def cvt_to_float(l):
    n = len(l)
    r = np.empty((n, 3), 'f')
    for i in range(n):
        r[i,:] = [float(x) for x in l[i].split()[:3]]
    return r

# Now get the surfaces (with normals) for both hemispheres as .asc files.

if verbose:
    msg("Generating %s surface meshes.\n" % surf)
dir = os.path.join(Subjdir, subjid, "surf")

d = {}
vert = {}
norm = {}
nv = {}
d['convert'] = os.path.join(FShome, "bin", "mris_convert")
for h in hemi:
    d['insurf'] = os.path.join(dir, "{}.{}".format(h, surf))
    d['asc'] = os.path.join(TempDir, "{}.{}.asc".format(h, surf))
    d['norm'] = os.path.join(TempDir, "{}.{}.normals.asc".format(h, surf))

    os.system("{convert} {insurf} {asc} > /dev/null".format(**d))
    os.system("{convert} -n {insurf} {norm} > /dev/null".format(**d))

    # Read them.

    try:
        v = open("{asc}".format(**d)).readlines()
        n = open("{norm}".format(**d)).readlines()
    except IOError:
        printerror("There was a problem generating the .asc files.")
        cleanup()

    nv[h] = int(v[1].split()[0])

    vert[h] = cvt_to_float(v[2:nv[h]+2])
    norm[h] = cvt_to_float(n[2:nv[h]+2])

# Read the cortex label files.

labeldir = os.path.join(Subjdir, subjid, "label")

ctxlabel = {}
for h in hemi:
    name = os.path.join(labeldir, "{}.cortex.label".format(h))
    l = nibabel.freesurfer.io.read_label(name)
    a = np.zeros(nv[h], 'i')
    for i in l:
        a[i] = 1
    ctxlabel[h] = a

# Read the annotation files.

if len(annots) > 0:
    if verbose:
        msg("Restricting to given annotation labels.\n")
    if aparcfile[-6:] != ".annot":
        aparcfile += ".annot"
    for h in hemi:
        name = os.path.join(labeldir, "{}.{}".format(h, aparcfile))
        try:
            annot = nibabel.freesurfer.io.read_annot(name)[0]
        except:
            printerror("There was a problem reading the annotation file {}".format(name))
            cleanup()
        a = ctxlabel[h]
        for i in range(nv[h]):
            if annot[i] not in annots:
                a[i] = 0

# Convert from LPI mm to (ortho) PRI m.

if verbose:
    msg("Converting to PRI coordinates.\n")
    if ortho:
        msg("Rotating to ortho space.\n")

for h in hemi:
    vs = vert[h]
    ns = norm[h]

    for i in range(nv[h]):
        v = vs[i]
        n = ns[i]

        # Translate back to the original origin (origin is LPI mm).

        v += origin

        # Convert from LPI to RAI.

        v[:] = -v[0], -v[1], v[2]
        n[:] = -n[0], -n[1], n[2]

        # Optionally rotate.

        if ortho:
            v[:] = rot.dot(v) + trans
            n[:] = rot.dot(n)

        # Convert from RAI mm to PRI m.

        v /= 1000.
        v[:] = -v[1], v[0], v[2]
        n[:] = -n[1], n[0], n[2]

        # Optionally move the vertex up along the normal.

        if fac:
            v[:] = v + fac * n

# Compute bounding box and centroid.

bbox = np.zeros((3, 2))
for h in hemi:
    for i in range(3):
        b = [vert[h][:,i].min(), vert[h][:,i].max()]
        bbox[i,:] = [min(bbox[i,0], b[0]), max(bbox[i,1], b[1])]
cent = np.zeros((3,))
for h in hemi:
    cent += vert[h].mean(axis = 0)
cent /= len(hemi)
bbox *= 100
if verbose:
    msg("bbox (cm):\n    {}\n    {}\n".format(bbox[:,0], bbox[:,1]))
    msg("centroid (cm):\n    {}\n".format(cent * 100))

# Optionally restrict based on orientation.

def normalize(v):
    return v / np.sqrt((v * v).sum())

if radfilt:
    for h in hemi:
        v = vert[h]
        n = norm[h]
        a = ctxlabel[h]

        # If the normal is too radial, ignore this vertex.

        for i in range(nv[h]):
            x = normalize(v[i] - cent).dot(n[i])
            if abs(x) > radfilt:
                a[i] = 0

# Further restrict by taking a random sample of vertices,
# and create the node index.

idx = {}
for h in hemi:
    a = ctxlabel[h]
    n = a.shape[0]
    i = np.arange(n)
    np.random.shuffle(i)
    ndel = int((1. - sfrac) * n + .5)
    a[i[:ndel]] = 0
    idx[h] = np.where(a == 1)[0]

# Write a GIFTI file.

if verbose:
    msg("Writing output GIFTI file.\n")
    s = ""
    for h in hemi:
        s += "{} {} vertices\n".format(len(idx[h]), h.upper())
    msg(s)

dalist = []
meta = {}
for h in hemi:
    # GIFTI standard metadata
    meta['AnatomicalStructurePrimary'] = 'CortexLeft' if h == 'lh' else 'CortexRight'
    meta['GeometricType'] = 'Anatomical'
    md = nibabel.gifti.GiftiMetaData().from_dict(meta)

    dalist.append(GiftiDataArray(idx[h], meta = md,
                                 encoding = 'GIFTI_ENCODING_ASCII',
                                 intent = 'NIFTI_INTENT_NODE_INDEX',
                                 datatype = 'NIFTI_TYPE_INT32'))
    dalist.append(GiftiDataArray(vert[h][idx[h]], meta = md,
                                 intent = 'NIFTI_INTENT_POINTSET',
                                 datatype = 'NIFTI_TYPE_FLOAT32'))
    dalist.append(GiftiDataArray(norm[h][idx[h]], meta = md,
                                 intent = 'NIFTI_INTENT_VECTOR',
                                 datatype = 'NIFTI_TYPE_FLOAT32'))

meta = { 'SubjectID': subjid, 'SurfaceID': surf,
         'centroid': ' '.join([str(x) for x in cent]),
         'atlascmd': ' '.join(sys.argv) }
md = nibabel.gifti.GiftiMetaData().from_dict(meta)

gim = nibabel.gifti.gifti.GiftiImage(meta = md, darrays = dalist)
nibabel.save(gim, outfile)

cleanup(True)
