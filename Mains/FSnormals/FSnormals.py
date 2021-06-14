#! /usr/bin/env python

from __future__ import print_function

import os, sys
import tempfile
import struct
import numpy as np
from scipy.spatial import KDTree

sys.path.append("@@libdir@@")
# these are from orthohull
from samutil import *
from thd_atr import afni_header_read

surf = "smoothwm"                   # default surface to use

valid_s = ["pial", "smoothwm", "white"]
surfmsg = "surface must be one of " + ', '.join(valid_s) + " (default %s)" % surf

usage("""[-q] [-a annotation] ... [-t ortho+orig] [-i fac] [-p parcfile]
    [-g "step radius"] [-r thresh] [-L|-R|-B] -o outputfile subject [surface]

Output annotated surface vertices and normal vectors for the FreeSurfer
subject 'subject'. The output has 7 columns: vx, vy, vz, nx, ny, nz, and
annot * hemi, where hemi is 1 for left and -1 for right. Coordinates are in
meters, in PRI order.

-q means quiet, default verbose.

-a restricts the output to the given annotation region. -a may be repeated;
   the value represents either left or right (or both) hemispheres.

-p specifies the parcellation file to use. Default is aparc.a2009s, which is
   the Destrieux atlas. Other options include aparc (the Desikan-Killiany
   atlas) and aparc.DKTatlas (the Desikan-Killiany-Tourville atlas).

-t rotates the output vertices and normals into the ortho+orig space;
   ortho+orig must be the output of 3dTagalign.

-i moves each vertex along its normal vector according to v' = v + fac * n
   (fac in mm, may be negative). Note that this conversion is done before
   resampling to a grid (if -g is used).

-g "step radius" restricts the output vertices to the grid defined by
   the ortho+orig MRI (you must use -t in this case). The step argument
   is an integer specifying the number of MRI voxels between output grid
   points. At each grid point, all the normal vectors within a ball of
   the given radius (mm) are averaged together to form the output
   normal. For example -g "3 1.5" means one grid point every 3 MRI
   voxels, and an averaging ball with a 3 mm (2 x 1.5) diameter.

-r thresh filters out vertices which have near radial normals. The
   threshold is the absolute value of the dot product of the vertex's
   position vector (relative to the centroid) with the normal vector, so
   ranges from 0 to 1. Large values allow more radial normals, small values
   (e.g. -r .2) pass only more tangential normals.

-L or -R means just output that hemisphere, -B (the default) means both.

-o specifies the output filename (this is a non-optional option).

%s

You must set the FREESURFER_HOME environment variable, and optionally,
SUBJECTS_DIR (default $FREESURFER_HOME/subjects).
""" % surfmsg)

verbose = True
annots = []
aparcfile = 'aparc.a2009s'
hemi = 'B'
ortho = None
gridStep = None
fac = 0
radfilt = 0
outfile = None

optlist, args = parseargs("qLRBa:p:t:g:i:r:o:")
for opt, arg in optlist:
    if opt == '-a':
        annots.append(int(arg))
    elif opt == '-p':
        aparcfile = arg
    elif opt == '-o':
        outfile = arg
    elif opt == '-R':
        hemi = 'R'
    elif opt == '-L':
        hemi = 'L'
    elif opt == '-q':
        verbose = False
    elif opt == '-t':
        ortho = arg
    elif opt == '-g':
        l = arg.split()
        if len(l) != 2:
            printerror("-g requires 2 arguments, step and radius.")
            sys.exit(1)
        gridStep = int(l[0])
        gridRadius = float(l[1])
    elif opt == '-i':
        fac = float(arg)
    elif opt == '-r':
        radfilt = float(arg)
        if radfilt == 0:
            radfilt = -1.

if len(args) < 1 or len(args) > 2:
    printusage()
    sys.exit(1)

if outfile is None:
    printerror("You must specify an output filename with -o.")
    sys.exit(1)

subj = args[0]
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

if gridStep and ortho is None:
    printerror("You must specify an ortho MRI with -t when using -g.")
    sys.exit(1)

# Get the transform to ortho space.

if ortho:
    h = afni_header_read(ortho)

    m = np.array(h['TAGALIGN_MATVEC'])  # RAI mm
    m.shape = (3, 4)

    rot = m[0:3, 0:3]   # rotation
    trans = m[0:3, 3]   # translation (mm)

    if gridStep:
        # We'll also need the intrinsic coordinate system of the ortho volume.

        toLPI = np.array(h['IJK_TO_DICOM_REAL'])    # RAI mm
        toLPI.shape = (3, 4)

        orthoShape = h['DATASET_DIMENSIONS'][:3]

# Read FreeSurfer annotation files.

def read_annot(afile, nv, ctab):
    f = open(afile, 'rb')

    # The first int is the number of elements.

    fmt = ">1i"
    h = f.read(struct.calcsize(fmt))
    l = struct.unpack(fmt, h)

    N = l[0]

    # Next come N pairs of ints. Vertex number and annotation.
    # The annotation is a color map value. Turn it into its index.

    fmt = ">%di" % (N * 2)
    b = f.read(struct.calcsize(fmt))
    l = struct.unpack(fmt, b)
    na = int(len(l) / 2)
    #print("na: %d, nv:%d" % (na,nv))
    annot = [0] * nv
    #print(na)
    for i in range(na):
        ii = i * 2
        v = l[ii]
        #print(i, end=',')
        #print("i: %d ii: %d v: %d" % (i, ii, v))
        if v < 0 or v >= nv:
            printerror("annotation vertex number out of range! (%d)" % v)
            cleanup()

        # This mess is here because in some of the parcellations the cmap value
        # for Unknown is 0 0 0 and in others it's 25 5 25, but the annot files
        # all reference 0.

        key = l[ii+1]
        if key == 0:
            if key in ctab:
                annot[v] = ctab[key]
            else:
                annot[v] = 0
        else:
            annot[v] = ctab[key]

    f.close()
    return annot

def read_annot_ctab(name):
    f = open(name)
    ctab = {}
    for l in f:
        s = l.split()
        i = int(s[0])
        a = int(s[4]) * 65536
        a += int(s[3]) * 256
        a += int(s[2])
        ctab[a] = i
    return ctab

# Create a temp directory to work in.

OrigDir = os.getcwd()
TempDir = tempfile.mkdtemp()
assert len(TempDir) != 0

def cleanup(normal = False):
    os.system("rm -f %s/info.txt" % TempDir)
    os.system("rm -f %s/*.asc" % TempDir)
    os.system("rmdir %s" % TempDir)
    if not normal:
        sys.exit(1)
    sys.exit(0)

# We need the origin offset of the original MRI ...

dir = os.path.join(Subjdir, subj, "mri", "orig")
try:
    os.chdir(dir)
except OSError:
    printerror("Error locating subject directory %s" % dir)
    cleanup()

try:
    cmd = os.path.join(FShome, "bin", "mri_info")
    os.system("%s --cras 001.mgz > %s/info.txt" % (cmd, TempDir))
    l = open("%s/info.txt" % TempDir).readline().split()
    origin = np.array(list(map(float, l)))
except:
    printerror("Can't find origin offset!")
    cleanup()

# Now switch to the dir with the surfaces in it.

dir = os.path.join(Subjdir, subj, "surf")
try:
    os.chdir(dir)
except OSError:
    printerror("Error locating subject directory %s" % dir)
    cleanup()

# Generate the .asc files for both hemispheres.

if verbose:
    msg("Generating %s surface meshes ...\n" % surf)
cmd = os.path.join(FShome, "bin", "mris_convert")
t = (cmd, surf, TempDir, surf)
os.system("%s rh.%s %s/rh.%s.asc > /dev/null" % t)
os.system("%s -n rh.%s %s/rh.%s.normals.asc > /dev/null" % t)
os.system("%s lh.%s %s/lh.%s.asc > /dev/null" % t)
os.system("%s -n lh.%s %s/lh.%s.normals.asc > /dev/null" % t)

# Read them.

try:
    rhv = open("%s/rh.%s.asc" % (TempDir, surf)).readlines()
    rhn = open("%s/rh.%s.normals.asc" % (TempDir, surf)).readlines()
    lhv = open("%s/lh.%s.asc" % (TempDir, surf)).readlines()
    lhn = open("%s/lh.%s.normals.asc" % (TempDir, surf)).readlines()
except IOError:
    printerror("There was a problem generating the .asc files.")
    cleanup()

rh_nv  = int(rhv[1].split()[0])
rhn_nv = int(rhn[1].split()[0])
lh_nv  = int(lhv[1].split()[0])
lhn_nv = int(lhn[1].split()[0])

if (rh_nv != rhn_nv) or (lh_nv != lhn_nv):
    printerror("Inconsistent numbers of vertices and normals!")
    cleanup()

# Now switch to the label dir.

dir = os.path.join(Subjdir, subj, "label")
try:
    os.chdir(dir)
except OSError:
    printerror("Error locating subject directory %s" % dir)
    cleanup()

# Read the cortex label files.

try:
    rhlab = open("rh.cortex.label").readlines()
    lhlab = open("lh.cortex.label").readlines()
    if verbose:
        msg("Restricting to cortex ...\n")
    ctxmap_rh = [False] * rh_nv
    ctxmap_lh = [False] * lh_nv
    for s in rhlab[2:]:
        i = int(s[:s.index(' ')])       # convert a space delimited int
        ctxmap_rh[i] = True             # mark this vert as being cortex
    for s in lhlab[2:]:
        i = int(s[:s.index(' ')])
        ctxmap_lh[i] = True

except IOError:
    printerror("Warning: there was a problem reading the label files.")
    ctxmap_rh = [True] * rh_nv
    ctxmap_lh = [True] * lh_nv

# Read the annotation files.

aparclist = aparcfile.split('.')
if len(aparclist) == 1 and aparclist[0] == 'aparc':
    ctabfile = "aparc.annot.ctab"
    rhannot = "rh.aparc.annot"
    lhannot = "lh.aparc.annot"
else:
    ctabfile = "aparc.annot.%s.ctab" % (aparclist[1])
    rhannot = "rh.aparc.%s.annot" % (aparclist[1])
    lhannot = "lh.aparc.%s.annot" % (aparclist[1])

try:
    ctab = read_annot_ctab(ctabfile)
    #print(rhannot, rh_nv, ctab)
    rhannot = read_annot(rhannot, rh_nv, ctab)
    lhannot = read_annot(lhannot, lh_nv, ctab)
except IOError:
    printerror("Warning: there was a problem reading the annotation files.")

if len(annots) > 0:
    if verbose:
        msg("Restricting to given annotation labels ...\n")
    for i in range(rh_nv):
        if rhannot[i] not in annots:
            ctxmap_rh[i] = False
    for i in range(lh_nv):
        if lhannot[i] not in annots:
            ctxmap_lh[i] = False

def counth(nv, ctxmap):
    c = 0
    for i in range(nv):
        if not ctxmap[i]:
            continue
        c += 1
    return c

def norm(v):
    return v / np.sqrt((v * v).sum())

# Convert each vertex and normal.

v = np.zeros(3)
n = np.zeros(3)
cent = np.array([0., 0., 50.]) # default centroid in LPI mm

def mkLPIsurf(nv, hv, hn, origin, ctxmap, annotmap, hemi):
    global v, n

    surf = []

    for i in range(nv):
        if not ctxmap[i]:
            continue

        # Get the vert and the normal in LPI (ras) order (mm).

        v[:] = [float(x) for x in hv[i+2].split()[:3]]
        n[:] = [float(x) for x in hn[i+2].split()[:3]]

        # Translate back to the original origin (origin is LPI mm).

        v += origin

        # Optionally rotate.

        if ortho:
            # convert from LPI mm to RAI mm
            v[:] = -v[0], -v[1], v[2]       # @@@ could modify rot and trans
            n[:] = -n[0], -n[1], n[2]
            # transform
            v[:] = rot.dot(v) + trans
            n[:] = rot.dot(n)
            # convert from RAI mm to LPI mm
            v[:] = -v[0], -v[1], v[2]
            n[:] = -n[0], -n[1], n[2]

        # Optionally move the vertex up along the normal.

        if fac:
            v[:] = v + fac * n

        # If the normal is too radial, optionally ignore this vertex.

        if radfilt:
            x = norm(v - cent).dot(n)
            if abs(x) > radfilt:
                continue

        surf.append((v[0], v[1], v[2], n[0], n[1], n[2], annotmap[i] * hemi))

    return surf

crh = counth(rh_nv, ctxmap_rh)
clh = counth(lh_nv, ctxmap_lh)
if verbose:
    if hemi in ['B', 'R']:
        msg("Note: %d rh verts\n" % crh)
    if hemi in ['B', 'L']:
        msg("Note: %d lh verts\n" % clh)
    if ortho:
        msg("Rotating into %s space.\n" % ortho)
    if fac:
        msg("Inflating verts by %g\n" % fac)
    if gridStep:
        msg("Resampling to grid.\n")

# Create LPI mm surfaces.

rhsurf = lhsurf = None

if hemi in ['B', 'R']:
    rhsurf = mkLPIsurf(rh_nv, rhv, rhn, origin, ctxmap_rh, rhannot, -1)
if hemi in ['B', 'L']:
    lhsurf = mkLPIsurf(lh_nv, lhv, lhn, origin, ctxmap_lh, lhannot, 1)

# Combine into one.

if hemi in ['B', 'R']:
    surf = np.array(rhsurf)
if hemi in ['B', 'L']:
    if hemi == 'L':
        surf = np.array(lhsurf)
    else:
        surf = np.vstack((surf, np.array(lhsurf)))

# Optionally (but this is the usual case) restrict to a grid.

if gridStep:

    # Make a KD tree of just the verts.

    k = KDTree(surf[:,0:3])

    # Create a grid which we will transform into LPI space.

    Nx, Ny, Nz = orthoShape
    Rx = range(0, Nx, gridStep)
    Ry = range(0, Ny, gridStep)
    Rz = range(0, Nz, gridStep)

    # For each gridpoint, convert to LPI and find the nearest
    # vertices. Average those normals together.

    out = []
    for x in Rx:
        for y in Ry:
            for z in Rz:
                pos = [x, y, z, 1]
                a = np.dot(toLPI, pos)
                idx = k.query_ball_point(a[0:3], gridRadius)
                if len(idx) > 0:
                    n = surf[idx][:,3:6].mean(axis = 0)
                    v = a[0:3]
                    l = [v[0], v[1], v[2], n[0], n[1], n[2], surf[idx[0]][6]]
                    out.append(l)

    surf = np.array(out)

# Convert to PRI m and write everything out.

os.chdir(OrigDir)
f = open(outfile, "w")

nv = surf.shape[0]
for i in range(nv):

    # Convert to meters, and normalize the normal.

    v = surf[i, 0:3] / 1000.
    n = surf[i, 3:6]
    n = n / np.sqrt((n * n).sum())

    # Write in PRI order.

    print(v[1], -v[0], v[2], n[1], -n[0], n[2], int(surf[i, 6]), file = f)

f.close()

if verbose:
    msg("Wrote %d vertices.\n" % nv)

cleanup(True)
