#! /usr/bin/env python

import sys, os.path, getopt, re
from numpy import array, hypot

# Add the SAMsrcV5 lib directory to the Python path,
# so we can find thd_atr.py and the others.

sys.path.append("@@libdir@@")
from thd_atr import afni_header_read
from fid import *
from samutil import *

usage("""[-m] name
Return the interfiducial distances for name.
If name is an AFNI brik, the tagset is examined.
If name is an MEG dataset, the .hc file is examined (the head coil locations),
if -m is used, the .hdm file is examined (the output of localSpheres).""")

HDM = 1
HC = 2
BRIK = 3
mflg = None
type = None

optlist, args = parseargs("m")

for opt, arg in optlist:
    if opt == '-m':
        mflg = True

if len(args) != 1:
    printusage()
    sys.exit(1)

name = args[0]
if name[-1] == '/':
    name = name[:-1]    # remove trailing slash
dirname = os.path.expanduser(os.path.dirname(name))
filename = os.path.basename(name)

# If the argument is a .ds directory, get the corresponding .hc or .hdm file.

base, ext = os.path.splitext(filename)
if ext == '.ds':
    if mflg:
        s = "default.hdm"
        type = HDM
    else:
        s = "%s.hc" % base
        type = HC
    msg("using %s\n" % s)
    filename = os.path.join(dirname, filename, s)
    x = open(filename)

# If it's not an MEG dataset, assume it's a BRIK and read the header.

else:
    filename = os.path.join(dirname, filename)
    h = afni_header_read(filename)
    if not h.get('TAGSET_NUM'):
        printerror("%s has no tags!" % filename)
        sys.exit(1)
    type = BRIK

def get_coord(f):
    l = next(f)
    return float(l.split()[2])

# HC
def coord(f):
    x = get_coord(f)
    y = get_coord(f)
    z = get_coord(f)
    return array((x, y, z))

# HDM
def coord2(s):
    return array(list(map(int, s.split()[-3:]))) * .1

# BRIK
def coord3(s):
    x = array(list(map(fuzz, tl))[0:3]) * .1
    # convert from RAI to PRI
    return array((-x[1], x[0], x[2]))

def fuzz(t):
    if abs(t) < 1.e-8:
        return 0.
    return t

if type == HC:
    nasion = re.compile('measured nasion .* head')
    left = re.compile('measured left .* head')
    right = re.compile('measured right .* head')

    for s in x:
        if nasion.match(s):
            n = coord(x)
        elif left.match(s):
            l = coord(x)
        elif right.match(s):
            r = coord(x)

elif type == HDM:
    nasion = re.compile('.*NASION:')
    left = re.compile('.*LEFT_EAR:')
    right = re.compile('.*RIGHT_EAR:')
    sres = re.compile('.*SAGITTAL:')
    cres = re.compile('.*CORONAL:')
    ares = re.compile('.*AXIAL:')
    sr = None

    for s in x:
        if nasion.match(s):
            n = coord2(s)
        if left.match(s):
            l = coord2(s)
        if right.match(s):
            r = coord2(s)
        if sres.match(s):
            sr = float(s.split()[-1])
        if cres.match(s):
            cr = float(s.split()[-1])
        if ares.match(s):
            ar = float(s.split()[-1])
    if sr == None:
        msg("no resolution found, assuming 1 mm/voxel\n")
        sr = cr = ar = 1.
    res = array((sr, cr, ar))
    n *= res
    l *= res
    r *= res
    msg("using slice coordinates: SCA\n")

elif type == BRIK:
    ntags, pertag = h['TAGSET_NUM']
    f = h['TAGSET_FLOATS']
    lab = h['TAGSET_LABELS']
    d = {}
    for i in range(ntags):
        tl = f[i * pertag : (i+1) * pertag]
        d[lab[i]] = coord3(tl)
    n = d[NASION]
    l = d[LEAR]
    r = d[REAR]

def length(d):
    return hypot.reduce(d)

print('nasion: %.3f %.3f %.3f' % tuple(n))
print('left ear: %.3f %.3f %.3f' % tuple(l))
print('right ear: %.3f %.3f %.3f' % tuple(r))
print('left - right: %.3f cm' % length(l - r))
print('nasion - left: %.3f cm' % length(n - l))
print('nasion - right: %.3f cm' % length(n - r))

if l[1] < 0:
    msg("Warning: left / right flip detected.\n")
