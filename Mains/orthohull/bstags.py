#! /usr/bin/env python3

import os, sys

if len(sys.argv) != 2:
    print("""usage: {} file.txt
Where file.txt is the saved electrode location output from Brainsight.""".format(sys.argv[0]))
    sys.exit(1)

# The tags have to be output in this order ...

fiducials0 = [ "Nasion", "Left Ear", "Right Ear" ]
fiducials = [ "Nasion", "Left Ear", "Right Ear", "LPA", "RPA" ]
alt_fid = { "LPA" : "Left Ear", "RPA" : "Right Ear" }

# Write blah.tag if the input was blah.txt, otherwise just append .tag

txtname = sys.argv[1]
name, ext = os.path.splitext(txtname)
if ext != ".txt":
    name = txtname
tagname = "{}.tag".format(name)

# Check that the file has the proper coordinates.

scanner = False
ll = open(txtname).readlines()
for l in ll:
    l = l.strip()
    if l[0] == '#':
        if l == "# Coordinate system: NIfTI:Scanner":
            scanner = True

if not scanner:
    print("{} does not appear to be in NIfTI:Scanner coordinates".format(txtname))

# Find the tags and write them in the correct format for AFNI.

tags = {}
for l in ll:
    l = l.strip()
    if l[0] != '#':
        l = l.split('\t')
        if l[0] in fiducials:
            x, y, z = [float(l[i]) for i in [3, 4, 5]]
            if l[0] in alt_fid:
                l[0] = alt_fid[l[0]]
            # the coordinates are LPI, convert them to RAI
            tags[l[0]] = "'{}' {} {} {}".format(l[0], -x, -y, z)

print("writing {}".format(tagname))
f = open(tagname, 'w')
for tag in fiducials0:
    print(tags[tag], file = f)
f.close()
