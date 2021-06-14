#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:14:05 2019

@author: nugenta
"""
from __future__ import print_function
import subprocess
import sys

sys.path.append("@@libdir@@")
master = "@@libdir@@/master+orig"
from samutil import *
from runscript import runcmd
from runscript import runscript

usage("""[-p surfacename] [-x prefix] [-i <inf>] [-t] [-m] [-q] [-c] [-o] [anat+orig])

Automated script for MRI processing. This script takes the MRI anatomical
image with marked fiducial points, performs skull stripping and Talairaching,
and computes the transform from the ortho space (plane of the fiducials) to
Talairach space and back. Also generates the hull.shape file corresponding
to the brain surface, with or without inflation. The original anatomy image
defaults to anat+orig, or it may be specified on the command line.

-p       Specifies which surface to use, the default is brainhull.ply.
         Other options are innerskull.ply, or outerskull.ply. See
         3dSkullStrip -help for more information.
-x       Used to specify a prefix for the output files; defaults to none
-i <inf> Inflates the surface by a given factor, defaults to -i.005, which
         inflates the surface by 0.5cm, preventing clipping. Use -i0 for no
         inflation
-t       Do the Talaraich transformation, using @auto_tlrc. The transform is
         stored by default in brain_tlrc
-m       Same as -t, but uses the MNI template instead
-o       By default, the Talairach transform is performed on the ortho brain
         image. Sometimes this fails. This option requests that the tranform
         is derived from the original brain image, producing an ORTHO_to_TLRC.1D
         transform file that can be applied to the output MEG images.
-q       Forces the output hull to be convex (for problem surfaces only)
-c       Cleanup, just delete all output files if present and exit. Note that if
         you want to clean up files generated with -x and -o, you'll need to
         specify those.
""")

doclean=0
doqhull=0
dotlrc=0
domni=0
doinf=0.005
doxformorig=0
prefix=''
prefix_rscr=''
ply='brainhull.ply'

optlist, args = parseargs("cqmtop:i:x:")
for opt, arg in optlist:
    if opt == '-c':
        doclean=1
    elif opt == '-q':
        doqhull=1;
    elif opt == '-m':
        domni=1
        dotlrc=1
        base='MNI_caez_N27+tlrc'
    elif opt == '-t':
        dotlrc=1
        base='TT_N27+tlrc'
    elif opt == '-i':
        doinf=float(arg)
    elif opt == '-x':
        prefix_no_undrscr=arg
        prefix = "%s_" % arg
    elif opt == '-p':
        ply = arg
    elif opt == '-o':
        doxformorig = 1

if doclean==1:
    if prefix=='':
        cmd = 'rm -f ortho+* mask+orig.* brain+*'
        subprocess.call(cmd,shell=True);
        cmd = 'rm -f anat*.ply ortho*.ply smooth*ply CreateIco.*'
        subprocess.call(cmd,shell=True);
        cmd = 'rm -f hull.shape multisphere.shape* tagalign_xform.1D'
        subprocess.call(cmd,shell=True);
        if doxformorig==0:
          cmd = 'rm -f ortho_WarpDrive.log ortho.Xaff12.1D ortho.Xat.1D'
        else:
          cmd = 'rm -f brain_WarpDrive.log brain.Xaff12.1D brain.Xat.1D'
        subprocess.call(cmd,shell=True);
    else:
        cmd = 'rm -f %sortho+* %smask+orig.* %sbrain+*' % (prefix, prefix, prefix)
        subprocess.call(cmd,shell=True);
        cmd = 'rm -f %s*.ply CreateIco.* ' % (prefix_no_undrscr)
        subprocess.call(cmd,shell=True);
        cmd = 'rm -f %shull.shape %s.shape* %stagalign_xform' % (prefix, prefix_no_undrscr,prefix)
        subprocess.call(cmd,shell=True);
        if doxformorig==0:
          cmd = 'rm -f %sortho_WarpDrive.log %sortho.Xaff12.1D %sortho.Xat.1D' % (prefix, prefix, prefix)
        else:
          cmd = 'rm -f %sbrain_WarpDrive.log %sbrain.Xaff12.1D %sbrain.Xat.1D' % (prefix, prefix, prefix)
        subprocess.call(cmd,shell=True);
    sys.exit('Done cleaning up')


if len(args) < 1 or len(args) > 2:
    printusage()
    sys.exit(1)

image = args[0]

# defining names for 3dSkullStrip

mask_outfileprefix = "%smask" % prefix
mask_outfile = "%smask+orig" % prefix
if prefix=='':
  ply_outfile="anat_%s" % (ply)
else:
  ply_outfile = "%s%s" % (prefix,ply)
brain_outfileprefix = "%sbrain" % prefix
brain_outfile = "%sbrain+orig" % prefix
brain_outfile_tlrc = "%sbrain+tlrc" % prefix

# Skullstrip the anatomical file

print("Skull Stripping the Anatomical Image")
if prefix=='':
  cmd = "3dSkullStrip -overwrite -input %s -prefix %s -mask_vol -skulls -o_ply anat" % (image,mask_outfileprefix)
else:
  cmd = "3dSkullStrip -overwrite -input %s -prefix %s -mask_vol -skulls -o_ply %s" % (image,mask_outfileprefix,prefix_no_undrscr)
subprocess.call(cmd,shell=True);
print("Masking the orig image")
cmd = """3dcalc -overwrite -a %s -b %s -expr 'b*step(a-2.9)' -prefix %s""" % (mask_outfile,image,brain_outfileprefix)
print(cmd)
subprocess.call(cmd,shell=True);

# copy tagset information

print("Copying tagset info from anat to brain");
cmd= "3drefit -overwrite -saveatr -atrcopy %s TAGSET_NUM -atrcopy %s TAGSET_FLOATS -atrcopy %s TAGSET_LABELS %s" % (image, image, image, brain_outfile)
subprocess.call(cmd,shell=True);

# defining names for 3dTagAlign

brain_ortho_outfileprefix = "%sortho" % (prefix)
brain_ortho_outfile = "%sortho+orig" % (prefix)
orthoply_outfile = "%sortho_%s" % (prefix,ply)
tagalign_xform = "%stagalign_xform.1D" % (prefix)

# Align to fiducial plane

print("Aligning the stripped brain file and surface ply to the fiducial plane")
cmd="3dTagalign -overwrite -prefix %s -master %s -matvec %s %s" % (brain_ortho_outfileprefix,master,tagalign_xform,brain_outfile)
subprocess.call(cmd,shell=True);
cmd="ConvertSurface -overwrite -i %s -o %s -xmat_1D %s" % (ply_outfile,orthoply_outfile,tagalign_xform)
subprocess.call(cmd,shell=True);

# defining names for SpharmDeco meshnorm and ply2fid

smooth_orthoply_prefix = "%ssmooth_ortho_%s" % (prefix,ply)
smooth_orthoply = "%ssmooth_ortho_%s" % (prefix,ply)
hull_shape_filename = "%shull.shape" % (prefix)
hull_name = "%shull" % (prefix)

# smooth the hull and make headmodel files

print("Smoothing the hull and running meshnorm")
cmd="CreateIcosahedron -rad 1 -ld 30 -tosphere"
subprocess.call(cmd,shell=True);
cmd="SpharmDeco -overwrite -i CreateIco.asc -unit_sph CreateIco -i %s -l 6 -o_ply %s" % (orthoply_outfile,smooth_orthoply_prefix)
subprocess.call(cmd,shell=True);
if doqhull == 1:
  cmd="meshnorm -i %s -q %s" % (doinf,smooth_orthoply_prefix)
  fp=open(hull_shape_filename,"w")
  subprocess.call(cmd,shell=True,stdout=fp);
else:
  cmd="meshnorm -i %s %s" % (doinf,smooth_orthoply_prefix)
  fp=open(hull_shape_filename,"w")
  subprocess.call(cmd,shell=True,stdout=fp);
if prefix=='':
    prefix_no_undrscr='multisphere'
cmd="ply2fid %s %s %s" % (brain_ortho_outfile,smooth_orthoply,prefix_no_undrscr)
subprocess.call(cmd,shell=True);

cmd="rm CreateIco.asc CreateIco.spec"
subprocess.call(cmd,shell=True)

# if Talairaching requested...

if dotlrc==1:

  if doxformorig==0:

    print("Talairaching the ortho/MEG space brain image")
    cmd="\@auto_tlrc -overwrite -no_ss -base %s -input %s -suffix NONE -pad_base 60" % (base,brain_ortho_outfile)
    subprocess.call(cmd,shell=True)
    print("Talairach transform stored in header of ortho image")

  elif doxformorig==1:

    print("Talairaching the orig space brain image, and computing ORTH_to_TLRC transform")
    cmd="\@auto_tlrc -overwrite -no_ss -base %s -input %s -suffix NONE" % (base,brain_outfile)
    subprocess.call(cmd,shell=True)

    # define names for transform files

    talairach_xform = "%sbrain.Xat.1D" % (prefix)
    talairach_xform_inv = "%sbrain_inv.Xat.1D" % (prefix)
    ORTH_to_TLRC = "%sORTHO_to_TLRC.1D" % (prefix)

    # make the MEG to TLRC transform

    fp=open(talairach_xform_inv,'w')
    cmd="cat_matvec %s -I" % (talairach_xform)
    subprocess.call(cmd,shell=True,stdout=fp);
    fp=open(ORTH_to_TLRC,'w')
    cmd="cat_matvec %s -I %s" % (tagalign_xform,talairach_xform_inv)
    subprocess.call(cmd,shell=True,stdout=fp);

message = """
  If you want to use multisphere, you'll have to run:
  localspheres -s multisphere.shape -d DSdirectory
  to create the model.

  To look at your surface, either use:
    surf=ortho_brainhull
    suma -i_ply $surf.ply2fid
  or
    afni -niml -dset ortho+orig mask+orig &
    suma -niml -i_ply $surf.ply -sv ortho+orig -novolreg
    <press 't' in suma window>"""

print(message)
