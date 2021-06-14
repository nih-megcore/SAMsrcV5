#! /usr/bin/env python

import sys, os
from pylab import *
from tics import scale1

def plotst(s, srate, **kwargs):
    """
    Display a time-frequency plot as computed by st().
    Various normalizations are possible. None, versus a control plot,
    or versus the time average. The plot is logarithmic unless you say
    -e, in which case it straight power (i.e., based on zero); use this if you
    want to emphasize power increases.

    Other options:

	baseline=[t0 t1] Calculate a time average for the given time interval.

	title=str       The title for the plot. Defaults to "".

	range=[lo, hi]  Upper and lower bounds for the color bar.

	cmap=colormap   Default: cm.jet

	save=file       Save the plot in the named file. Extension defaults to .png.
    """

    titlestr = None
    norm = 'default'
    logscale = True
    rlo = None
    plotfile = None
    start = None
    end = None
    lo = None
    cmap = cm.jet

    for opt, arg in kwargs.items():
	if opt == 'start':
	    start = float(arg) / srate
	elif opt == 'end':
	    end = float(arg) / srate
	elif opt == 'lo':
	    lo = float(arg)
	elif opt == 'hi':
	    hi = float(arg)
	elif opt == 'title':
	    titlestr = arg
	elif opt == 'norm':
	    norm = arg
	elif opt == 'logscale':
	    logscale = arg
	elif opt == 'range':
	    if len(arg) != 2:
		raise RuntimeError, "plotst:range needs a pair of numbers"
	    rlo = float(arg[0])
	    rhi = float(arg[1])
	elif opt == 'save':
	    plotfile = arg
	elif opt == 'cmap':
	    cmap = arg

    if start is None:
	start = 0
    if end is None:
	end = float(s.shape[1]) / srate

    #[marker, start, end, srate, lo, hi, xdim, ydim, isstat] = l

    dual = False
    """
    brikname = args[0]
    l = read_brik(brikname)
    if len(args) == 2:
	    dual = True
	    contbrikname = args[1]
	    l2 = read_brik(contbrikname)
	    if l[3:] != l2[3:]:
		    printerror("brik parameters do not match")

    if titlestr is None:
	    titlestr = marker
	    if titlestr == "None":
		    titlestr = brikname
	    if dual:
		    titlestr = "%s / %s" % (marker, l2[0])
    """

    def avg_over_time(s):
	"""Average over time and return an array of the same size
	where each column (time point) is the same."""

	x = add.reduce(s, 1) / s.shape[1]
	x.shape = (x.shape[0], 1)
	x = repeat(x, s.shape[1], 1)
	return x

    # Normalization

    if 'default'.startswith(norm):
	pass
	#if not logscale:
	#    s = exp(s)
    elif 'norm'.startswith(norm):
	s /= avg_over_time(s)

    if logscale:
	s = log(s)

    figure()

    n = min(minimum.reduce(s))
    m = max(maximum.reduce(s))

    # Sanity check.
    if n == m:
	    print "no data"
	    sys.exit()

    if rlo is not None:
	    n = rlo
	    m = rhi

    nlevels = 40
    clevel = linspace(n, m, nlevels)
    ticks, mticks = scale1(clevel[0], clevel[-1])
    hold(True)
    time = linspace(start, end, s.shape[1])
    if lo is None:
	lo = 0
	hi = s.shape[0]
    fr = linspace(lo, hi, s.shape[0])
    contourf(time, fr, s, clevel, cmap = cmap)
    ax = gca()
    ax.set_xlim(start, end)
    ax.set_ylim(lo, hi)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Frequency (Hz)')
    colorbar(format = '%.2g', ticks = ticks)
#    nlevels = 7
#    clevel = linspace(n, m, nlevels)
#    contour(time, fr, s, clevel, colors = 'black')
    title(titlestr, x = 0, horizontalalignment = 'left', fontsize = 15)

    #rlab = "range [%.3g..%.3g]" % (n, m)
    #text(time[-1], fr[-1] + (fr[1] - fr[0]), rlab, fontsize = 10,
    #     horizontalalignment = 'right', verticalalignment = 'bottom')

    if plotfile:
	    (root, ext) = os.path.splitext(plotfile)
	    if len(ext) > 0:
		    savefig(plotfile)
	    else:
		    savefig("%s%s%s" % (plotfile, os.path.extsep, "png"))

