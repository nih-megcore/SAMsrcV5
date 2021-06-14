from math import pow, floor, ceil, log10, fmod
from numpy import arange

def scale1(lo, hi, ticint = None):

    """tl, mticint = scale1(lo, hi[, ticint = <nint>])
    Given data in an interval from lo to hi, calculate nice tic marks.
    ticint, if specified, gives the initial number of intervals to try.
    The tics are returned as a list in tl, and a suggested number of minor
    tic intervals is returned in mticint."""

    if ticint is None:
	ticint = 5

    lo = float(lo)
    hi = float(hi)

    # Reference: Lewart, C. R., "Algorithms SCALE1, SCALE2, and
    # SCALE3 for Determination of Scales on Computer Generated
    # Plots", Communications of the ACM, 10 (1973), 639-640.
    # Also cited as ACM Algorithm 463.

    a, b = magform((hi - lo) / ticint)
    if a < 1.41:    # sqrt(2)
	x = 1
    elif a < 3.16:  # sqrt(10)
	x = 2
    elif a < 7.07:  # sqrt(50)
	x = 5
    else:
	x = 10
    if b < 0:
	sep = x * pow(10., b)
    else:
	sep = float(x * (10 ** b))
    if x == 10:
	x = 1

    # The following guarantees that if zero is in the range, it will be
    # included as a tic.

    a = floor(lo / sep)
    b = ceil(hi / sep)
    return arange(a, b + 1.) * sep, x

def magform(x):

    """magform(x)
    Returns (a, b), where x = a * 10^b, a >= 1., and b is integral."""

    if x == 0:
	return 0., 0
    l = log10(abs(x))
    r = fmod(l, 1.)
    a, b = pow(10., r), int(l - r)
    if a < 1.:
	a, b = a * 10., b - 1
    if x < 0.:
	a = -a
    return a, b
