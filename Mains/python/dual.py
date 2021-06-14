# Compute the dual of a set of basis vectors.

from numpy import *
from numpy.linalg import *

def normalize(v):
	return v / sqrt(dot(v, v))

def dual(vv, norm = 1):
	"""dual(v) returns an array whose rows form the dual basis of the
normalized rows of v. dual(v, 0) doesn't normalize first."""
	vv = asarray(vv, 'd')
	n = vv.shape[0]
	nv = zeros(vv.shape, 'd')
	for i in range(n):
		if norm:
			nv[i] = normalize(vv[i])
		else:
			nv[i] = vv[i]
	a = zeros((n, n), 'd')
	for i in range(n):
		for j in range(i+1):
			a[i, j] = a[j, i] = dot(nv[i], nv[j])
	return dot(inv(a), nv)
