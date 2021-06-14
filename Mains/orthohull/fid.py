# Compute the FID basis.

import numpy, numpy.linalg

# FID names:
NASION = 'Nasion'
LEAR = 'Left Ear'
REAR = 'Right Ear'

def cross(u, v):
  # x     y      z      z      y
    x = u[1] * v[2] - u[2] * v[1]
    y = u[2] * v[0] - u[0] * v[2]
    z = u[0] * v[1] - u[1] * v[0]
    return numpy.array([x, y, z])

def normalize(x):
    return x / numpy.sqrt(numpy.dot(x, x))

def fid(nasn, lear, rear):
    """m = fid(nasn, lear, rear)
Construct the rotation into the FID coordinate system. The FID coordinate
basis is constructed from the fiducial points as follows: The origin is the
center of the line from the left to the right preauricular points. The +x
axis passes through the nasion. The +z axis is normal to the plane
containing the fiducial points, and the +y axis is orthogonal to both x and
z. The transform is returned as a 4x4 matrix which includes the translation
to the origin, used with fid_transform(). Scaling must be done
separately."""

    nasn = numpy.asarray(nasn)
    lear = numpy.asarray(lear)
    rear = numpy.asarray(rear)

    # Calculate the origin.

    o = rear + (lear - rear) / 2.

    # The X axis.

    x = normalize(nasn - o)

    # The Z axis.

    z = normalize(cross(rear - o, x))

    # The Y axis.

    y = normalize(cross(z, x))

    # The transformation matrix from the standard basis to the FID basis.

    m = numpy.identity(4, 'd')
    m[0:3, 0] = x
    m[0:3, 1] = y
    m[0:3, 2] = z
    m[0:3, 3] = o

    # Invert it to get the transformation into the FID basis.

    return numpy.linalg.inv(m)

def fid_transform(m, v):
    """fid_transform(m, v)
Returns a vector rotated and translated by the 4x4 transform m."""

    r = m[0:3, 0:3]
    t = m[0:3, 3]
    return numpy.inner(r, v) + t

def fid_rotate(m, v):
    """fid_rotate(m, v)
Returns a vector only rotated by m, useful for normals."""

    r = m[0:3, 0:3]
    return numpy.inner(r, v)

