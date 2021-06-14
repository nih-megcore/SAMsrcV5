// cme.h -- indices & definitions for current multipole expansion
//	for dipole + quadrupole terms. For realistic (non-spherical)
//	conductor, there are 3 dipole terms + 8 unique quadrupole terms.
//
//	realistic (non-spherical) conducting model
//	Qx		Qy		Qz
//	dQx/dx		dQx/dy		dQx/dz
//	dQy/dx		dQy/dy		dQy/dz
//	dQz/dx		dQz/dy		<dQz/dz>
//
//	Author: Stephen E. Robinson (with much gratitude to Guido Nolte)
//

#define DELTA		1.0e-05		// 10 um derivative on either side of expansion origin

// current multipoles for non-spherical conductor (3 dipole + 8 quadrupole terms)
#define N3D		11			// number of independent multipole terms in realistic head model
#define QX_		0			// Qx dipole moment
#define	QY_		1			// Qy dipole moment
#define QZ_		2			// Qz dipole moment
#define dQXX_		3			// dQx/dx quadrupole moment
#define dQXY_		4			// dQx/dy quadrupole moment
#define dQXZ_		5			// dQx/dz quadrupole moment
#define dQYX_		6			// dQy/dx quadrupole moment
#define dQYY_		7			// dQy/dy quadrupole moment
#define dQYZ_		8			// dQy/dz quadrupole moment
#define dQZX_		9			// dQz/dx quadrupole moment
#define dQZY_		10			// dQz/dy quadrupole moment
