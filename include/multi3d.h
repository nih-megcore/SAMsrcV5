// multi3d.h -- indices & definitions for current multipole expansion
//	for dipole + quadrupole terms. For the single homogeneous
//	sphere model there are 3 dipole terms + 2 unique quadrupole
//	terms derived from the orthogonal dipole vector. These are shown
//	here (terms in brackets are eliminated):
//
//	Author: Stephen E. Robinson (with much gratitude to Guido Nolte)
//

#define DELTA		1.0e-05		// 10 um derivative on either side of expansion origin
#define M3D			5			// 5 moment terms

// current multipole terms for realistic conductor (3 dipoles + 2 quadrupole terms)
#define DX_		0					// x-dipole forward solution
#define DY_		1					// y-dipole forward solution
#define DZ_		2					// z-dipole forward solution
#define Q1_		3					// quadrupole 1
#define Q2_		4					// quadrupole 2
