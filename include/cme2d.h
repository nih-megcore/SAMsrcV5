// cme.h -- indices & definitions for current multipole expansion
//	for dipole + quadrupole terms. For the single homogeneous
//	sphere model there are 2 dipole terms + 5 unique quadrupole
//	terms. These are shown here (terms in brackets are either zero
//	or are dependent):
//
//		single homogeneous spherical conductor
//	Qx'			Qy'			<Qr = 0>
//	dQx'/dx'		dQx'/dy'		dQx'/dr
//	dQy'/dx'		<dQy'/dy'>		dQy'/dr
//	<dQr/dx' ~ Qx'>		<dQr/dy' ~ Qy'>		<dQr/dr = 0>
//
//	Author: Stephen E. Robinson (with much gratitude to Guido Nolte)
//		MEG Core Facility
//		NIMH
//

#define DELTA		1.0e-05		// 10 um derivative on either side of expansion origin

// current multipoles for spherical conductor (2 dipole + 5 quadrupole terms)
#define N2D		7		// number of multipole terms
#define QU_		0		// Qx' dipole moment
#define	QV_		1		// Qy' dipole moment
#define dQUU_		2		// dQx'/dx' quadrupole moment
#define dQUV_		3		// dQx'/dy' quadrupole
#define dQUR_		4		// dQx'/dr quadrupole
#define dQVU_		5		// dQy'/dx' quadrupole
#define dQVR_		6		// dQy'/dr quadrupole
