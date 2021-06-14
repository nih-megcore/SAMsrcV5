// cme2d,4.h -- indices & definitions for current multipole expansion
//	for dipole + quadrupole terms. For the single homogeneous
//	sphere model there are 2 dipole terms + 2 unique quadrupole
//	terms derived from the orthogonal dipole vector. These are shown
//	here (terms in brackets are eliminated):
//
//			single homogeneous spherical conductor
//				Qx'					Qy'						<Qr = 0>
//				dQx'/dx'			dQx'/dy'				dQx'/dr
//				dQy'/dx'			<dQy'/dy'>				dQy'/dr
//				<dQr/dx' ~ Qx'>		<dQr/dy' ~ Qy'>			<dQr/dr = 0>
//
//	Author: Stephen E. Robinson (with much gratitude to Guido Nolte)
//

#define DELTA		1.0e-05		// 10 um derivative on either side of expansion origin

// current multipoles for spherical conductor (2 dipole + 2 quadrupole terms)
#define N4D			4			// number of multipole terms
#define QP_			0			// Qp dipole moment
#define	QO_			1			// Qo dipole moment
#define dQOP_		2			// dQo/dp quadrupole moment
#define dQOR_		3			// dQo/dr quadrupole