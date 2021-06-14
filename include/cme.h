// cme.h -- indices & definitions for current multipole expansion
//	for dipole + quadrupole terms. For realistic (non-spherical)
//	conductor, there are 3 dipole terms + 8 unique quadrupole terms.
//	For the single homogeneous sphere model there are 2 dipole terms
//	+ 5 unique quadrupole terms. These are shown here (terms in
//	brackets are either zero or are dependent):
//
//			single homogeneous spherical conductor
//				Qx'					Qy'						<Qr = 0>
//				dQx'/dx'			dQx'/dy'				dQx'/dr
//				dQy'/dx'			<dQy'/dy'>				dQy'/dr
//				<dQr/dx' ~ Qx'>		<dQr/dy' ~ Qy'>			<dQr/dr = 0>
//
//			realistic (non-spherical) conducting model
//				Qx					Qy						Qz
//				dQx/dx				dQx/dy					dQx/dz
//				dQy/dx				dQy/dy					dQy/dz
//				dQz/dx				dQz/dy					<dQz/dz>
//
//	Author: Stephen E. Robinson (with much gratitude to Guido Nolte)
//

#ifndef H_CME
#define H_CME

// current multipoles for non-spherical conductor (3 dipole + 9 quadrupole terms)
#define NM			12			// number of multipole terms
#define N3D			11			// number of independent multipole terms in realistic head model
#define QX_			0			// Qx dipole moment
#define	QY_			1			// Qy dipole moment
#define QZ_			2			// Qz dipole moment
#define dQXX_		3			// dQx/dx quadrupole moment
#define dQXY_		4			// dQx/dy quadrupole moment
#define dQXZ_		5			// dQx/dz quadrupole moment
#define dQYX_		6			// dQy/dx quadrupole moment
#define dQYY_		7			// dQy/dy quadrupole moment
#define dQYZ_		8			// dQy/dz quadrupole moment
#define dQZX_		9			// dQz/dx quadrupole moment
#define dQZY_		10			// dQz/dy quadrupole moment
#define dQZZ_		11			// dQz/dz quadrupole moment

// current multipoles for spherical conductor (2 dipole + 5 quadrupole terms)
#define N2D			7			// number of multipole terms
#define QU_			0			// Qx' dipole moment
#define	QV_			1			// Qy' dipole moment
#define dQUU_		2			// dQx'/dx' quadrupole moment
#define dQUV_		3			// dQx'/dy' quadrupole
#define dQUR_		4			// dQx'/dr quadrupole
#define dQVU_		5			// dQy'/dx' quadrupole
#define dQVR_		6			// dQy'/dr quadrupole

// defines for 7 quadrupole terms in spherical conductor
#define NQ7			7
#define dQuu_		0			// dQx'/dx' quadrupole
#define dQuv_		1			// dQx'/dy' quadrupole
#define dQur_		2			// dQx'/dr quadrupole
#define dQvu_		3			// dQy'/dx' quadrupole
#define dQvr_		4			// dQy'/dr quadrupole
#define dQru_		5			// dQr/dx' quadrupole
#define dQrv_		6			// dQr/dy' quadrupole

// defines for 3 dipole terms for realistic non-spherical conductor
#define ND3			3			// number of independent dipole terms
#define Qx_			0			// Qx dipole moment
#define Qy_			1			// Qy dipole moment
#define Qz_			2			// Qz dipole moment

// defines for 8 quadrupole terms for realistic non-spherical conductor
#define NQ3			8			// number of independent quadrupole terms
#define dQxx_		0			// dQx/dx quadrupole moment
#define dQxy_		1			// dQx/dy quadrupole moment
#define dQxz_		2			// dQx/dz quadrupole moment
#define dQyx_		3			// dQy/dx quadrupole moment
#define dQyy_		4			// dQy/dy quadrupole moment
#define dQyz_		5			// dQy/dz quadrupole moment
#define dQzx_		6			// dQz/dx quadrupole moment
#define dQzy_		7			// dQz/dy quadrupole moment
#define dQzz_		8			// dQz/dz quadrupole moment

// defines for 5-component multipoles in non-spherical conductor
#define NQ5			5			// number of independent terms
#define CQXX_		0			// dQx/dx quadrupole moment
#define CQYY_		1			// dQy/dy quadrupole moment
#define	CMx_		2			// Mx (identical to dQy/dz - dQz/dy)
#define CMy_		3			// My (identical to dQx/dz - dQz/dz)
#define CMz_		4			// Mz (identical to dQx/dy - dQy/dx)

// defines for 11-component multipoles in non-spherical conductor, using symmetric & antisymmetric terms
#define NM11		11			// number of independent terms
#define	AX_			0			// Qx dipole moment
#define	AY_			1			// Qy dipole moment
#define AZ_			2			// Qz dipole moment
#define AsXX_		3			// dQx/dx quadrupole moment
#define AsYY_		4			// dQy/dy quadrupole moment
#define AsYZ_		5			// (dQy/dz + dQz/dy) / 2 -- antisymmetric term in yz plane
#define AsXZ_		6			// (dQx/dz + dQz/dx) / 2 -- antisymmetric term in xz plane
#define AsXY_		7			// (dQx/dy + dQy/dx) / 2 -- antisymmetric term in xy plane
#define	AaYZ_		8			// (dQy/dz - dQz/dy) / 2 -- symmetric term in yz plane, identical to Mx
#define AaXZ_		9			// (dQx/dz - dQz/dx) / 2 -- symmetric term in xz plane, identical to My
#define AaXY_		10			// (dQx/dy - dQy/dx) / 2 -- symmetric term in xy plane, identical to Mz

// defines for 8-component multipoles in non-spherical conductor
#define NM8			8
#define _Qx_		0			// Qx dipole moment
#define _Qy_		1			// Qy dipole moment
#define _Qz_		2			// Qz dipole moment
#define _Qxx_		3			// dQx/dx quadrupole moment
#define _Qyy_		4			// dQy/dy quadrupole moment
#define _Mx_		5			// Mx (identical to dQy/dz - dQz/dy)
#define _My_		6			// My (identical to dQx/dz - dQz/dz)
#define _Mz_		7			// Mz (identical to dQx/dy - dQy/dx)


#endif	// H_CME
