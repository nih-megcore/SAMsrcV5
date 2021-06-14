// CtoS() -- convert cartesian to spherical coordinates
//
// enter with x, y, and z
//
// exit with radius 'rho' in units of x, y, and z
// and spherical angles theta and phi
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <math.h>
#include <geoms.h>

void	CtoS(
	double	*c,		// Cartesian vector
	double	*s		// spherical vector
) {
    s[RHO_] = sqrt(c[X_] * c[X_] + c[Y_] * c[Y_] + c[Z_] * c[Z_]);
    s[THETA_] = atan2(c[Y_], c[X_]);
    s[PHI_] = atan2(sqrt(c[X_] * c[X_] + c[Y_] * c[Y_]), c[Z_]);
}
