// StoC() -- convert spherical coordinates to cartesian
//
// entr with rho & theta, phi (radians)
//
// enter with x, y, and z
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <math.h>
#include <geoms.h>

void    StoC(
    double  *s, // spherical vector
    double	*c  // Cartesian vector
) {
    c[X_] = s[RHO_] * sin(s[PHI_]) * cos(s[THETA_]);
    c[Y_] = s[RHO_] * sin(s[PHI_]) * sin(s[THETA_]);
    c[Z_] = s[RHO_] * cos(s[PHI_]);
}
