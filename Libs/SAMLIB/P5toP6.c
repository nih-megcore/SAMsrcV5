// P5toP6() -- convert 5-parameter dipole representation
//	(rho, theta, phi, psi, q) to 6-parameter
//	(x, y, z, i, j, k).  Note that 6-parameter is ALL
//	Cartesian, and 5-parameter is ALL spherical.
//	There are no mixed notations.  Angles are in
//	radians.
//
//	Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdlib.h>
#include <math.h>
#include <samlib.h>


void	P5toP6(
    FIELD   *p5,
    FIELD   *p6
) {
	double	cart[3];	// Cartesian vector
	double	sphr[3];	// spherical vector
	double	sinp;		// sin(phi)
	double	cosp;		// cos(phi)
	double	sint;		// sin(theta)
	double	cost;		// cos(theta)
	double	x;			// temporary x
	double	y;			// temporary y
	double	z;			// temporary z

	// convert position vector to Cartesian & move to output
	sphr[RHO_] = p5->p[RHO_];
	sphr[THETA_] = p5->p[THETA_];
	sphr[PHI_] = p5->p[PHI_];
	StoC(sphr, cart);
	p6->p[X_] = cart[X_];
	p6->p[Y_] = cart[Y_];
	p6->p[Z_] = cart[Z_];

	// compute trig parameters from spherical position vector (note sign)
	sinp = sin(-p5->p[PHI_]);
	cosp = cos(-p5->p[PHI_]);
	sint = sin(-p5->p[THETA_]);
	cost = cos(-p5->p[THETA_]);

	// convert orientation vector to Cartesian
	sphr[RHO_] = p5->v[Q_];
	sphr[THETA_] = M_PI_2;
	if(p5->v[PSI_] >= 0.)
		sphr[PHI_] = p5->v[PSI_] - M_PI;
	else
		sphr[PHI_] = p5->v[PSI_] + M_PI;
	StoC(sphr, cart);

	// rotate vector around y-axis by minus theta radians
	x = cost * cart[X_] - sint * cart[Z_];
	y = cart[Y_];
	z = cost * cart[Z_] + sint * cart[X_];

	// rotate vector around z-axis by minus phi radians
	p6->v[X_] = cosp * x + sinp * y;
	p6->v[Y_] = cosp * y - sinp * x;
	p6->v[Z_] = z;
}
