// LtoB() -- compute field at observation point due to
//	current through a filamentary wire loop
//
//	Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_ellint.h>
#include <geoms.h>

#define MU0         (4.0e-07 * M_PI)
#define MU0_2PI     2.0e-07 // mu0/(2*pi)

void    LtoB(
    FIELD       *Loop,      // circular loop origin & orientation
    double      Pnt[3],     // field observation point
    double      a,          // loop radius (meters)
    double      current,    // current (Amperes)
    double      B[3]        // vector field at observation point
)

{
    double  A[3];           // position difference vector between loop & sensor
    double  R[3];           // radial unit vector
    double  Z[3];           // axial unit vector (orientation of Loop)
    double  Bo;             // field at center of loop at r=0, z=0
    double  Br;             // radial field magnitude
    double  Bz;             // axial field magnitude
    double  alpha;          // r/a
    double  beta;           // z/a
    double  gamma;          // z/r
    double  Q;              // (1+alpha)^2+beta^2
    double  r;              // radial displacement
    double  z;              // axial displacement
    double  a2;             // a^2 -- coil radius^2
    double  r2;             // r^2 -- radial displacement^2
    double  x2;
    double  z2;             // z^2 -- axial displacement^2
    double  Km;             // K(m) -- complete elliptic integral of the 1st kind
    double  Em;             // E(m) -- complete elliptic integral of the 2nd kind
    double  k2;             // k^2 -- argument to complete elliptic integrals
    int     v;              // vector-index

    // Z[] is the axial unit vector (loop orientation)
    for (v=X_, r2=0.; v<=Z_; v++)
        r2 += Loop->v[v] * Loop->v[v];
    r = sqrt(r2);
    for (v=X_; v<=Z_; v++)
        Z[v] = Loop->v[v] / r;

    // compute axial distance to observation point
    for (v=X_, x2=z=0.; v<=Z_; v++) {
        A[v] = Pnt[v] - Loop->p[v];         // A[] is vector pointing from loop to observation positions
        z += A[v] * Z[v];                   // axial displacement (note sign of dot product)
        x2 += A[v]* A[v];                   // x2 is the square of the hypotenuse (distance between loop and observation point)
    }
    z2 = z * z;
    a2 = a * a;

    // determine the radial distance & its unit vector R[]
    r2 = x2 - z2;
    r = sqrt(r2);                           // radial distance (at least, according to Pythagorus)
	if(z > 0.)                              // radial normal vector from observation point to loop
		for(v=X_; v<=Z_; v++)
			R[v] = (Pnt[v] - (Loop->p[v] - z * Z[v])) / r;
	else
		for(v=X_; v<=Z_; v++)
			R[v] = (Pnt[v] - (Loop->p[v] + z * Z[v])) / r ;

    // compute axial & radial field
    //if (z == 0.) {                          // radial field vanishes
    if (fabs(z) < 1.0e-06) {
        //printf("Avoiding Singularity\n");
        Br = 0.;
        Bz = MU0 * current * a2 / (2. * pow((a2 + z2), 1.5));
    } else {                                // finite axial & radial fields
        Bo = current * MU0 / (2. * a);
        alpha = r / a;
        beta = z / a;
        gamma = z / r;
        Q = (1. + a) * (1. + a) + beta * beta;
        k2 = 4. * a / Q;
        Km = gsl_sf_ellint_Kcomp(k2, 1.0e-12);
        Em = gsl_sf_ellint_Ecomp(k2, 1.0e-12);
        Bz = Bo / (M_PI * sqrt(Q)) * (Em * (1. - alpha * alpha - beta * beta) / (Q - 4. * alpha) + Km);
        Br = Bo * gamma / (M_PI * sqrt(Q)) * (Em * (1. + alpha * alpha + beta * beta) / (Q - 4. * alpha) - Km);
    }

    // project radial & axial field into cartesian space
    for (v=X_; v<=Z_; v++)
        B[v] = Br * R[v] + Bz * Z[v];
}
