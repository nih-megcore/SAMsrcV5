// LtoB() -- compute field at observation point due to
//	current through filamentary wire loop. The complete
//  elliptic integrals of the 1st & 2nd kind are computed
//  using Gauss's method (arithmetic mean, geometric mean).
//
//	Author:	Stephen E. Robinson
//			Neuromagnetism Laboratory
//			Dept.of Neurology
//			Henry Ford Hospital
//			Detroit, MI
//			Copyright (c) 2007
//
//  Revised: 24 February, 2021
//  Revised: 2 July, 2021 -- scale current by loop magnitude (needed for FitDip to avoid zero partials)
//  Revised: 5 August, 2021 -- use GSL to compute complete elliptic integrals
//  Revised: 18 August 2021 -- GSL requires k -- not k^2 (thanks, Amaia!)
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_ellint.h>
#include <geoms.h>

#define mu0_2pi 2.0e-07     // mu0/(2*pi)


void    LtoB(
    FIELD   *Loop,      // circular loop origin & orientation
    double  Pnt[3],     // field observation point
    double  a,          // loop radius (meters)
    double  current,    // current (Amperes)
    double  B[3]        // vector field at observation point
)
{
    double  A[3];       // position difference vector
    double  R[3];       // radial unit vector
    double  Z[3];       // axial unit vector
    double  Zpnt[3];    // position on z-axis displaced by z
    double  adz;        // a<dot>z
    double  x2;         // hypotenuse^2
    double  r;          // radial displacement
    double  z;          // axial displacement
    double  a2;         // a^2
    double  r2;         // r^2
    double  z2;         // z^2
    double  m;          // magnitude of orientation vector
    double  m2;         // m^2
    double  Km;         // K(m)
    double  Em;         // E(m)
    double  k;          // argument to GSL elliptic integrals
    double  Br;         // radial field
    double  Bz;         // axial field
    double  q0, q1, q2, q3, q4;
    int     v;          // vector-index

    // make Loop orientation a unit axial vector
    for (v=X_, m2=0.; v<=Z_; v++)
        m2 += Loop->v[v] * Loop->v[v];
    m = sqrt(m2);
    if (m == 0.)
        Cleanup("zero-length orientation vector");
    for (v=X_; v<=Z_; v++)
        Z[v] = Loop->v[v] / m;              // derive loop axial unit vector
    current *= m;                           // scale the current by the coil vector magnitude

    // determine axial & radial distance from loop origin to observation point
    for (v=X_, x2=adz=0.; v<=Z_; v++) {
        A[v] = Pnt[v] - Loop->p[v];         // vector pointing from loop to observation positions
        adz += A[v] * Z[v];                 // dot product of loop normal w/ observation point vector
        x2 += A[v] * A[v];                  // x^2 -- hypotenuse^2
    }
    z = adz;                                // axial displacement from loop
    r = sqrt(x2 - z * z);                   // radial displacement from loop (at least, according to Pythagorus)

    // determine the radial displacement & its unit vector R[]
    for (v=X_; v<=Z_; v++)                  // signed coordinate along z-axis (displacement from Loop)
        if (adz < 0.)
            Zpnt[v] = Loop->p[v] -  z * Z[v];
        else
            Zpnt[v] = Loop->p[v] +  z * Z[v];
    for (v=X_; v<=Z_; v++)
        R[v] = (Pnt[v] - Zpnt[v]) / r;

    // compute useful squared values
    a2 = a * a;
    r2 = r * r;
    z2 = z * z;

    // compute k^2 & complete elliptic integrals E(m) & K(m)
    //
    //		    4 * a * r
    //	k^2 = --------------- (use square root for GSL)
    //		  (a + r)^2 + z^2
    //
    q0 = a + r;
    q1 = q0 * q0 + z2;          // (a + r)^2 + z^2
    k = sqrt(4 * a * r / q1);   // this is the argument k
    Km = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
    Em = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);

    // compute radial & axial field components
    //
    //                  mu0 * I * z                    (a^2 + r^2 + z^2) * Em
    //  Br = ------------------------------------- * (------------------------ - Km)
    //	      2 * pi * r * ((a + r)^2 + z^2)^1/2          (a - r)^2 + z^2
    //
    //                    mu0 * I                  (a^2 - r^2 - z^2) * Em
    //  Bz = --------------------------------- * (------------------------ + Km)
    //	      2 * pi * ((a + r)^2 + z^2)^1/2          (a - r)^2 + z^2
    //
    q2 = sqrt(q1);
    q0 = a - r;
    q3 = q0 * q0 + z2;          // (a - r)^2 + z^2
    q3 = Em / q3;
    q4 = mu0_2pi * current;
    Br = q4 * z / (r * q2) * ((a2 + r2 + z2) * q3 - Km);
    Bz = q4 / q2 * ((a2 - r2 - z2) * q3 + Km);

    // project cylindrical field components into cartesian space
    for (v=X_; v<=Z_; v++)
        B[v] = Bz * Z[v] + Br * R[v];
}
