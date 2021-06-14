// LtoB() -- compute field at observation point due to
//  current through filamentary wire loop
//
//  Author: Stephen E. Robinson
//          Neuromagnetism Laboratory
//          Dept.of Neurology
//          Henry Ford Hospital
//          Detroit, MI
//          Copyright (c) 2007
//

#include <math.h>
#include <geoms.h>

void    LtoB(
    FIELD       *Loop,          // circular loop origin (vector length is proportional to current)
    double      Pnt[3],         // field observation point
    double      a,              // loop radius (meters)
    double      current,        // current
    double      B[3]                // magnetic field at observation point (Tesla)
)
{
    double      A[3];           // position difference vector
    double      R[3];           // radial unit vector
    double      Z[3];           // axial unit vector
    double      i;              // loop current
    double      xdy;            // x<dot>y
    double      x2;
    double      r;              // radial displacement
    double      z;              // axial displacement
    double      a2;             // a^2
    double      r2;             // r^2
    double      z2;             // z^2
    double      Km;             // K(m)
    double      Em;             // E(m)
    double      k2;             // k^2 -- argument to elliptic integrals
    double      Br;             // radial field
    double      Bz;             // axial field
    double      tmp;
    double      tmp2;
    double      q0, q1, q2, q3, q4;
    double      a0, a1, b0, b1, c0, c1, e0;
    static double   mu0_2pi = 2.0e-07;  // mu0/(2*pi)
    register int    v;          // vector-index
#define MU0         (4.0e-07 * M_PI)


    // extract loop current from normal vector, & make unit axial vector
    for(v=X_, tmp2=0.; v<=Z_; v++)
        tmp2 += Loop->v[v] * Loop->v[v];    // i^2 -- loop normal vector^2 (current^2)
    tmp = sqrt(tmp2);               // current (Ampere-turns) -- magnitude of loop axial vector
    for(v=X_; v<=Z_; v++)
        Z[v] = Loop->v[v] / tmp;        // derive loop axial unit vector

    i = current;

    // determine axial & radial distance from loop origin to observation point
    for(v=X_, x2=xdy=0.; v<=Z_; v++) {
        A[v] = Pnt[v] - Loop->p[v]; // vector pointing from loop to observation positions
        xdy += A[v] * Z[v];     // dot product of loop normal w/ observation point vector
        x2 += A[v] * A[v];      // x^2 -- hypotenuse^2
    }
    z = xdy;                // axial distance
    r = sqrt(x2 - z * z);           // radial distance (at least, according to Pythagorus)
    if(xdy > 0.)                // radial normal vector from observation point to loop
        for(v=X_; v<=Z_; v++)
            R[v] = (Pnt[v] - (Loop->p[v] - z * Z[v])) / r;
    else
        for(v=X_; v<=Z_; v++)
            R[v] = (Pnt[v] - (Loop->p[v] + z * Z[v])) / r ;

    // compute useful squared values
    a2 = a * a;
    r2 = r * r;
    z2 = z * z;

    // compute k^2 & complete elliptic integrals E(m) & K(m)
    //
    //          4 * a * r
    //  k^2 = -------------------
    //      (a + r)^2 + z^2
    //
    q0 = a + r;
    q1 = q0 * q0 + z2;
    k2 = 4 * a * r / q1;        // this is the argument k^2
    a0 = 1.;
    b0 = sqrt(1. - k2);
    c0 = k2;
    e0 = .5;
    Em = 1. - e0 * c0;
    do {                // this will probably converge in ~4 passes
        a1 = 0.5 * (a0 + b0);
        b1 = sqrt(a0 * b0);
        c1 = a1 * a1 - b1 * b1;
        e0 *= 2.;
        Em -= e0 * c1;
        a0 = a1;
        b0 = b1;
        c0 = c1;
    } while(c0 > 1.0e-15);      // or whatever resolution we require
    Km = M_PI / (2. * a0);
    Em *= Km;

    // compute radial & axial field components
    //
    //          mu0 * i * z              (a^2 + r^2 + z^2) * Em
    //  Br = ------------------------------------- * (-Km + ------------------------)
    //        2 * pi * r * ((a + r)^2 + z^2)^-1/2            (a - r)^2 + z^2
    //
    //          mu0 * i             (a^2 - r^2 - z^2) * Em
    //  Bz = --------------------------------- * (Km + ------------------------)
    //        2 * pi * ((a + r)^2 + z^2)^-1/2          (a - r)^2 + z^2
    //
    q2 = sqrt(q1);
    q0 = a - r;
    q3 = q0 * q0 + z2;
    q3 = Em / q3;
    q4 = mu0_2pi * i;
    Br = q4 * z / (r * q2) * ((a2 + r2 + z2) * q3 - Km);
    Bz = q4 / q2 * ((a2 - r2 - z2) * q3 + Km);

    // project cylindrical field components into cartesian space
    for(v=X_; v<=Z_; v++)
        B[v] = Bz * Z[v] + Br * R[v];
}
