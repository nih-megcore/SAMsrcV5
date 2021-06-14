// MtoB() - compute field due to magnetic dipole
//  near-field correction is sqrt(area)/r = ~5.406e-03/r
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <geoms.h>
#include <lfulib.h>
#include <jig.h> // for AM2 definition

double  MtoB(               // returns scalar field
    FIELD       *Dipole,    // magnetic dipole position & vector
    FIELD       *Sensor,    // sensor position & vector
    double      current      // magnitude of coil current (Amperes)
)	

{
    FIELD       M;          // scaled magnetic dipole
    double      B[3];       // vector magnetic field at observation point
    double      R[3];       // unit vector r -- between dipole & measurement point
    double      Bs;         // scalar field
    double      r;          // |r|
    double      r2;         // r^2
    double      r3;         // r^3
    double      mdr;        // M<dot>R
    int         v;          // vector-index

    // copy M & scale moment vector
    for (v=X_; v<=Z_; v++) {
        M.p[v] = Dipole->p[v];
        M.v[v] = AM2 * current * Dipole->v[v];
    }

    // compute vector R & its magnitude
    for (v=X_, r2=0.; v<=Z_; v++) {
        R[v] = Sensor->p[v] - M.p[v];
        r2 += R[v] * R[v];
    }
    r = sqrt(r2);
    r3 = r * r2;
    for (v=X_; v<=Z_; v++)
        R[v] /= r;                          // make R a unit vector

    // compute B
    if (r > 0.) {                           // to evade the singularity at r == 0...
        for (v=X_, mdr=0.; v<=Z_; v++)      // compute M<dot>R
            mdr += M.v[v] * R[v];
        for (v=X_; v<=Z_; v++)
            B[v] = MU0_4PI * ((3. * R[v] * mdr - M.v[v]) / r3);
        for (v=X_, Bs=0.; v<=Z_; v++)
            Bs += Sensor->v[v] * B[v];
    } else {
        Bs = 0.;
    }
    return(Bs);
}
