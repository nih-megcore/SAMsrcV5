// MtoB_plus() - compute field due to magnetic dipole
//  near-field correction is sqrt(area)/r = ~5.406e-03/r
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <geoms.h>
#include <lfulib.h>
#include <samlib.h>
#include <jig.h> // for AM2 definition

double  MtoB_plus(          // returns scalar field
    FIELD       *Dipole,    // magnetic dipole position & vector
    OPMFIELD    *Sensor,    // sensor position & vector
    double      current     // magnitude of coil current (Amperes)
)
{
    FIELD       M;          // scaled magnetic dipole
    double      B[3];       // vector magnetic field at observation point
    double      R[3];       // unit vector r -- between dipole & measurement point
    double      xyz[3];     // cartesian vesrion of orientation
    double      rtp[3];     // spherical version of orientation
    double      Bs;         // scalar field
    double      r;          // |r|
    double      r2;         // r^2
    double      r3;         // r^3
    double      mdr;        // M<dot>R
    double      gain;       // gain of the sensor
    int         v;          // vector-index

    // copy M & scale the dipole moment vector
    for (v=X_; v<=Z_; v++) {
        M.p[v] = Dipole->p[v];
        M.v[v] = AM2 * current * Dipole->v[v];
    }

    rtp[0]=1.;
    rtp[1]=Sensor->v[0];
    rtp[2]=Sensor->v[1];
    StoC(rtp, xyz);
    gain = Sensor->g;

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
    return(Bs*gain);
}
