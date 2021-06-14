// FwdOPM() -- compute field sensed by OPM magnetometer due to current dipole
//  in homogeneous spherical volume conductor
//
// Author:  Stephen E. Robinson (after Sarvas)
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <math.h>
#include <geoms.h>
#include <samlib.h>

double FwdOPM (
    FIELD   *Q,             // location & orientation of current dipole
    OPM     *OPMSensor      // OPM magnetometer structure
) {
    double      R[3];       // vector from Sp to B
    double      R0[3];      // vector from Sp to Q
    double      A[3];       // vector from B to Q
    double      V[3];
    double      DELF[3];
    double      r2;
    double      a2;
    double      r;
    double      a;
    double      r0dr;
    double      adr;
    double      vdr;
    double      f;
    double      k0;
    double      k1;
    double      k2;
    double      k3;
    double      bs;         // total field sensed by all coils
    int         v;          // vector index

    // compute vectors A, R0, & R, & scalars R0<dot>R, A<dot>R, |R|^2, & |A|^2
    for (v=X_, r2=a2=r0dr=adr=0.; v<=Z_; v++) {
        R[v] = OPMSensor->origin.p[v] - OPMSensor->LocalSphere[v];
        R0[v] = Q->p[v] - OPMSensor->LocalSphere[v];
        A[v] = R[v] - R0[v];
        r0dr += R0[v] * R[v];
        adr += A[v] * R[v];
        r2 += R[v] * R[v];
        a2 += A[v] * A[v];
    }

    // compute V = Q<cross>R0
    V[X_] = Q->v[Y_] * R0[Z_] - Q->v[Z_] * R0[Y_];
    V[Y_] = Q->v[Z_] * R0[X_] - Q->v[X_] * R0[Z_];
    V[Z_] = Q->v[X_] * R0[Y_] - Q->v[Y_] * R0[X_];

    // compute |R|, |A|, f, & coefficients for <del>f
    r = sqrt(r2);
    a = sqrt(a2);
    f = a * (r * a + r2 - r0dr);
    k0 = adr / a + 2. * r;
    k1 = a2 / r + 2. * a + k0;
    k2 = a + k0;
    k3 = MU0_4PI / (f * f);

    // compute V<dot>R & <del>f
    for (v=X_, vdr=0.; v<=Z_; v++) {
        vdr += V[v] * R[v];
        DELF[v] = k1 * R[v] - k2 * R0[v];
    }

    // compute scalar field response
    for (v=X_, bs=0.; v<=Z_; v++)
        bs += k3 * (f * V[v] - vdr * DELF[v]) * OPMSensor->origin.v[v];

    return bs;
}
