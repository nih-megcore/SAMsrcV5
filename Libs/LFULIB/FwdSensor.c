// FwdSensor() -- compute field sensed by arbitrary gradiometer due to current dipole
//  This version features programmable integration order from 1st to 7th
//
// Author:  Stephen E. Robinson (after Sarvas, Roth, & Sato)
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <math.h>
#include <geoms.h>
#include <samlib.h>

double FwdSensor(
    FIELD       *Q,         // location & orientation of current dipole
    SENSOR      *Sensor,    // sensor channel (magnetometer, gradiometer) structure
    int         Order       // integration order: 1, 2, 3, 4, 5, & 7 are implemented
) {
    double      R[3];       // vector from Sp to B
    double      R0[3];      // vector from Sp to Q
    double      A[3];       // vector from B to Q
    double      V[3];
    double      B0[3];
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
    double      bc;         // magnetic flux sensed by coil
    double      bs;         // total field sensed by all coils
    int         c;          // coil index
    int         i;          // integration index
    int         v;          // vector index

    // sum field detected by all coils
    for (c=0, bs=0.; c<Sensor->NumCoils; c++) {

        // integrate coil computed field at each integration point
        for (i=0, bc=0.; i<Sensor->Coil[c].NumInt; i++) {

            // compute vectors A, R0, & R, & scalars R0<dot>R, A<dot>R, |R|^2, & |A|^2
            for (v=X_, r2=a2=r0dr=adr=0.; v<=Z_; v++) {
                R[v] = Sensor->Coil[c].B[i].p[v] - Sensor->LocalSphere[v];
                R0[v] = Q->p[v] - Sensor->LocalSphere[v];
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

            // compute B
            for (v=X_; v<=Z_; v++)
                B0[v] = k3 * (f * V[v] - vdr * DELF[v]);
            for (v=X_; v<=Z_; v++)
                bc += Sensor->Coil[c].w[i] * Sensor->Coil[c].origin.v[v] * B0[v];
        }
        bs += bc * (double)Sensor->Coil[c].SenseTurns;
    }
    return bs / (double)abs(Sensor->Coil[0].SenseTurns);
}
