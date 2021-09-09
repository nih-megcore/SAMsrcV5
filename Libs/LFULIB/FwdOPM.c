// FwdOPM() -- compute field sensed by OPM magnetometer due to current dipole
//  in homogeneous spherical volume conductor. The sensing region for the
//  FieldLine OPMs is a cylinder ~1.5 mm diameter and 3 mm long (on the
//  sensing axis). This routine does 3-point integration along the sensing
//  axis,
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
    double      B[3];	    // vector field
    double      R[3];       // vector from sphere origin to B
    double      R0[3];      // vector from sphere origin to Q
    double      A[3];       // vector from B to Q
    double      V[3];	    // Q cross R0
    double      Z[3];       // sensor unit orientation vector
    double      Zpnt[3];    // observation points along sensing axis
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
    double      bs;         // scalar field sensed by OPM 
    double      bi;         // integrated sensitivity
    int         v;          // vector index
    int         w;          // axial integration index

    
    // assure unit sensing axis
    for (v=X_, r2=0.; v<=Z_; v++)
        r2 += OPMSensor->origin.v[v] * OPMSensor->origin.v[v];
    r = sqrt(r2);
    for (v=X_; v<=Z_; v++)
        Z[v] = OPMSensor->origin.v[v] / r;
    
    //  perform 3-point integration along sensing axis
    for (w=0, bi=0.; w<3; w++) {

      // displace the measurement point along the sensor orientation axis
      for (v=X_; v<=Z_; v++)
            Zpnt[v] = OPMSensor->origin.p[v] + (double)w * 0.001 * Z[v];

      // compute vectors A, R0, & R, & scalars R0<dot>R, A<dot>R, |R|^2, & |A|^2
      for (v=X_, r2=a2=r0dr=adr=0.; v<=Z_; v++) {
          R[v] = Zpnt[v] - OPMSensor->LocalSphere[v];
          R0[v] = Q->p[v] - OPMSensor->LocalSphere[v];
          A[v] = R[v] - R0[v];
          r0dr += R0[v] * R[v];
          adr += A[v] * R[v];
          r2 += R[v] * R[v];
          a2 += A[v] * A[v];
      }

      // compute |R|, |A|, f, & coefficients for <del>f
      r = sqrt(r2);
      a = sqrt(a2);
      f = a * (r * a + r2 - r0dr);
      k0 = adr / a + 2. * r;
      k1 = a2 / r + 2. * a + k0;
      k2 = a + k0;
      k3 = MU0_4PI / (f * f);

      // compute V = Q<cross>R0
      V[X_] = Q->v[Y_] * R0[Z_] - Q->v[Z_] * R0[Y_];
      V[Y_] = Q->v[Z_] * R0[X_] - Q->v[X_] * R0[Z_];
      V[Z_] = Q->v[X_] * R0[Y_] - Q->v[Y_] * R0[X_];

      // compute V<dot>R & <del>f
      for (v=X_, vdr=0.; v<=Z_; v++) {
          vdr += V[v] * R[v];
          DELF[v] = k1 * R[v] - k2 * R0[v];
      }

      // compute vector field at sensor
      for (v=X_; v<=Z_; v++)
            B[v] = k3 * (f * V[v] - vdr * DELF[v]);

      // compute sensor scalar field response
      for (v=X_, bs=0.; v<=Z_; v++)
          bs += B[v] * Z[v];

      // sum integration values
      bi += bs / 3.;
    
    }

    return bi;
}
