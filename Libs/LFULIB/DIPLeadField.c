// DIPLeadField() -- compute L x 3 lead field matrix L for current dipole in
//  sphere or multisphere model, for an array of sensors. Note that each
//  sensor has a sphere origin that is set to either an identical (single)
//  or independent (multiple) sphere origins.
//
//  Author: Stephen E. Robinson (after Nolte, Sarvas, & Hamalainen)
//          MEG Core Facility
//          NIMH
//

#include <stdlib.h>
#include <lfulib.h>
#include <math.h>
#include <model.h>

void DIPLeadField(
    double              Vc[3],          // dipole origin
    gsl_matrix          *L,             // L - Mx3 lead field matrix
    HeaderInfo          *Header,        // header information, including primary & reference sensor indices
    ChannelInfo         *Channel,       // Channel[M] -- channel structure with primary and reference sensors
    int                 Order           // integration order: 1, 2, 3, 4, 5, & 7 are implemented
) {
    FIELD               Q6;             // test dipole
    double              *B;             // B[M] -- forward solution for one lead-field component
    int                 m;              // sensor index
    int                 M;              // number of sensors

    // allocate memory
    M = Header->NumPri;
    B = new_arrayE(double, M, "B[]");

    // set dipole origin for all dipole terms
    Q6.p[X_] = Vc[X_];  Q6.p[Y_] = Vc[Y_];  Q6.p[Z_] = Vc[Z_];

    // (1) compute Qx
    Q6.v[X_] = 1.;  Q6.v[Y_] = 0.;  Q6.v[Z_] = 0.;
    FwdVec(&Q6, B, Header, Channel, Order);
    for (m=0; m<M; m++)
        gsl_matrix_set(L, m, X_, B[m]);

    // (2) compute Qy
    Q6.v[X_] = 0.;  Q6.v[Y_] = 1.;  Q6.v[Z_] = 0.;
    FwdVec(&Q6, B, Header, Channel, Order);
    for (m=0; m<M; m++)
        gsl_matrix_set(L, m, Y_, B[m]);

    // (3) compute Qz
    Q6.v[X_] = 0.;  Q6.v[Y_] = 0.;  Q6.v[Z_] = 1.;
    FwdVec(&Q6, B, Header, Channel, Order);
    for (m=0; m<M; m++)
        gsl_matrix_set(L, m, Z_, B[m]);

    free(B);
}
