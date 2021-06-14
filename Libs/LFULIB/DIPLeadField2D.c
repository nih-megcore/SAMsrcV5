// DIPLeadField2D() -- compute lead field matrix L for 2 tangential dipoles
//	for an array of sensors.
//
//	Note: origin/position is given in _headframe_coordinate_ &
//	not local sphere coordinates.
//
//	Author:	Stephen E. Robinson (after Nolte, Sarvas, & Hamalainen)
//          MEG Core Group
//          NIMH
//

#include <stdlib.h>
#include <math.h>
#include <model.h>
#include <lfulib.h>
#include <samlib.h>
#include <samutil.h>


void    DIPLeadField2D(
    double      Vc[3],      // expansion origin (headframe coordinates)
    gsl_matrix  *L,         // Mx2 lead-field matrix (external array allocation)
	HeaderInfo  *Header,    // header information, including primary & reference sensor indices
	ChannelInfo *Channel,   // Channel[M] -- channel structure with primary and reference sensors
	int         Order       // integration order: 1, 2, 3, 4, 5, & 7 are implemented for SQUID sensors, order 0 selects OPMsensor
) {
    FIELD   Q6;             // test dipole
    double  Sp[3];          // local sphere origin
    double  Vx[3];          // vector in the x' tangential direction
    double  Vy[3];          // vector in the y' tangential direction
    double  Vr[3];          // vector in the radial direction
    double  *B;             // forward solution of one lead-field component
    double  r;              // magnitude of r
    double  r2;             // r^2
    int     M;              // sensor count
    int     m;              // sensor index
    int     v;              // vector index

    // set local sphere origin for selecting tangential plane
    m = Header->PsIndex[0];
    for (v=X_; v<=Z_; v++)
        if (Order == 0)
            Sp[v] = Channel[m].Geom.OPMSensor.LocalSphere[v];
        else
            Sp[v] = Channel[m].Geom.MEGSensor.LocalSphere[v];

    // allocate memory
    M = L->size1;
    B = new_arrayE(double, M, "B[]");

	// find normalized position vector relative to sphere origin
	for (v=X_, r2=0.; v<=Z_; v++) {
	    Vc[v] -= Sp[v];
	    r2 += Vc[v] * Vc[v];
	}
	r = sqrt(r2);

    // derive radial vector, & solve tangential plane
    for (v=X_; v<=Z_; v++)
        Vr[v] = Vc[v] / r;
    OrthoPlane(Vr, Vx, Vy);

    // set dipole position to Vc
    for (v=X_; v<=Z_; v++)
        Q6.p[v] = Vc[v];

    // (1) compute field for Qx'
    for (v=X_; v<=Z_; v++)
        Q6.v[v] = Vx[v];
    FwdVec(&Q6, B, Header, Channel, Order);
    for (m=0; m<M; m++)
        gsl_matrix_set(L, m, X_, B[m]);

    // (2) compute field for Qy'
    for (v=X_; v<=Z_; v++)
        Q6.v[v] = Vy[v];
    FwdVec(&Q6, B, Header, Channel, Order);
    for (m=0; m<M; m++)
        gsl_matrix_set(L, m, Y_, B[m]);

    free(B);
}
