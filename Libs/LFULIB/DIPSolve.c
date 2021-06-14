// DIPSolve() - solve for beamformer optimal forward solution
//  for current dipole in single or multiple local spheres model.
//  Note that individual sensor local sphere origins are set in
//  the calling program & may be single or multiple origins.
//
//  Author: Stephen E. Robinson
//      MEG Core Facility
//      NIMH
//

#include <stdlib.h>
#include <math.h>
#include <model.h>
#include <geoms.h>
#include <voxel.h>
#include <samlib.h>
#include <lfulib.h>
#include <siglib.h>

void DIPSolve(
    VOXELINFO       *Voxel,     // dipole position, orientation, & flags
    gsl_vector      *Bp,        // Bp - Mx1 optimal target forward solution
    gsl_matrix      *Cinv,      // Cinv - MxM inverse covariance matrix
    HeaderInfo      *Header,    // header information with sensor list
    ChannelInfo     *Channel,   // sensor information
    int             Order,      // field integration order
    double          *span       // eigenvalue span of 3x3 solution (return value)
) {
    gsl_matrix      *L;         // L- MxQV primary sensor lead field
    double          Vec[3];     // moment vector
    double          tmp;        // you don't want to know what this is...
    int             m;          // sensor index
    int             v;          // vector index
    int             M;          // primary sensor rank

    // allocate memory
    M = Cinv->size1;
    L = gsl_matrix_alloc(M, 3);

    // compute DIP lead field matrix - compute forward solution matrix for unit dipole terms
    DIPLeadField(Voxel->p, L, Header, Channel, Order);

    // linear solution of dipole moment vector
    if (Voxel->Solve) {

        // solve generalized eigensystem for moment vector
        SolveMoment(Cinv, L, Vec, span);

        // copy moment vector to voxel
        for (v=X_; v<=Z_; v++)
            Voxel->v[v] = Vec[v];
    } else {
        *span = 0.;
    }

    // compute dipole forward solution
    for (m=0; m<M; m++) {
        for (v=X_, tmp=0.; v<=Z_; v++)
            tmp += Voxel->v[v] * gsl_matrix_get(L, m, v);
        gsl_vector_set(Bp, m, tmp);
    }

    // free allocations
    gsl_matrix_free(L);
}
