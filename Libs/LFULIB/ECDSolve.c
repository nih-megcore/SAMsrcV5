// ECDSolve() - solve for beamformer optimal forward solution
//  for current dipole in realistic (Nolte) model. See DIPSolve.c
//  for the multisphere & single sphere models, respectively.
//
//  Author: Stephen E. Robinson
//      MEG Core Facility
//      NIMH
//

#include <stdlib.h>
#include <math.h>
#include <lfulib.h>
#include <voxel.h>
#include <siglib.h>

void ECDSolve(
    VOXELINFO       *Voxel,     // dipole position, orientation, & flags
    gsl_vector      *B,         // B M optimal target forward solution
    gsl_matrix      *Cinv,      // Cinv - MxM inverse covariance matrix
    HeaderInfo      *Header,    // header information with sensor list
    ChannelInfo     *Channel,   // sensor information
    HULL            *Hull,      // cortical hull model
    double          *span       // eigenvalue span of 3x3 solution (return value)
) {
    gsl_matrix      *L;         // L - MxQV primary sensor lead field
    double          Vec[3];     // Vec[3] - moment vector
    double          tmp;        // accumulator
    int             m;          // sensor index
    int             v;          // vector index
    int             M;          // primary sensor rank

    // allocate memory
    M = Cinv->size1;
    L = gsl_matrix_alloc(M, 3);

    // compute ECD lead field matrix - compute forward solution matrix for unit dipole terms
#if ELEKTA
    ECDLeadFieldElekta(Voxel->p, L. Header, Channel, Hull);
#else
    ECDLeadField(Voxel->p, L, Header, Channel, Hull);
#endif

    // if req'd, solve generalized eigensystem for moment vector
    if (Voxel->Solve) {
        SolveMoment(Cinv, L, Vec, span);
        for (v=X_; v<=Z_; v++)      // replace voxel moment vector w/ Vec
            Voxel->v[v] = Vec[v];
    }

    // compute forward solution
    for (m=0; m<M; m++) {
        for (v=X_, tmp=0.; v<=Z_; v++)
            tmp += Voxel->v[v] * gsl_matrix_get(L, m, v);
        gsl_vector_set(B, m, tmp);
    }

    // free allocations
    gsl_matrix_free(L);
}

// temporary version until everything is converted over to multithread

void ECDSolve0(
    VOXELINFO       *Voxel,     // dipole position, orientation, & flags
    gsl_vector      *B,         // B M optimal target forward solution
    gsl_matrix      *Cinv,      // Cinv - MxM inverse covariance matrix
    HeaderInfo      *Header,    // header information with sensor list
    ChannelInfo     *Channel,   // sensor information
    HULL            *Hull,      // cortical hull model
    double          *span,      // eigenvalue span of 3x3 solution (return value)
    COEFFS          *coeff      // integration points and basis coefficients
) {
    gsl_matrix      *L;         // L - MxQV primary sensor lead field
    double          Vec[3];     // Vec[3] - moment vector
    double          tmp;        // accumulator
    int             m;          // sensor index
    int             v;          // vector index
    int             M;          // primary sensor rank

    // allocate memory
    M = Cinv->size1;
    L = gsl_matrix_alloc(M, 3);

    // compute ECD lead field matrix - compute forward solution matrix for unit dipole terms
#if ELEKTA
    ECDLeadFieldElekta(Voxel->p, L. Header, Channel, Hull);
#else
    ECDLeadField0(Voxel->p, L, Header, Channel, Hull, coeff);
#endif

    // if req'd, solve generalized eigensystem for moment vector
    if (Voxel->Solve) {
        SolveMoment(Cinv, L, Vec, span);
        for (v=X_; v<=Z_; v++)      // replace voxel moment vector w/ Vec
            Voxel->v[v] = Vec[v];
    }

    // compute forward solution
    for (m=0; m<M; m++) {
        for (v=X_, tmp=0.; v<=Z_; v++)
            tmp += Voxel->v[v] * gsl_matrix_get(L, m, v);
        gsl_vector_set(B, m, tmp);
    }

    // free allocations
    gsl_matrix_free(L);
}
