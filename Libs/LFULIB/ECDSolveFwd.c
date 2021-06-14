// ECDSolve() - solve for constrained & unconstrained forward solutions
//  plus unconstrained moment vector
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdlib.h>
#include <lfulib.h>
#include <samlib.h>
#include <siglib.h>
#include <math.h>
#include <voxel.h>

#define QV      3           // 3 dipole moment solution

// Convenience function to call when the sensor frame hasn't changed.

void ECDSolveFwd(
    VOXELINFO   *Voxel,     // dipole position, orientation, & flags
    gsl_vector  *Bc,        // M anatomically constrained forward solution
    gsl_vector  *Bu,        // M unconstrained forward solution
    gsl_matrix  *C,         // MxM covariance matrix
    gsl_matrix  *Cinv,      // MxM inverse covariance matrix
    HeaderInfo  *Header,    // header information with sensor list
    ChannelInfo *Channel,   // sensor information
    HULL        *Hull,      // cortical hull model
    double      *span       // eigenvalue span of 3x3 solution (return value)
) {
    ECDSolveFwd0(Voxel, Bc, Bu, C, Cinv, Header, Channel, Hull, span, NULL);
}

void ECDSolveFwd0(
    VOXELINFO   *Voxel,     // dipole position, orientation, & flags
    gsl_vector  *Bc,        // M anatomically constrained forward solution
    gsl_vector  *Bu,        // M unconstrained forward solution
    gsl_matrix  *C,         // MxM covariance matrix
    gsl_matrix  *Cinv,      // MxM inverse covariance matrix
    HeaderInfo  *Header,    // header information with sensor list
    ChannelInfo *Channel,   // sensor information
    HULL        *Hull,      // cortical hull model
    double      *span,      // eigenvalue span of 3x3 solution (return value)
    COEFFS      *coeff      // COEFFS or NULL to indicate which ECDLeadField to call
) {
    gsl_matrix  *L;         // MxQV primary sensor lead field
    double      Vec[QV];    // moment vector
    double      tmp;        // accumulator
    int         i, j;       // indices
    int         M;          // primary sensor rank

    // allocate memory
    M = C->size1;
    L = gsl_matrix_alloc(M, QV);

    // compute ECD lead field matrix -
    // compute forward solution matrix for unit (1. A-m) dipole in each cardinal axis
#if ELEKTA
    ECDLeadFieldElekta(Voxel->p, L, Header, Channel, Hull);
#else
    if (coeff) {
        ECDLeadField0(Voxel->p, L, Header, Channel, Hull, coeff);
    } else {
        ECDLeadField(Voxel->p, L, Header, Channel, Hull);
    }
#endif

    // compute constrained forward solution
    for (i=0; i<M; i++) {
        for (j=0, tmp=0.; j<QV; j++)
            tmp += Voxel->v[j] * gsl_matrix_get(L, i, j);
        gsl_vector_set(Bc, i, tmp);
    }

    // solve generalized eigensystem for unconstrained moment vector
    SolveMoment(Cinv, L, Vec, span);
    for (i=0; i<QV; i++)
        Voxel->v[i] = Vec[i];       // make sure to return unconstrained moment vector!

    // compute unconstrained forward solution
    for (i=0; i<M; i++) {
        for (j=0, tmp=0.; j<QV; j++)
            tmp += Voxel->v[j] * gsl_matrix_get(L, i, j);
        gsl_vector_set(Bu, i, tmp);
    }

    gsl_matrix_free(L);
}
