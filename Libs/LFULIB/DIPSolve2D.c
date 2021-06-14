// DIPSolve2D() - solve for beamformer optimal forward solution
//  for current dipole in single local sphere model (i.e. two
//  dipole vector components
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <model.h>
#include <lfulib.h>


void    DIPSolve2D(
    VOXELINFO       *Voxel,     // dipole position, orientation, & flags
    gsl_vector      *Bp,        // Mx1 target forward solution
    gsl_matrix      *Cinv,      // MxM inverse covariance matrix
    HeaderInfo      *Header,    // header information with sensor list
    ChannelInfo     *Channel,   // sensor information
    int             Order,      // field integration order
    double          *span       // eigenvalue span of 2x2 solution (return value)
)
{
    gsl_matrix  *L;             // Mx2 lead-field matrix
    gsl_matrix  *K1;            // Mx2 L' C^(-1) L -- source power
    gsl_matrix  *K2;            // Mx2 L' C^(-2) L -- noise power (unscaled)
    gsl_matrix  *A;             // Mx2 intermediate matrix
    gsl_matrix  *Cinv2;         // MxM Cinv^2
    gsl_matrix  *eVec;          // 2x2 eigenvectors
	gsl_vector  *eVal;          // 2x1 eigenvalues
    double      tmp;            // accumulator
    int         M;              // number of primary sensors
    int         i;              // i-index
    int         j;              // j-index

    // allocate memory
    M = Cinv->size1;
    L = gsl_matrix_alloc(M, 2);
    Cinv2 = gsl_matrix_alloc(M, M);
    K1 = gsl_matrix_alloc(2, 2);
    K2 = gsl_matrix_alloc(2, 2);
    A = gsl_matrix_alloc(M, 2);
    eVec = gsl_matrix_alloc(2, 2);
    eVal = gsl_vector_alloc(2);

    // compute C^-2
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv, Cinv, 0., Cinv2);

    // compute DIP lead field matrix - compute forward solution matrix for unit dipole terms
    DIPLeadField2D(Voxel->p, L, Header, Channel, Order);
    
    // compute K1 = L' C^(-1) L
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv, L, 0., A);     // A = C^(-1) * L -- MxM * Mx2 = Mx2
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., L, A, 0., K1);         // K1 = L' * A -- 2xM * Mx2 = 2x2

    // compute K2 = L' C^(-2) L
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv2, L, 0., A);    // A = C^(-2) * L -- MxM * Mx2 = Mx2
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., L, A, 0., K2);         // K2 = L' * A -- 2xM * Mx2 = 2x2

	// solve generalized eigensystem problem for SNR
    geig22(K1, K2, eVec, eVal);
    *span = gsl_vector_get(eVal, 0) / gsl_vector_get(eVal, 1);

	// compute forward solution
	for(i=0; i<M; i++) {
		for(j=0, tmp=0.; j<2; j++)
			tmp += gsl_matrix_get(eVec, 0, j) * gsl_matrix_get(L, i, j);    // eigenvectors are row vectors in the 'geig22()' solution
        gsl_vector_set(Bp, i, tmp);
    }

    gsl_matrix_free(L);
    gsl_matrix_free(Cinv2);
    gsl_matrix_free(K1);
    gsl_matrix_free(K2);
    gsl_matrix_free(A);
    gsl_matrix_free(eVec);
    gsl_vector_free(eVal);
}
