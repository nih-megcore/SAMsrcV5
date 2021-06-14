// NSolvewts() - solve for independent beamformer weights
//	Forward solutions are orthogonalized. Call NSolveFwd()
//  for computing Bfwd
//
//  WARNING!: Do not call this subroutine multiple times
//  without first freeing the gsl weight matrix!
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <stdlib.h>
#include <lfulib.h>
#include <siglib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

void    NSolveWts(
    gsl_matrix  *W,                 // W - MxN beamformer coefficient matrix (externally allocated)
	gsl_matrix	*Bfwd,              // Bfwd - MxN N dipole forward solutions (input)
	gsl_matrix	*C                  // C - MxM solution covariance matrix (input)
)

{
	gsl_matrix	        *Cinv;		// Cinv - MxM inverse of solution covariance matrix
	gsl_matrix	        *K1;		// K1 - MxN K1 = C^-1 B
	gsl_matrix	        *K2;		// K2 - NxN K2 = (B' C^-1 B)^-1 = (B' K1)^-1
	gsl_matrix	        *K2inv;		// K2inv - NxN denominator matrix
	gsl_vector	        *Bt;        // Bt - Mx1 test forward solution
	int			        i;			// i-index
	int			        nzeros;		// number of zero forward solutions
	int			        M;			// number of primary sensors
	int			        N;			// number of sources

	// check dimensions
	if(C->size1 != Bfwd->size1)
		cleanup("number of sensors in B != dimension of C");
	M = Bfwd->size1;
	N = Bfwd->size2;

	// allocate W, Cinv & invert covariance
	Cinv = gsl_matrix_alloc(M, M);
	pinv(C, Cinv);

	// allocate solution matrices
	Bt = gsl_vector_alloc(M);
	K1 = gsl_matrix_alloc(M, N);
	K2 = gsl_matrix_alloc(N, N);
	K2inv = gsl_matrix_alloc(N, N);

	// check if any of the N-forward solution is zero
	for(i=nzeros=0; i<N; i++) {
		gsl_matrix_get_col(Bt, Bfwd, i);
		if(gsl_vector_isnull(Bt) == 1)
			nzeros++;
	}
	if(nzeros == 0) {

		// compute K1 = C^-1 B
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv, Bfwd, 0., K1);		// K1 = C^(-1) * Bfwd -- MxM x MxN = MxN

		// compute K2 = (B' C^-1 B)^-1 = (B' K1)^-1
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., Bfwd, K1, 0., K2);			// K2 = Bfwd' * K1 -- NxM x MxN = NxN

		// invert K2
		gsl_matrix_memcpy(K2inv, K2);
        gsl_linalg_cholesky_decomp(K2inv);
        gsl_linalg_cholesky_invert(K2inv);

		// compute W = K1 * K2^-1
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., K1, K2inv, 0., W);

	} else {	// zero weights
		gsl_matrix_set_zero(W);
	}

	// free matrices
	gsl_matrix_free(K1);
	gsl_matrix_free(K2);
	gsl_matrix_free(K2inv);
	gsl_matrix_free(Cinv);
	gsl_vector_free(Bt);
}
