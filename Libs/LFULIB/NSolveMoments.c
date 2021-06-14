// NSolveMoments() - solve for moment vector for any
//	N current dipole in realistic head model.
//
//  WARNING!: Do not call this subroutine multiple times
//  without first freeing the gsl weight matrix!
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <siglib.h>
#include <lfulib.h>
#include <samlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <voxel.h>


void    NSolveMoments(
    gsl_matrix      *Vec,                       // Vec Nx3 source unit moment vectors (externally allocated)
    gsl_matrix      *Cinv,                      // Cinv - MxM C^-1 (input)
    N_LEADFIELD     *Lfld,                      // L[N] - array of Mx3 lead field matrix (input)
    int             N,                          // number of sources (input)
    double          *span                       // span of eigensystem solution (output)
)

{
	gsl_eigen_gensymmv_workspace	*w;         // generalized eigensystem workspace
	gsl_matrix                      *K1;        // source power matrix (numerator)
	gsl_matrix                      *K2;        // noise power matrix (denominator)
	gsl_vector                      *eval;      // eigenvalues
	gsl_matrix                      *evec;      // eigenvectors
	gsl_matrix                      *L;         // L - MxQV lead-fields
	gsl_matrix                      *A;         // A - MxQV intermediate matrix
	gsl_matrix                      *Cinv2;     // Cinv2 - MxM inverse orientation covariance matrix squared
	gsl_vector                      *xyz;       // xyz - 3x1 vector
	gsl_vector                      *mom;       // mom - QVx1 unnormalized moment vector
	gsl_vector_view                 vn;         // vn - VSizex1 moment vector (this a stack variable & not a pointer)
	double                          r;          // eigenvector normalization constant
	double                          tmp;        // accumulator
	int                             i;          // i-index
	int                             n;          // target index
	int                             m;          // sensor index
	int                             v;          // vector index
    int                             M;          // number of primary sensors
	int                             QV;         // number of components (3N)

    // get dimensions
    M = Cinv->size1;
    QV = 3 * N;

    // allocate memory & invert covariance
    Cinv2 = gsl_matrix_alloc(M, M);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv, Cinv, 0., Cinv2);
	K1 = gsl_matrix_alloc(QV, QV);
	K2 = gsl_matrix_alloc(QV, QV);
	L = gsl_matrix_alloc(M, QV);
	A = gsl_matrix_alloc(M, QV);
	evec = gsl_matrix_alloc(QV, QV);
	eval = gsl_vector_alloc(QV);
	w = gsl_eigen_gensymmv_alloc(QV);
	mom = gsl_vector_alloc(QV);
	xyz = gsl_vector_alloc(3);

    // assemble Mx3N lead-field matrix in order of n1-xyz, n2-xyz, ...N-1-xyz
    for(m=0; m<M; m++)
        for(n=0; n<N; n++) {
            gsl_matrix_get_row(xyz, Lfld[n].L, m);      // get the mth row of L(n) - a 3x1 vector
            for(v=X_; v<=Z_; v++) {
                i = n * 3 + v;
                gsl_matrix_set(L, m, i, gsl_vector_get(xyz, v));
            }
        }

	// compute numerator matrix -- 1/S^2: K1 = L' C^(-1) L
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv, L, 0., A);		// A = M*M * M*QV = MxQV
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., L, A, 0., K1);			// K1 = QVxM * M*QV = QVxQV

	// compute denominator matrix -- 1/N^2 (without sigma^2): K2 = L' C^(-2) L
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv2, L, 0., A);	// A = M*M * M*QV = MxQV
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., L, A, 0., K2);			// K2 = QVxM * M*QV = QVxQV

	// solve generalized eigensystem for eigenvector corresponding to maximum eigenvalue
	if(gsl_eigen_gensymmv(K1, K2, eval, evec, w) == GSL_SUCCESS) {

        // sort eigenvalues, get span, & select eigenvector
		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
        *span = gsl_vector_get(eval, 0) / gsl_vector_get(eval, QV-1);
		gsl_matrix_get_col(mom, evec, 0);	// get eigenvector -- source unnormalized moment vectors for N sources

		// parse N moment vectors & independently normalize them
		for(n=0; n<N; n++) {

            // parse mom into individual moment vectors vn
            vn = gsl_vector_subvector(mom, 3*n, 3);

            // normalize vn
            for(v=X_, tmp=0.; v<=Z_; v++)
                tmp += gsl_vector_get(&vn.vector, v) * gsl_vector_get(&vn.vector, v);
            r = 1. / sqrt(tmp);
            gsl_vector_scale(&vn.vector, r);

            // save vn in Vec
            gsl_matrix_set_row(Vec, n, &vn.vector);
        }

	} else {	// eigensystem solution failed, set all moments to zero
	    gsl_matrix_set_zero(Vec);
	    *span = 0.;
	}

	// free gsl vectors, matrices & workspace
	gsl_matrix_free(Cinv2);
	gsl_matrix_free(K1);
	gsl_matrix_free(K2);
	gsl_matrix_free(A);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);
	gsl_vector_free(mom);
	gsl_vector_free(xyz);
	gsl_eigen_gensymmv_free(w);
}
