// SolveMoment() -- solves eigensystem for moment vector associated with
//  maximum eigenvalue.
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <lfulib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

void SolveMoment(
    gsl_matrix      *Cinv,          // Cinv - NxN C^-1
    gsl_matrix      *L,             // L - Mx3 lead field matrix
    double          *Vec,           // Vec[M] - moment vector
    double          *span           // span of eigensystem solution
) {
    gsl_eigen_gensymmv_workspace    *w; // generalized eigensystem workspace
    gsl_matrix      *K1;            // L' C^(-1) L -- source power
    gsl_matrix      *K2;            // L' C^(-2) L -- noise power (unscaled)
    gsl_matrix      *A;             // intermediate matrix
    gsl_matrix      *evec;          // eigenvectors
    gsl_vector      *eval;          // eigenvalues
    gsl_vector      *vec;           // moment vector
    gsl_matrix      *Cinv2;         // Cinv^2
    int             v;              // vector index
    int             M;              // number of primary sensors

    // compute C^-2
    M = Cinv->size1;
    Cinv2 = gsl_matrix_alloc(M, M);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv, Cinv, 0., Cinv2);

    // allocate gsl vectors, matrices & workspace
    A = gsl_matrix_alloc(M, 3);
    K1 = gsl_matrix_alloc(3, 3);
    K2 = gsl_matrix_alloc(3, 3);
    evec = gsl_matrix_alloc(3, 3);
    eval = gsl_vector_alloc(3);
    vec = gsl_vector_alloc(3);
    w = gsl_eigen_gensymmv_alloc(3);

    // compute K1 = L' C^(-1) L
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv, L, 0., A);     // A = C^(-1) * L -- NxN * NxQV = NxQV
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., L, A, 0., K1);         // K1 = L' * A -- QVxN * NxQV = QVxQV

    // compute K2 = L' C^(-2) L
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., Cinv2, L, 0., A);    // A = C^(-2) * L -- NxN * NxQV = NxQV
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1., L, A, 0., K2);         // K2 = L' * A -- QVxN * NxQV = QVxQV

    // solve generalized eigensystem for eigenvector corresponding to maximum eigenvalue
    if (gsl_eigen_gensymmv(K1, K2, eval, evec, w) == GSL_SUCCESS) {
        gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
        gsl_matrix_get_col(vec, evec, 0);   // get eigenvector for maximum eigenvalue -- source moment vector
        for (v=X_; v<=Z_; v++)
            Vec[v] = gsl_vector_get(vec, v);
        *span = gsl_vector_get(eval, 0) / gsl_vector_get(eval, 2);
    } else {                                // eigensystem solution failed -- solution is zero
        for (v=X_; v<=Z_; v++)
            Vec[v] = 0.;
    }

    // free gsl vectors, matrices & workspace
    gsl_matrix_free(A);
    gsl_matrix_free(K1);
    gsl_matrix_free(K2);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_vector_free(vec);
    gsl_eigen_gensymmv_free(w);
    gsl_matrix_free(Cinv2);
}
