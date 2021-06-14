// pinv() -- pseudoinverse of square matrix
//
// Author:  Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "siglib.h"
#include "samutil.h"

void pinv(
    gsl_matrix  *A,         // matrix A (input)
    gsl_matrix  *Ainv       // inverse matrix A (output)
) {
    pinvN(A, Ainv, 0);
}

void pinvN(
    gsl_matrix  *A,         // matrix A (input)
    gsl_matrix  *Ainv,      // inverse matrix A (output)
    int ndim                // number of small singular values to remove
) {
    gsl_matrix  *U;
    gsl_matrix  *V;
    gsl_vector  *S;
    gsl_vector  *Work;
    gsl_matrix  *Sinv;
    gsl_matrix  *VSinv;
    double      min;
    int         N;          // input matrix, full rank
    int         i;

    // get matrix dimensions
    if (A->size1 != A->size2)
        cleanup("matrix must be square");
    N = A->size1;

    // allocate matrices & vectors for SVD
    U = gsl_matrix_alloc(N, N);
    V = gsl_matrix_alloc(N, N);
    S = gsl_vector_alloc(N);
    Work = gsl_vector_alloc(N);

    // copy A to U & compute SVD
    gsl_matrix_memcpy(U, A);
    gsl_linalg_SV_decomp(U, V, S, Work);
    gsl_vector_free(Work);

    // invert S & reduce rank to math precision
    Sinv = gsl_matrix_calloc(N, N);
    min = gsl_vector_get(S, 0) * 1.0e-15;
    for (i=0; i<N; i++) {
        if (gsl_vector_get(S, i) > min)
            gsl_matrix_set(Sinv, i, i, 1. / gsl_vector_get(S, i));
    }

    // further reduce rank if requested
    if (0 < ndim && ndim < N) {
msg("[pinv] truncating for ndim at %g\n", gsl_vector_get(S, N - ndim));
        for (i=N-ndim; i<N; i++)
            gsl_matrix_set(Sinv, i, i, 0.);
    }
    gsl_vector_free(S);

    // compute pseudoinverse -- Ainv = V S(-1) U'
    VSinv = gsl_matrix_alloc(N, N);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Sinv, 0., VSinv);     // VSinv = V Sinv
    gsl_matrix_free(Sinv);
    gsl_matrix_free(V);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., VSinv, U, 0., Ainv);       // Ainv = V Sinv U'
    gsl_matrix_free(VSinv);
    gsl_matrix_free(U);
}
