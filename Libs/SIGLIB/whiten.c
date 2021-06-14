// whiten() -- use eigendecomposition of covariance to whiten data
//
//  Author: Stephen E. Robinson
//      MEG Core Facility
//      NIMH
//

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <samlib.h>


void whiten(
    double  **x,            // x[N][T] -- input data
    double  **y,            // y[N][T] -- output (whitened) data
    int     N,              // matrix rank
    int     *K,             // effective data rank
    int     T               // number of samples
)
{
    gsl_eigen_symm_workspace *w;    // eigensystem workspace
    gsl_matrix  *C;         // C - NxN covariance matrix (and eigenvectors)
    gsl_vector  *eval;      // eigenvalues (unsorted)
    gsl_vector  *sval;      // sorted eigenvalues
    double      **S;        // S - NxN sphering matrix
    double      delta;
    double      min;
    double      sum;        // integrator
    double      tmp;        // what else could this be?
    int         i;
    int         j;
    int         k;
    int         t;          // sample index

    // allocate memory
    C = gsl_matrix_alloc(N, N);
    eval = gsl_vector_alloc(N);
    sval = gsl_vector_alloc(N);
    w = gsl_eigen_symm_alloc(N);
    if((S = (double **)malloc(N * sizeof(double *))) == NULL)
        allocfailed("S[]");
    for(i=0; i<N; i++)
        if((S[i] = (double *)malloc(N * sizeof(double))) == NULL)
            allocfailed("S[][]");

    // compute covariance of x
    for(i=0; i<N; i++)
        for(j=0; j<N; j++) {
            for(t=0, sum=0.; t<T; t++)
                sum += x[i][t] * x[j][t];
            sum /= (double)T;
            gsl_matrix_set(C, i, j, sum);
        }

    // compute eigendecomposition of covariance
    gsl_eigen_symm(C, eval, w);
    gsl_vector_memcpy(sval, eval);

    // look at derivatives - find K for minimum derivative
    gsl_sort_vector(sval);      // eigenvalues are now in ascending order
    for(k=2, min=HUGE; k<N-3; k++) {
        delta = gsl_vector_get(sval, k+2) - gsl_vector_get(sval, k-2);
        if (delta < min) {
            min = delta;
            j = k;
        }
    }
    *K = j;     // output effective rank

    // compute sphering matrix using unsorted eigenvalues & eigenvectors
    for(i=0; i<N; i++) {
        tmp = sqrt(gsl_vector_get(eval, i));
        for(j=0; j<N; j++)
            S[i][j] = gsl_matrix_get(C, i, j) / tmp;
    }

    // compute y(t) = L^(-1) x(t)
    for(t=0; t<T; t++) {
        for(i=0; i<N; i++) {
            for(j=0, sum=0.; j<N; j++) 
                sum += S[i][j] * x[j][t];
            y[i][t] = sum;
// printf("x[%d][%d]=%e, y[%d][%d]=%e\n", i, t, x[i][t], i, t, y[i][t]);
        }
    }

    // free allocations
    gsl_matrix_free(C);
    gsl_vector_free(eval);
    gsl_vector_free(sval);
    gsl_eigen_symm_free(w);
    for(i=0; i<N; i++)
        free((void *)S[i]);
    free((void *)S);
}
