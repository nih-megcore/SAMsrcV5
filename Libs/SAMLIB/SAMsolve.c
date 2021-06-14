// SAMsolve() -- all-in-one SAM solution
//
//  Author: Stephen E. Robinson
//      MEG Core Facility
//      NIMH
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>

void SAMsolve(
    gsl_matrix  *C,         // C - MxM unaugmented covariance matrix (input)
    gsl_matrix  *Cinv,      // Cinv - MxM inverse covariance matrix (input)
    gsl_vector  *B,         // B - Mx1 forward solution vector (input)
    gsl_vector  *Wgt,       // W - M SAM coefficients (output)
    double      noise,      // mean sensor noise (input)
    double      *s2,        // SAM source power (output)
    double      *n2         // SAM noise power (output)
) {
    gsl_vector  *A;
    double      tmp;
    double      norm;
    int         M;

    // allocate gsl memory
    M = C->size1;
    A = gsl_vector_alloc(M);

    // solve Wgt = C^-1 B / B' C^-1 B
    gsl_blas_dgemv(CblasNoTrans, 1., Cinv, B, 0., Wgt);     // Wgt (unnormalized) = C(-1) B -- MxM x Mx1 = Mx1
    gsl_blas_ddot(B, Wgt, &tmp);                            // tmp = B' Wgt -- M x M
    norm = 1. / tmp;
    gsl_vector_scale(Wgt, norm);                            // Wgt /= norm

    // solve source power & noise power-- s2 = Wgt' C Wgt
    gsl_blas_dgemv(CblasNoTrans, 1., C, Wgt, 0., A);        // A = C Wgt = MxM x Mx1 = Mx1
    gsl_blas_ddot(Wgt, A, s2);                              // s2 = Wgt' A -- Mx1 x Mx1
    gsl_blas_ddot(Wgt, Wgt, n2);                            // n2 (unscaled) = Wgt' Wgt = Mx1 x Mx1
    *n2 *= noise;                                           // n2 *= noise

    gsl_vector_free(A);
}
