// FitDip() -- fit one magnetic dipole to measurements from an array of OPM sensors
//
//  Author:	Stephen E. Robinson
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
#include <samlib.h>
#include <geoms.h>
#include <jig.h>

// #define DEBUG   TRUE

double  FitDip(
    int         D,                  // dipole coil index
    gsl_matrix  *Measured,          // calibrator measurements - S x D (input)
    gsl_matrix  *Variance,          // measured sensor noise variance - S x D (input)
    double      *Current,           // Current[N_DIP] - rms coil current (A) or moment (A-m^2) (input)
    FIELD       *Sensor,            // Sensor[N_OPM] - array of OPM sensors (input)
    FIELD       *Dipole             // Dipole[N_DIP] - array of magnetic dipoles (initial guess & output)
)
{
    gsl_matrix      *covar;         // covar - 6x6 augmented covariance of partial derivatives
    gsl_matrix      *A;             // copy of covariance for LU decomposition
    gsl_matrix      *alpha;         // alpha - 6x6 partial derivatives
    gsl_vector      *beta;          // beta - 6x1 vector
    gsl_vector      *da;            // input vector - 6x1
    gsl_vector      *x;             // solution vector - 6x1
    gsl_vector      *a;             // linear OPM parameters - 6x1
    gsl_vector      *Try;           // trial parameters - 6x1
    gsl_permutation *p;             // LU permutation
    double          lambda;         // covariance augmentation factor
    double          chi2;           // chi-square
    double          Ochi2;          // previous chi2
    int             signum;         // sign of LU decomposition
    int             its;            // iterations
    int             done = FALSE;

    // gsl allocations
    covar = gsl_matrix_alloc(6, 6);
    A = gsl_matrix_alloc(6, 6);
    alpha = gsl_matrix_alloc(6, 6);
    beta = gsl_vector_alloc(6);
    da = gsl_vector_alloc(6);
    x = gsl_vector_alloc(6);
    a = gsl_vector_alloc(6);
    Try = gsl_vector_alloc(6);
    p = gsl_permutation_alloc(6);

    // convert parameters of single dipole to vector
    gsl_vector_norm(a);                             // make orientation a unit vector
    FIELD2Vector(&Dipole[D], a);
    chi2 = DipCoeff(D, Measured, Variance, Current, Sensor, a, alpha, beta);    // solve goodness of fit (chi2) for initial parameters
    Ochi2 = chi2;

    // iterate using perturbed parameters until lambda is tiny
    its = 0;
    lambda = INIT_LAMBDA;
    do {

        // alter linearized fitting matrix by augmenting diagonal elements
        gsl_matrix_memcpy(covar, alpha);            // alpha --> covar
        gsl_matrix_mul_diag(covar, 1. + lambda);    // augment diagonal
        gsl_vector_memcpy(da, beta);                // beta --> da

        // matrix solution for perturbed trial vector
        gsl_matrix_memcpy(A, covar);                // covar --> A
        gsl_linalg_LU_decomp(A, p, &signum);
        gsl_linalg_LU_solve(A, p, da, x);

        // Try = a + x
        gsl_vector_memcpy(Try, a);                  // a --> Try
        gsl_vector_add(Try, x);                     // Try += x

        // did the trial succeed?
        chi2 = DipCoeff(D, Measured, Variance, Current, Sensor, Try, covar, da);    // test trial parameters
        if (chi2 > Ochi2) {                         // failure -- increase lambda & return
            lambda *= LAMBDA_UP;
        } else {                                    // success -- accept new solution
            lambda /= LAMBDA_DN;
            Ochi2 = chi2;
            gsl_matrix_memcpy(alpha, covar);        // covar --> alpha
            gsl_vector_memcpy(beta, da);            // da --> beta
            gsl_vector_memcpy(a, Try);              // Try --> a
        }
        its++;
        if (its > MAXITS || lambda > MAX_LAMBDA || lambda < END_LAMBDA)
            done = TRUE;;

    } while(!done);

    // return dipole parameters in type FIELD
    gsl_vector_norm(a);
    Vector2FIELD(a, &Dipole[D]);

    // free gsl allocations
    gsl_matrix_free(covar);
    gsl_matrix_free(A);
    gsl_matrix_free(alpha);
    gsl_vector_free(beta);
    gsl_vector_free(da);
    gsl_vector_free(x);
    gsl_vector_free(a);
    gsl_vector_free(Try);
    gsl_permutation_free(p);

    // return goodness of fit
    return(chi2);
}


// compute coefficients for one dipole
double  DipCoeff(                   // return chi-square
    int         D,                  // dipole coil index
    gsl_matrix  *Measured,          // OPM sensor measurements - S x D (input)
    gsl_matrix  *Variance,          // measured noise variance (sigma^2) for each OPM sensor - S x D (input)
    double      *Current,           // rms coil current (A) or moment (A-m^2) (input)
    FIELD       *Sensor,            // Sensor[S] - array of OPM sensors (input)
    gsl_vector  *a,                 // linear dipole parameter vector (input)
    gsl_matrix  *alpha,
    gsl_vector  *beta
)
{
    FIELD       Dipole;             // test dipole parameters
    gsl_vector  *Try;               // test parameters - 6x1
    double      dyda[6];            // partial derivatives
    double      dy[6];              // difference between measured and trial predicted
    double      measured;           // measured field
    double      predicted;          // predicted field
    double      sigma2;             // variance
    double      dY;                 // difference between measured & predicted
    double      chi2;               // chi-square
    double      wt;                 // dyda / sigma2
    double      tmp;
    double      delta = 1.0e-05;    // perturbations for position and orientation
    int         S;                  // number of OPM sensors
    int         i;
    int         j;
    int         s;

    // extract sensor array size
    S = Measured->size1;

    // allocate deltas
    Try = gsl_vector_alloc(6);

	// initialize (symmetric) alpha & beta
	gsl_matrix_set_zero(alpha);
	gsl_vector_set_zero(beta);

    // create Dipole with unit vector
    gsl_vector_norm(a);
    Vector2FIELD(a, &Dipole);

    // sum over all magnetic sensor locations
    for (s=0, chi2=0.; s<S; s++) {

        // compute difference between predicted & measured field & accumulate chi2
        measured = gsl_matrix_get(Measured, s, D);
#ifdef MTOB
        predicted = MtoB(&Dipole, &Sensor[s], Current[D]);
#else
        predicted = Coil2B(&Dipole, &Sensor[s], Current[D]);
#endif
        sigma2 = gsl_matrix_get(Variance, s, D);
        dY = measured - predicted;

        // compute dy/da
        for (i=0; i<6; i++) {

            // compute delta field
            for (j=0; j<6; j++)
                if (i == j)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) + delta);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));
            Vector2FIELD(Try, &Dipole);
#ifdef MTOB
            predicted = MtoB(&Dipole, &Sensor[s], Current[D]);
#else
            predicted = Coil2B(&Dipole, &Sensor[s], Current[D]);
#endif
            // compute partial derivative
            dy[i] = measured - predicted;
            dyda[i] = dy[i] / delta;
        }

        // accumulate alpha & beta (comments: J: J=dyda, W=1/sigma2)
        for (i=0; i<6; i++) {
            for (j=0; j<6; j++) {
                wt = dyda[i] / sigma2;                              // J'W
                tmp = gsl_matrix_get(alpha, i, j) + wt * dyda[j];   // J'WJ
                gsl_matrix_set(alpha, i, j, tmp);
            }
            tmp = gsl_vector_get(beta, i) + wt * dy[i];             // J'W dy
            gsl_vector_set(beta, i, tmp);
        }

        // scale alpha & beta
        for (i=0; i<6; i++) {
            gsl_matrix_scale(alpha, 0.01);
            gsl_vector_scale(beta, 0.01);
        }

        // accumulate chi2
        chi2 += dY * dY / sigma2;
    }
    gsl_vector_free(Try);

    // return chi-square
    return(chi2 / (double)S);
}
