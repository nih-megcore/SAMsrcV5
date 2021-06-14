// FitOPM() -- fit OPM a single sensor to measurements from an array of magnetic dipoles
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

#define P       7   //  7-parameter solution
#define GAIN_   6
// #define DEBUG   TRUE

double  FitOPM(
    int         S,                  // OPM sensor index
    gsl_matrix  *Measured,          // sensor measurements - S x D (input)
    gsl_matrix  *Variance,          // measured noise variance - S x D (input)
    double      *Current,           // rms coil drive current (A)
    double      *Scale,             // scale factor (sensor gain multiplier)
    FIELD       *Sensor,            // Sensor[S] - OPM sensor array parameters position & orientation (initial guess & outout)
    FIELD       *Dipole             // Dipole[D] - magnetic dipoles array (input)
)
{
    gsl_matrix      *covar;         // covar - PxP augmented covariance of partial derivatives
    gsl_matrix      *A;             // copy of covariance for LU decomposition
    gsl_matrix      *alpha;         // alpha - PxP partial derivatives
    gsl_vector      *beta;          // beta - Px1 vector
    gsl_vector      *da;            // input vector
    gsl_vector      *x;             // solution vector
    gsl_vector      *a;             // linear OPM parameters
    gsl_vector      *Try;           // trial parameters
    gsl_permutation *p;             // LU permutation
    double          lambda;         // covariance augmentation factor
    double          chi2;           // chi-square
    double          Ochi2;          // previous chi2
    int             signum;         // sign of LU decomposition
    int             its;            // iterations
    int             done = FALSE;

    // gsl allocations
    covar = gsl_matrix_alloc(P, P);
    A = gsl_matrix_alloc(P, P);
    alpha = gsl_matrix_alloc(P, P);
    beta = gsl_vector_alloc(P);
    da = gsl_vector_alloc(P);
    x = gsl_vector_alloc(P);
    a = gsl_vector_alloc(P);
    Try = gsl_vector_alloc(P);
    p = gsl_permutation_alloc(P);

    // convert OPM sensor parameters to vector
    gsl_vector_norm(a);                             // make orientation a unit vector
    FIELD2Vector(&Sensor[S], a);
    gsl_vector_set(a, GAIN_, 1.);                   // initialize unity gain
    chi2 = OPMCoeff(S, Measured, Variance, Current, Dipole, a, alpha, beta);    // solve goodness of fit (chi2) for initial parameters
    Ochi2 = chi2;

    // iterate using perturbed parameters until lambda is tiny
    its = 0;
    lambda = INIT_LAMBDA;
    do {

        // alter linearized fitting matrix by augmenting diagonal elements
        gsl_matrix_memcpy(covar, alpha);            // alpha --> covar
        gsl_matrix_mul_diag(covar, 1. + lambda);    // augment diagonal
        gsl_vector_memcpy(da, beta);                // beta --> da
#ifdef DEBUG
        printf("\ncovar:\n");
        gsl_matrix_fprintf(stdout, covar, "%e");
        printf("da:\n");
        gsl_vector_fprintf(stdout, da, "%e");
        fflush(stdout);
#endif
        // matrix solution for perturbed trial vector
        gsl_matrix_memcpy(A, covar);              // covar --> A
        gsl_linalg_LU_decomp(A, p, &signum);
        gsl_linalg_LU_solve(A, p, da, x);
        gsl_vector_scale(x, 0.5);
#ifdef DEBUG
        printf("x:\n");
        gsl_vector_fprintf(stdout, x, "%e");
        fflush(stdout);
#endif
        // Try = a + x
        gsl_vector_memcpy(Try, a);                  // a --> Try
        gsl_vector_add(Try, x);                     // Try += x

        // did the trial succeed?
        chi2 = OPMCoeff(S, Measured, Variance, Current, Dipole, Try, covar, da);    // test trial parameters
        if (chi2 > Ochi2) {                         // failure -- increase lambda
            lambda *= LAMBDA_UP;
        } else {                                    // success -- accept new solution
            lambda /= LAMBDA_DN;
            Ochi2 = chi2;
            gsl_matrix_memcpy(alpha, covar);        // covar --> alpha
            gsl_vector_memcpy(beta, da);            // da --> beta
            gsl_vector_memcpy(a, Try);              // Try --> a
        }
// printf("pass %d: lambda=%e\n", its, lambda);fflush(stdout);
        its++;
        if (its > MAXITS || lambda > MAX_LAMBDA || lambda < END_LAMBDA)
            done = TRUE;

    } while (!done);

    // return fitted sensor parameters in type FIELD & Scale
    *Scale = gsl_vector_get(a, GAIN_);
    gsl_vector_norm(a);
    Vector2FIELD(a, &Sensor[S]);

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

    // return goodness of fit (reduced chi^2)
    return(chi2);
}


// OPMcoeff()  -- compute coefficients for one OPM sensor
double  OPMCoeff(                   // return chi-square
    int         S,                  // OPM sensor index
    gsl_matrix  *Measured,          // calibrator measurements - S x D (input)
    gsl_matrix  *Variance,          // measured sensor noise variance - S x D (input)
    double      *Current,           // dipole coil currrent (A)
    FIELD       *Dipole,            // array of magnetic dipole positions & orientations (input)
    gsl_vector  *a,                 // OPM parameter vector (input)
    gsl_matrix  *alpha,
    gsl_vector  *beta
)
{
    FIELD       Sensor;             // test OPM sensor
    gsl_vector  *Try;               // test linear parameters            
    double      dyda[P];            // partial derivatives
    double      dy[P];              // difference between measured and trial predicted
    double      measured;           // measured field
    double      predicted;          // predicted field
    double      sigma2;             // variance
    double      gain;               // sensor gain  
    double      dY;                 // difference between measured & predicted
    double      chi2;               // chi-square
    double      wt;                 // dyda / sigma2
    double      tmp;
    double      delta[] = {
        1.0e-04, 1.0e-04, 1.0e-04,  // position deltas
        1.0e-04, 1.0e-04, 1.0e-04,  // orientation deltas
        1.0e-03                     // gain delta
//        1.0e-03, 1.0e-03, 1.0e-03,
//        1.0e-03, 1.0e-03, 1.0e-03,
//        1.0e-03
    };
    int         D;                  // number of dipoles
    int         d;                  // dipole index
    int         i;
    int         j;

    // extract dipole array size
    D = Measured->size2;

    // sanity checks for sensor radial distance & orientation
    tmp = gsl_vector_radius(a);
    if (tmp < MIN_OPM_RAD || tmp > MAX_OPM_RAD)
        return(HUGE);

	// initialize (symmetric) alpha & beta
	Try = gsl_vector_alloc(P);
	gsl_matrix_set_zero(alpha);
	gsl_vector_set_zero(beta);

    // generate 'Sensor' from 'a', with unit orientation
    gsl_vector_norm(a);
    Vector2FIELD(a, &Sensor);

    // sum over all magnetic dipole locations
    for (d=0, chi2=0.; d<D; d++) {

        // compute difference between measured & predicted, & accumulate chi^2 over all dipoles
        measured = gsl_matrix_get(Measured, S, d);
        gain = gsl_vector_get(a, GAIN_);
#ifdef MTOB
        predicted = MtoB(&Dipole[d], &Sensor, gain * Current[d]);
#else
        predicted = Coil2B(&Dipole[d], &Sensor, gain * Current[d]);
#endif
        sigma2 = gsl_matrix_get(Variance, S, d);
        dY = measured - predicted;
#ifdef DEBUG
        printf("dY=%e, fract=%f\n", dY, dY/measured);
        fflush(stdout);
#endif
        // compute dy & dy/da
        for (i=0; i<P; i++) {

            // compute trial predicted field
            for (j=0; j<P; j++)
                if (i == j)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) + delta[j]);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));

            // copy vector parameters to Sensor
            Vector2FIELD(Try, &Sensor);
#ifdef MTOB
            predicted = MtoB(&Dipole[d], &Sensor, gain * Current[d]);
#else
            predicted = Coil2B(&Dipole[d], &Sensor, gain * Current[d]);
#endif
            // compute partial derivative
            dy[i] = measured - predicted;
            dyda[i] = dy[i] / delta[i];
#ifdef DEBUG
            printf("dipole %d: measured=%e, predicted=%e, dy[%d]=%e, dyda[%d]=%e, gain=%e\n", d, measured, predicted, i, dy[i], i, dyda[i], gain);
            fflush(stdout);
#endif
        }

        // accumulate alpha & beta (comments: J: J=dyda, W=1/sigma2)
        for (i=0; i<P; i++) {
            wt = dyda[i] / sigma2;                                  // J'W
            for (j=0; j<P; j++) {
                tmp = gsl_matrix_get(alpha, i, j) + wt * dyda[j];   // J'WJ
                gsl_matrix_set(alpha, i, j, tmp);
            }
            tmp = gsl_vector_get(beta, i) + wt * dy[i];             // J'W dy
            gsl_vector_set(beta, i, tmp);
        }

        // scale alpha & beta
        for (i=0; i<P; i++) {
            gsl_matrix_scale(alpha, 0.01);
            gsl_vector_scale(beta, 0.01);
        }
 
        // accumulate chi2
        chi2 += dY * dY / sigma2;
    }
    free(Try);

    // return chi-square
    return (chi2 / (double)D);
}

