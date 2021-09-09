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

#define P   6       // 6-parameter solution
// #define DEBUG       TRUE
// #define SHOW_SOL    TRUE
// #define SHOW_ITS    TRUE
// #define DYDA        TRUE

double  FitDip(
    int         D,                  // dipole coil index
    gsl_matrix  *Measured,          // calibrator measurements - S x D (input)
    gsl_matrix  *Variance,          // measured sensor noise variance - S x D (input)
    double      *Current,           // Current[N_DIP] - rms coil current (A) or moment (A-m^2) (input)
    FIELD       *Sensor,            // Sensor[N_OPM] - array of OPM sensors (input)
    FIELD       *Dipole             // Dipole[N_DIP] - array of magnetic dipoles (initial guess & output)
)
{
    gsl_matrix      *covar;         // covar - PxP augmented covariance of partial derivatives
    gsl_matrix      *alpha;         // alpha - PxP partial derivatives
    gsl_vector      *beta;          // beta - Px1 vector
    gsl_vector      *da;            // input vector - Px1
    gsl_vector      *x;             // solution vector - Px1
    gsl_vector      *a;             // linear OPM parameters - Px1
    gsl_vector      *Try;           // trial parameters - Px1
    double          lambda;         // covariance augmentation factor
    double          chi2;           // chi-square
    double          Ochi2;          // previous chi2
    double          change;         // change in chi2
    int             its;            // iterations
    int             done = FALSE;

    // gsl allocations
    covar = gsl_matrix_alloc(P, P);
    alpha = gsl_matrix_alloc(P, P);
    beta = gsl_vector_alloc(P);
    da = gsl_vector_alloc(P);
    x = gsl_vector_alloc(P);
    a = gsl_vector_alloc(P);
    Try = gsl_vector_alloc(P);

    // convert parameters of single dipole to vector
    FIELD2Vector(&Dipole[D], a);
    gsl_vector_norm(a);                             // initialize with unit vector orientation

    // solve goodness of fit (chi2), alpha, & beta for initial parameters
    chi2 = DipCoeff(D, Measured, Variance, Current, Sensor, a, alpha, beta);

    // iterate using perturbed parameters until lambda is tiny
    its = 0;
    lambda = INIT_LAMBDA;
    Ochi2 = HUGE;
    do {

        // evaluate change in chi^2
        change = fabs(chi2 - Ochi2) / Ochi2;
        if (change < EPSILON)
            done = TRUE;

        // update Ochi2 & check if small enough for linear fit
        Ochi2 = chi2;
        if (Ochi2 < 10.) {
            lambda = 0.;
            done = TRUE;
        }

        // alter linearized fitting matrix by augmenting diagonal elements
        gsl_matrix_memcpy(covar, alpha);            // alpha --> covar
        gsl_matrix_mul_diag(covar, 1. + lambda);    // augment diagonal
        gsl_vector_memcpy(da, beta);                // beta --> da
#ifdef SHOW_SOL
        printf("\ncovar:\n");
        gsl_matrix_fprintf(stdout, covar, "%e");
        printf("da:\n");
        gsl_vector_fprintf(stdout, da, "%e");
        fflush(stdout);
#endif
        // matrix solution for perturbed trial vector
        gsl_linalg_cholesky_decomp1(covar);
        gsl_linalg_cholesky_solve(covar, da, x);
#ifdef SHOW_SOL
        printf("x:\n");
        gsl_vector_fprintf(stdout, x, "%e");
        fflush(stdout);
#endif
        // Try = a + x
        gsl_vector_memcpy(Try, a);                  // a --> Try
        gsl_vector_add(Try, x);                     // Try += x

        // did the trial succeed?
        chi2 = DipCoeff(D, Measured, Variance, Current, Sensor, Try, covar, da);    // test trial parameters
        if (chi2 <= Ochi2) {                        // try succeeded -- accept new solution
            gsl_matrix_memcpy(alpha, covar);        // covar --> alpha
            gsl_vector_memcpy(beta, da);            // da --> beta
            gsl_vector_memcpy(a, Try);              // Try --> a
            if (chi2 > 10.)
                lambda /= LAMBDA_DN;
// printf("\nsuccess: lambda=%e\n", lambda);fflush(stdout);
        } else {                                    // try failed -- increase lambda
            if (lambda < MAX_LAMBDA)
                lambda *= LAMBDA_UP;
// printf("\nfailure: lambda=%e\n", lambda);fflush(stdout);
        }
#ifdef SHOW_ITS
        printf("\npass %d: chi2=%e, lambda=%e\n", its, chi2, lambda);
        fflush(stdout);
#endif
        its++;
        if (its > MAXITS || lambda > MAX_LAMBDA)
            done = TRUE;;

    } while(!done);

    // return dipole parameters in type FIELD
    gsl_vector_norm(a);
    Vector2FIELD(a, &Dipole[D]);

    // free gsl allocations
    gsl_matrix_free(covar);
    gsl_matrix_free(alpha);
    gsl_vector_free(beta);
    gsl_vector_free(da);
    gsl_vector_free(x);
    gsl_vector_free(a);
    gsl_vector_free(Try);

    // return goodness of fit
    return(chi2);
}


// compute coefficients for one dipole
double  DipCoeff(                   // return chi-square
    int         d,                  // dipole coil index
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
    gsl_vector  *Try;               // test parameters - Px1
    double      *dyda;              // dyda[P] -- partial derivatives
    double      dy;                 // difference between measured & predicted
    double      measured;           // measured field
    double      predicted;          // predicted field
    double      Bp;                 // positive delta field
    double      Bn;                 // negative delta field
    double      sigma2;             // variance
    double      chi2;               // chi-square
    double      wt;                 // dyda / sigma2
    double      delta = 1.0e-04;    // parameter displacement
    double      tmp;
    int         S;                  // number of OPM sensors
    int         i;
    int         j;
    int         s;

    // extract sensor array size
    S = Measured->size1;

    // allocate trial array & partials
    Try = gsl_vector_alloc(P);
    dyda = new_arrayE(double, P, "dyda[P]");

	// initialize (symmetric) alpha & beta
	gsl_matrix_set_zero(alpha);
	gsl_vector_set_zero(beta);

    // generate Dipole from 'a'
    Vector2FIELD(a, &Dipole);

    // sum over all magnetic sensor locations
    for (s=0, chi2=0.; s<S; s++) {

        // compute difference between predicted & measured field & accumulate chi2
        measured = gsl_matrix_get(Measured, s, d);
        sigma2 = gsl_matrix_get(Variance, s, d);
#ifdef MTOB
        predicted = MtoB(&Dipole, &Sensor[s], Current[d]);
#else
        predicted = Coil2B(&Dipole, &Sensor[s], Current[d]);
#endif
        dy = measured - predicted;
#ifdef DYDA
        printf("\nsensor %d: measured=%e, predicted=%e, dy=%e\n", s, measured, predicted, dy);
        fflush(stdout);
#endif
        // compute dy/da
        for (i=0; i<P; i++) {

            // compute (+) delta predicted field
            for (j=0; j<P; j++)
                if (i == j)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) + delta);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));

            // copy 'Try' vector parameters to Dipole
            Vector2FIELD(Try, &Dipole);
#ifdef MTOB
            Bp = MtoB(&Dipole, &Sensor[s], Current[d]);
#else
            Bp = Coil2B(&Dipole, &Sensor[s], Current[d]);
#endif

            // compute (-) delta predicted field
            for (j=0; j<P; j++)
                if (i == j)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) - delta);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));

            // copy 'Try' vector parameters to Dipole
            Vector2FIELD(Try, &Dipole);
#ifdef MTOB
            Bn = MtoB(&Dipole, &Sensor[s], Current[d]);
#else
            Bn = Coil2B(&Dipole, &Sensor[s], Current[d]);
#endif

            // compute partial derivative
            dyda[i] = (Bp - Bn) / (2. * delta);
#ifdef DYDA
            printf("\tdyda[%d]=%e\n", i, dyda[i]);
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
            tmp = gsl_vector_get(beta, i) + wt * dy;                // J'W dy
            gsl_vector_set(beta, i, tmp);
        }

        // accumulate chi2
        chi2 += dy * dy / sigma2;
    }

    // scale alpha & beta to lambda range
    gsl_matrix_scale(alpha, 1.0e-12);
    gsl_vector_scale(beta, 1.0e-12);

    gsl_vector_free(Try);
    free(dyda);

    // return reduced chi-square
    return (chi2 / (double)(S - P + 1));
}
