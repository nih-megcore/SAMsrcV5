// MagFitSensor() -- best fit of sensor parameters to multiple
//  magnetic loop data -- solve for position, & orientation,
//  of a sensor, given an array of magnetic loops. This is a
//  Levenberg-Marquardt fit for which the measured and predicted
//  field patterns are normalized.
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <samlib.h>
#include <lfulib.h>
#include <samutil.h>
#include <jig.h>

// #define SHOW_SOL    TRUE
// #define SHOW_PASS   TRUE


// Levenberg-Marquardt nonlinear fit -- input with initial OPM parameters
double  MagFitSensor(       // return chi^2                                 (output)
    FIELD       *Sensor,	// single OPM sensor                            (input & output)
    FIELD       *Dipole,    // Dipole[D] -- magnetic dipole array           (input)               
    gsl_matrix  *Measured,  // OPM measurement -- SxD                       (input)
    gsl_matrix  *Variance,  // OPM measurement variance -- SxD              (input)
    int         S           // sensor index                                 (input)
)
{
    gsl_vector  *a;         // sensor parameters -- Px1
    gsl_vector  *Xm;        // measurements for OPM sensor -- Dx1
    gsl_vector  *Xn;        // measurement variance for OPM sensor -- Dx1
    double      chi2;       // chi^2
    double      Ochi2;      // old chi^2
    double      InitChi2;   // initial chi^2
    double      lambda;     // lambda
    double      diff;       // chi^2 change
    double      change;     // change in chi^2
    double      norm;       // vector magnitude
    double      tmp;        // I dunno what else to call it
    int         pass;       // pass number
    int         test;       // convergence test
    int         D;          // number of magnetic dipoles
    int         d;          // dipole index
    double      MagSensorPass(gsl_vector *a, FIELD *Dipole, gsl_vector *Xm, gsl_vector *Xn, double *lambda);

    // (1) allocate memory
    D = Measured->size2;
    Xm = gsl_vector_alloc(D);
    Xn = gsl_vector_alloc(D);
    a = gsl_vector_alloc(P);

    // (2) extract row vectors for Measured & Variance matrices for single sensor S
    gsl_matrix_get_row(Xm, Measured, S);
    gsl_matrix_get_row(Xn, Variance, S);

    // (3) normalize both measurements & variance
    norm = 1. / gsl_blas_dnrm2(Xm);
    gsl_vector_scale(Xm, norm);
    D = Xm->size;
    for (d=0; d<D; d++) {                   // normalize sigma -- not sigma^2!
        tmp = sqrt(gsl_vector_get(Xn, d)) * norm;
        gsl_vector_set(Xn, d, tmp * tmp);   // save as sigma^2
    }

    // (4) convert sensor parameters to linear array
    FIELD2Vector(Sensor, a);

    // (5) call 'MagSensorPass()' with initial sensor parameters & lambda < 0
    lambda = -1.;       // initialize with lambda < 0
    InitChi2 = MagSensorPass(a, Dipole, Xm, Xn, &lambda);
	chi2 = InitChi2;

	// (6) call 'MagSensorPass()' repeatedly until convergence
    pass = test = 0;
    do {
        Ochi2 = chi2;
        chi2 = MagSensorPass(a, Dipole, Xm, Xn, &lambda);
        diff = Ochi2 - chi2;
        change = fabs(diff) / Ochi2;
#ifdef SHOW_PASS
        printf("pass %d: lambda=%e, chi2=%e, diff=%e, change=%e\n", pass, lambda, chi2, diff, change);
        fflush(stdout);
#endif
        if (diff >= 0.) {        // x^2 getting smaller -- solution converging
            if (change < 1.0e-12)
                test++;         // x^2 making small changes -- getting close
        } else {                // x^2 getting larger -- solution diverging
            if (fabs(diff) > EPSILON)
                test = 0;       // x^2 growing significantly -- no problem yet
		}
        pass++;
    } while (pass <= MAXITS && test <= MAXTEST);

    // (7) if chi^2 decreased, the fit succeeded -- else it failed
    if (chi2 < InitChi2) {      // chi^2 improved -- return new sensor parameters
        Vector2FIELD(a, Sensor);
        gsl_vector_free(a);
        gsl_vector_free(Xm);
        gsl_vector_free(Xn);
        return (chi2);
    } else {					// chi-square didn't decrease -- fit failed
        gsl_vector_free(a);
        gsl_vector_free(Xm);
        gsl_vector_free(Xn);
        return (InitChi2);
    }
}


// MagSensorPass() -- Marquardt-Levenberg non-linear fit
double  MagSensorPass(          // return chi^2                             (output)
    gsl_vector  *a,             // sensor parameters -- Px1                 (input)
    FIELD       *Dipole,        // Dipole[D] -- magnetic dipole array       (input)
    gsl_vector  *Xm,            // normalized measurements -- Dx1           (input)
    gsl_vector  *Xn,            // normalized measurement variances -- Dx1  (input)
    double      *lambda         // damping factor                           (input)
)
{
    gsl_matrix  *covar;         // covariance matrix -- PxP
    gsl_matrix  *alpha;         // alpha matrix -- PxP
    gsl_vector  *beta;          // beta -- Px1
    gsl_vector  *da;            // da -- Px1
    gsl_vector  *x;             // x -- Px1
    gsl_vector  *Try;           // trial parameters -- Px1
    double      chi2;           // chi^2
    double      Ochi2;          // old chi^2
    double      MagSensorCoeff(gsl_vector *a, FIELD *Dipole, gsl_vector *Xm, gsl_vector *Xn, gsl_matrix *alpha, gsl_vector *beta);

    // (1) allocate memory
    covar = gsl_matrix_alloc(P, P);
    alpha = gsl_matrix_alloc(P, P);
    beta = gsl_vector_alloc(P);
    da = gsl_vector_alloc(P);
    x = gsl_vector_alloc(P);
    Try = gsl_vector_alloc(P);

    // (2) initialization
    if (*lambda < 0.)
        *lambda = INIT_LAMBDA;
    chi2 = MagSensorCoeff(a, Dipole, Xm, Xn, alpha, beta);
    Ochi2 = chi2;

    // (3) alter linearized fitting matrix by augmenting diagonal elements
    gsl_matrix_memcpy(covar, alpha);            // alpha --> covar
    gsl_matrix_mul_diag(covar, 1. + *lambda);   // augment diagonal
    gsl_vector_memcpy(da, beta);                // beta --> da
#ifdef SHOW_SOL
        printf("\ncovar:\n");
        gsl_matrix_fprintf(stdout, covar, "%e");
        printf("da:\n");
        gsl_vector_fprintf(stdout, da, "%e");
        fflush(stdout);
#endif
	// (4) matrix solution 
    gsl_linalg_cholesky_decomp1(covar);
    gsl_linalg_cholesky_solve(covar, da, x);
#ifdef SHOW_SOL
        printf("x:\n");
        gsl_vector_fprintf(stdout, x, "%e");
        fflush(stdout);
#endif
	// (5) test if the trial succeeded
	gsl_vector_memcpy(Try, a);              // a --> Try
	gsl_vector_add(Try, x);                 // Try += x
    chi2 = MagSensorCoeff(Try, Dipole, Xm, Xn, alpha, beta);
    if (chi2 < Ochi2) {                     // success -- accept new solution
        *lambda /= LAMBDA_DN;
		Ochi2 = chi2;
        gsl_vector_memcpy(a, Try);          // Try --> a
// printf("\nTrial succeeded: lambda=%e\n", *lambda);fflush(stdout);
    } else {                                // failure -- increase 'lambda' & return
        if (*lambda < MAX_LAMBDA)
            *lambda *= LAMBDA_UP;
// printf("\nTrial failed: lambda=%e\n", *lambda);fflush(stdout);
    }

    // free memory
    gsl_matrix_free(covar);
    gsl_matrix_free(alpha);
    gsl_vector_free(beta);
    gsl_vector_free(da);
    gsl_vector_free(x);
    gsl_vector_free(Try);

    // return reduced chi^2
	return (chi2);
}


// MagSensorCoeff() -- Levenberg-Marquardt coefficients, returns alpha, beta, & chi^2
double  MagSensorCoeff(         // return chi^2                             (output)
    gsl_vector  *a,             // OPM sensor parameters                    (input)
    FIELD       *Dipole,        // Dipole[D] -- magnetic dipole array       (input)
    gsl_vector  *Xm,            // normalized measurements -- Dx1           (input)
    gsl_vector  *Xn,            // normalized measurement variance -- Dx1   (input)
    gsl_matrix  *alpha,         // covariance matrix -- DxD                 (output)
    gsl_vector  *beta           // measured-predicted -- Dx1                (output)
)
{
    gsl_vector  *Xp;            // predicted field for each magnetic dipole (Dx1)
    gsl_vector  *Try;           // trial OPM sensor parameters
    FIELD       Sensor;         // trial OPM sensor
    double      *dyda;          // dyda[P] -- partial derivatives
    double      *dy;            // dy[D] -- observed - predicted
    double      delta = 1.0e-06;    // delta for partial derivatives
    double      wt;             // temporary weight
    double      Pdy;            // positive delta
    double      Ndy;            // negative delta
    double      chi2;           // chi^2
    double      norm;           // normalization factor for predicted value
    double      Bp;             // predicted field
    int         D;              // number of magnetic dipoles
    int         d;              // magnetic dipole index
    int         i;              // i-index
    int         j;              // j-index

    // (1) allocate memory
    D = Xm->size;
    Xp = gsl_vector_alloc(D);
    Try = gsl_vector_alloc(P);
    dyda = new_arrayE(double, P, "dyda[P]");
    dy = new_arrayE(double, D, "dy[D]");

	// (2) initialize (symmetric) alpha & beta
	gsl_matrix_set_zero(alpha);
	gsl_vector_set_zero(beta);

    // (3) make current sensor, compute predicted field, & normalize
    gsl_vector_norm(a);
    Vector2FIELD(a, &Sensor);
    for (d=0; d<D; d++) {
#ifdef MTOB
        Bp = MtoB(&Dipole[d], &Sensor, 1.);
#else
        Bp = Coil2B(&Dipole[d], &Sensor, 1.);
#endif
        gsl_vector_set(Xp, d, Bp);
    }
    norm = 1. / gsl_blas_dnrm2(Xp);
    gsl_vector_scale(Xp, norm);

    // accumulate current chi^2
    for (d=0, chi2=0.; d<D; d++) {
        dy[d] = gsl_vector_get(Xm, d) - gsl_vector_get(Xp, d);
        chi2 += dy[d] * dy[d] / gsl_vector_get(Xn, d);;
    }

    // (4) summation loop over all magnetic dipoles
    for (d=0; d<D; d++) {

        // (4a) compute 'dyda[]'
        for (i=0; i<P; i++) {
            for (j=0; j<P; j++)         // perturb one parameter at a time
                if (j == i)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) + delta);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));
            Vector2FIELD(Try, &Sensor);
#ifdef MTOB
            Pdy = MtoB(&Dipole[d], &Sensor, 1.) * norm;
#else
            Pdy = Coil2B(&Dipole[d], &Sensor, 1.) * norm;
#endif
            for (j=0; j<P; j++)          // perturb one parameter at a time
				if(j == i)
				    gsl_vector_set(Try, j, gsl_vector_get(a, j) - delta);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));
            Vector2FIELD(Try, &Sensor);
#ifdef MTOB
            Ndy = MtoB(&Dipole[d], &Sensor, 1.) * norm;
#else
            Ndy = Coil2B(&Dipole[d], &Sensor, 1.) * norm;
#endif
            dyda[i] = (Pdy - Ndy) / (2. * delta);
        }

        // (4b) accumulate 'alpha' & 'beta'
        for (i=0; i<P; i++) {
            wt = dyda[i] / gsl_vector_get(Xn, d);
            for (j=0; j<P; j++)
                gsl_matrix_plus_equal(alpha, i, j, wt * dyda[j]);   // alpha += wt * dyda
            gsl_vector_plus_equal(beta, i, wt * dy[d]);             // beta += wt * dy
		}
	}

    // scale alpha & beta to lambda range
    gsl_matrix_scale(alpha, 1.0e-08);
    gsl_vector_scale(beta, 1.0e-08);

    // (5) free memory
    gsl_vector_free(Try);
    gsl_vector_free(Xp);
    free(dyda);
    free(dy);

    // (6) return reduced chi^2
	return (chi2 / (double)(D - 1));
}
