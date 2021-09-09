// MagFitDipole() -- best fit of magnetic dipole parameters to
//  multiple OPM sensor data. This fit is independent of the
//  dipole magnitude (gain) because it fits the normalized
//  measurements to the normalized predicted fields.
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


double  MagFitDipole(       // return chi^2                                 (output)
    FIELD       *Dipole,	// single magnetic dipole                       (input & output)
    FIELD       *Sensor,    // Sensor[S] -- array of OPM sensors            (input)
    gsl_matrix  *Measured,  // measurements for OPM sensor -- SxD           (input)
    gsl_matrix  *Variance,  // measurement variance for OPM sensor -- SxD   (input)
    int         D           // magnetic dipole index                        (input)
)
{
    gsl_vector  *a;         // sensor parameters -- Px1
    gsl_vector  *Xm;        // measurements for magnetic dipole -- Sx1
    gsl_vector  *Xn;        // measurement variance for magnetic dipole -- Sx1
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
    int         S;          // number of magnetic dipoles
    int         s;          // dipole index
    double      MagDipolePass(gsl_vector *a, FIELD *Dipole, gsl_vector *Xm, gsl_vector *Xn, double *lambda);

    // (1) allocate memory
    S = Measured->size1;
    Xm = gsl_vector_alloc(S);
    Xn = gsl_vector_alloc(S);
    a = gsl_vector_alloc(P);

    // (2) extract column vectors for Measured & Variance matrices for single dipole D
    gsl_matrix_get_col(Xm, Measured, D);
    gsl_matrix_get_col(Xn, Variance, D);

    // (3) normalize both measurements & variance
    norm = 1. / gsl_blas_dnrm2(Xm);
    gsl_vector_scale(Xm, norm);
    S = Xm->size;
    for (s=0; s<S; s++) {                   // normalize sigma -- not sigma^2!
        tmp = sqrt(gsl_vector_get(Xn, s)) * norm;
        gsl_vector_set(Xn, s, tmp * tmp);   // save as sigma^2
    }

    // (4) convert magnetic dipole parameters to linear array
    FIELD2Vector(Dipole, a);

    // (5) call 'MagDipolePass()' with initial sensor parameters & lambda < 0
    lambda = -1.;       // initialize with lambda < 0
    InitChi2 = MagDipolePass(a, Sensor, Xm, Xn, &lambda);
	chi2 = InitChi2;

	// (6) call 'MagDipolePass()' repeatedly until convergence
    pass = test = 0;
    do {
        Ochi2 = chi2;
        chi2 = MagDipolePass(a, Sensor, Xm, Xn, &lambda);
        diff = Ochi2 - chi2;
        change = fabs(Ochi2 - chi2) / Ochi2;
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
    } while(pass < MAXITS && test <= MAXTEST);

    // (7) free memory


    // (7) if chi^2 decreased, the fit succeeded -- else it failed
    if (chi2 < InitChi2) {      // chi^2 improved -- return new sensor parameters
        Vector2FIELD(a, Dipole);
        gsl_vector_free(a);
        gsl_vector_free(Xm);
        gsl_vector_free(Xn);
        return (chi2);
    } else						// chi-square didn't decrease -- fit failed
        gsl_vector_free(a);
        gsl_vector_free(Xm);
        gsl_vector_free(Xn);
        return (InitChi2);
}


// MagDipolePass() -- Marquardt-Levenberg non-linear fit
double  MagDipolePass(          // return chi^2                             (output)
    gsl_vector  *a,             // sensor parameters -- Px1                 (input)
    FIELD       *Sensor,        // Sensor[S] -- OPM sensor array            (input)
    gsl_vector  *Xm,            // normalized measurements -- Sx1           (input)
    gsl_vector  *Xn,            // normalized measurement variances -- Sx1  (input)
    double      *lambda         // damping factor                           (input)
)
{
    gsl_matrix  *covar;         // PxP -- covariance matrix
    gsl_matrix  *alpha;         // PxP -- alpha matrix
    gsl_vector  *beta;          // Px1 -- beta
    gsl_vector  *da;            // Px1 -- da
    gsl_vector  *x;             // Px1 -- x
    gsl_vector  *Try;           // Px1 -- trial parameters
    double      chi2;           // chi^2
    double      Ochi2;          // old chi^2
    double      MagDipoleCoeff(gsl_vector *a, FIELD *Sensor, gsl_vector *Xm, gsl_vector *Xn, gsl_matrix *alpha, gsl_vector *beta);

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
    chi2 = MagDipoleCoeff(a, Sensor, Xm, Xn, alpha, beta);
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
    chi2 = MagDipoleCoeff(Try, Sensor, Xm, Xn, alpha, beta);
    if (chi2 < Ochi2) {                     // success -- accept new solution
        *lambda /= LAMBDA_DN;
		Ochi2 = chi2;
        gsl_vector_memcpy(a, Try);          // Try --> a
// printf("Trial succeeded: lambda=%e\n", *lambda);fflush(stdout);
    } else {                                // failure -- increase 'lambda' & return
        if (*lambda < MAX_LAMBDA)
            *lambda *= LAMBDA_UP;
// printf("Trial failed: lambda=%e\n", *lambda);fflush(stdout);
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


// MagDipoleCoeff() -- Marquardt-Levenberg non-linear least-squares fit
double  MagDipoleCoeff(                 // return chi^2                             (output)
    gsl_vector  *a,             // OPM sensor parameters                    (input)
    FIELD       *Sensor,        // Sensor[S] -- OPM sensor array            (input)
    gsl_vector  *Xm,            // normalized measurements -- Sx1           (input)
    gsl_vector  *Xn,            // normalized measurement variance -- Sx1   (input)
    gsl_matrix  *alpha,         // covariance matrix -- SxS                 (output)
    gsl_vector  *beta           // measured-predicted -- Sx1                (output)
)
{
    gsl_vector  *Xp;            // predicted field for each magnetic dipole (Dx1)
    gsl_vector  *Try;           // trial magnetic dipole parameters
    FIELD       Dipole;         // trial magnetic dipole
    double      *dyda;          // dyda[P] -- partial derivatives
    double      *dy;            // dy[S] -- obseved-predicted delta
    double      delta = 1.0e-06;    // delta for partial derivatives
    double      wt;             // temporary weight
    double      Pdy;            // positive delta
    double      Ndy;            // negative delta
    double      chi2;           // chi^2
    double      norm;           // normalization factor for predicted value
    double      Bp;             // predicted field
    int         S;              // number of OPM sensors
    int         s;              // OPM sensor index
    int         i;              // i-index
    int         j;              // j-index

    // (1) allocate memory
    S = Xm->size;
    Xp = gsl_vector_alloc(S);
    Try = gsl_vector_alloc(P);
    dyda = new_arrayE(double, P, "dyda[P]");
    dy = new_arrayE(double, S, "dy[S]");

	// (2) initialize (symmetric) alpha & beta
	gsl_matrix_set_zero(alpha);
	gsl_vector_set_zero(beta);

    // (3) make current dipole, compute predicted field, & normalize
    gsl_vector_norm(a);
    Vector2FIELD(a, &Dipole);
    for (s=0; s<S; s++) {
#ifdef MTOB
        Bp = MtoB(&Dipole, &Sensor[s], 1.);
#else
        Bp = Coil2B(&Dipole, &Sensor[s], 1.);
#endif
        gsl_vector_set(Xp, s, Bp);
    }
    norm = 1. / gsl_blas_dnrm2(Xp);
    gsl_vector_scale(Xp, norm);

    // accumulate current chi^2
    for (s=0, chi2=0.; s<S; s++) {
        dy[s] = gsl_vector_get(Xm, s) - gsl_vector_get(Xp, s);
        chi2 += dy[s] * dy[s] / gsl_vector_get(Xn, s);;
    }

    // (4) summation loop over all OPM sensors
    for (s=0; s<S; s++) {

        // (4a) compute 'dyda[]'
        for (i=0; i<P; i++) {
            for (j=0; j<P; j++)         // perturb one parameter at a time
                if (j == i)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) + delta);
				else
					gsl_vector_set(Try, j, gsl_vector_get(a, j));
            Vector2FIELD(Try, &Dipole);
#ifdef MTOB
            Pdy = MtoB(&Dipole, &Sensor[s], 1.) * norm;
#else
            Pdy = Coil2B(&Dipole, &Sensor[s], 1.) * norm;
#endif
            for (j=0; j<P; j++)          // perturb one parameter at a time
				if(j == i)
				    gsl_vector_set(Try, j, gsl_vector_get(a, j) - delta);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));
            Vector2FIELD(Try, &Dipole);
#ifdef MTOB
            Ndy = MtoB(&Dipole, &Sensor[s], 1.) * norm;
#else
            Ndy = Coil2B(&Dipole, &Sensor[s], 1.) * norm;
#endif
            dyda[i] = (Pdy - Ndy) / (2. * delta);
        }

        // (4b) accumulate 'alpha' & 'beta'
        for (i=0; i<P; i++) {
            wt = dyda[i] / gsl_vector_get(Xn, s);
            for (j=0; j<P; j++)
                gsl_matrix_plus_equal(alpha, i, j, wt * dyda[j]);   // alpha += wt * dyda
            gsl_vector_plus_equal(beta, i, wt * dy[s]);             // beta +- wt * dy
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
	return (chi2 / (double)(S - 1));
}
