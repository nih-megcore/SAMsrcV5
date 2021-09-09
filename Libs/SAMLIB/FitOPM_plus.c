// FitOPM_plus() -- fit OPM a single sensor to measurements from an array of magnetic dipoles
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//  Revisor:    Allison C. Nugent - now fits two angles for orientation and gain

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

#define P           6               //  6-parameter solution
// #define DEBUG       TRUE
// #define SHOW_SOL    TRUE
// #define SHOW_ITS    TRUE
// #define DYDA        TRUE

double  FitOPM_plus(
    int         S,                  // OPM sensor index
    gsl_matrix  *Measured,          // sensor measurements - S x D (input)
    gsl_matrix  *Variance,          // measured noise variance - S x D (input)
    double      *Current,           // rms coil drive current (A)
    OPMFIELD    *Sensor,            // Sensor[S] - OPM sensor array parameters position & orientation (initial guess & outout)
    FIELD       *Dipole,            // Dipole[D] - magnetic dipoles array (input)
    double	param_min[6],	    // parameter minimums
    double      param_max[6]        // parameter maximums
)
{
    gsl_matrix      *covar;         // covar - PxP augmented covariance of partial derivatives
    gsl_matrix      *alpha;         // alpha - PxP partial derivatives
    gsl_vector      *beta;          // beta - Px1 vector
    gsl_vector      *da;            // input vector
    gsl_vector      *x;             // solution vector
    gsl_vector      *a;             // linear OPM parameters
    gsl_vector      *Try;           // trial parameters
    double          lambda;         // covariance augmentation factor
    double          chi2;           // chi-square
    double          Ochi2;          // previous chi2
    double          change;         // change in chi2
    int             its;            // iterations
    int             done = FALSE;
    int 	    i;              // standard counter variable

    gsl_matrix      *U;             // SVD U - 6x6 covariance copy
    gsl_matrix      *V;             // SVD V - 6x6
    gsl_vector      *Sval;          // SVD singular values - 6x1
    gsl_vector      *work;          // SVD working vector - 6x1

    // gsl allocations
    covar = gsl_matrix_alloc(P, P);
    alpha = gsl_matrix_alloc(P, P);
    beta = gsl_vector_alloc(P);
    da = gsl_vector_alloc(P);
    x = gsl_vector_alloc(P);
    a = gsl_vector_alloc(P);
    Try = gsl_vector_alloc(P);

    U = gsl_matrix_alloc(P, P);
    V = gsl_matrix_alloc(P, P);
    Sval = gsl_vector_alloc(P);
    work = gsl_vector_alloc(P);

    // convert OPM sensor parameters to vector
    OPMFIELD2Vector(&Sensor[S], a);

    // solve initial goodness of fit (chi2), alpha, & beta for initial parameters
    chi2 = OPMCoeff_plus(S, Measured, Variance, Current, Dipole, a, alpha, beta);

    printf("Initial chi2 = %e, beginning itertions\n",chi2);

    // iterate using perturbed parameters until lambda is tiny
    its = 0;
    lambda = INIT_LAMBDA;
    Ochi2 = HUGE;
    do {

        // Stopping conditions:
        // Change in chi2 < epsilon

        // evaluate change in chi^2
        change = fabs(chi2 - Ochi2) / Ochi2;
        if (change < EPSILON) {
            done = TRUE;
            printf("change in chi2 < EPSILON, gonna stop after one last calc of new params\n");
        }

        // Ochi2
        // update Ochi2 & check if small enough for linear fit
        Ochi2 = chi2;
        if (Ochi2 < 10.) {
            lambda = 0.;
            done = TRUE;
            printf("Ochi2 < 10, gonna stop after one last calc of new params\n");
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

        // Cholesky matrix solution for perturbed trial vector
        // gsl_linalg_cholesky_decomp1(covar);
        // gsl_linalg_cholesky_solve(covar, da, x);

        // SVD matrix solution for perturbed trial vector
        gsl_matrix_memcpy(U, covar);      	   // covar --> u
        gsl_linalg_SV_decomp(U, V, Sval, work);
        gsl_linalg_SV_solve(U, V, Sval, da, x);

#ifdef SHOW_SOL
        printf("x:\n");
        gsl_vector_fprintf(stdout, x, "%e");
        fflush(stdout);
#endif
        // Try = a + x
        gsl_vector_memcpy(Try, a);                  // a --> Try
        gsl_vector_add(Try, x);                     // Try += x

        // make sure that the new trial parameters don't go outside the maximum and minimum values

        for(i=0; i<P; i++) {
          if (gsl_vector_get(Try,i) > param_max[i]) {
             printf("trial param %d a=%e is larger than param_max=%e \n",i,gsl_vector_get(Try,i),param_max[i]);
             gsl_vector_set(Try, i, param_max[i]);
          } else if (gsl_vector_get(Try,i) < param_min[i]) {
             printf("trial param %d a=%e is smaller than param_min=%e \n",i,gsl_vector_get(Try,i),param_min[i]);
             gsl_vector_set(Try, i, param_min[i]);
          }
        } 

        printf("Calculating chi2 for new params, iteration %d\n",its);
        chi2 = OPMCoeff_plus(S, Measured, Variance, Current, Dipole, Try, covar, da);    // test trial parameters

        // test trial parameters -- did the trial succeed?

        if (chi2 <= Ochi2) {                        // try succeded -- accept new solution
            gsl_matrix_memcpy(alpha, covar);        // covar --> alpha
            gsl_vector_memcpy(beta, da);            // da --> beta
            gsl_vector_memcpy(a, Try);              // Try --> a
            lambda /= LAMBDA_DN;
            printf("Success, chi2=%e, Ochi2=%e, decreasing lambda to %e\n", chi2, Ochi2, lambda);

        } else {                                    // try failed -- increase lambda
            if (lambda < MAX_LAMBDA)
                lambda *= LAMBDA_UP;
            printf("Failure, chi2=%e, Ochi2=%e, increasing lambda to %e\n", chi2, Ochi2, lambda);
        }

#ifdef SHOW_ITS
        printf("\npass %d: chi2=%e, lambda=%e\n", its, chi2, lambda);
        fflush(stdout);
#endif

        its++;
        if (its > MAXITS || lambda > MAX_LAMBDA) {
            printf("Stopping, its=%d greater than MAXITS or lambda=%e > MAX_LAMBDA\n",its,lambda);
            done = TRUE;
        }

    } while (!done);

    // return fitted sensor parameters in type FIELD & Scale
    //*Scale = gsl_vector_norm2(a);

    Vector2OPMFIELD(a, &Sensor[S]);
    printf("final a[0]: %f, a[1]: %f,a[2]: %f,a[3]: %f,a[4]: %f,a[5]: %f\n", gsl_vector_get(a,0),gsl_vector_get(a,1),gsl_vector_get(a,2),gsl_vector_get(a,3),gsl_vector_get(a,4),gsl_vector_get(a,5));


    // free gsl allocations
    gsl_matrix_free(covar);
    gsl_matrix_free(alpha);
    gsl_vector_free(beta);
    gsl_vector_free(da);
    gsl_vector_free(x);
    gsl_vector_free(a);
    gsl_vector_free(Try);
    
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(Sval);
    gsl_vector_free(work);

    // return goodness of fit (reduced chi^2)
    return(chi2);
}


// OPMcoeff_plus()  -- compute coefficients for one OPM sensor
double  OPMCoeff_plus(              // return chi-square
    int         s,                  // OPM sensor index
    gsl_matrix  *Measured,          // calibrator measurements - S x D (input)
    gsl_matrix  *Variance,          // measured sensor noise variance - S x D (input)
    double      *Current,           // dipole coil currrent (A)
    FIELD       *Dipole,            // array of magnetic dipole positions & orientations (input)
    gsl_vector  *a,                 // OPM parameter vector (input)
    gsl_matrix  *alpha,
    gsl_vector  *beta
)
{
    OPMFIELD       Sensor;             // test OPM sensor
    gsl_vector  *Try;               // test linear parameters            
    double      *dyda;              // dyda[P] -- partial derivatives
    double      measured;           // measured field
    double      predicted;          // predicted field
    double      Bp;                 // positive delta field
    double      Bn;                 // negative delta field
    double      dy;                 // difference between measured & predicted
    double      sigma2;             // variance
    double      chi2;               // chi-square
    double      wt;                 // dyda / sigma
    double      tmp;
    //double      delta = 1.0e-04;
    double      delta[] = {         // perturbations for position and orientation
        1.0e-8, 1.0e-8, 1.0e-8,  // position delta
        1.0e-8, 1.0e-8, 1.0e-4  // orientation delta + gain delta
     };
    int         D;                  // number of dipoles
    int         d;                  // dipole index
    int         i;
    int         j;

    // extract dipole array
    D = Measured->size2;

    printf("\nOPMMCoeff: Calculating fits for sensor S=%d\n",s);

    // allocate trial array & partials
	Try = gsl_vector_alloc(P);
	dyda = new_arrayE(double, P, "dyda[P]");

    // initialize (symmetric) alpha & beta
	gsl_matrix_set_zero(alpha);
	gsl_vector_set_zero(beta);

    // generate 'Sensor' from 'a'
    Vector2OPMFIELD(a, &Sensor);

    // sum over all magnetic dipole locations
    for (d=0, chi2=0.; d<D; d++) {

        // compute difference between measured & predicted, & accumulate chi^2 over all dipoles

        measured = gsl_matrix_get(Measured, s, d);
        sigma2 = gsl_matrix_get(Variance, s, d);

#ifdef MTOB
        predicted = MtoB_plus(&Dipole[d], &Sensor, Current[d]);
#else
        predicted = Coil2B_plus(&Dipole[d], &Sensor, Current[d]);
#endif

        dy = measured - predicted;
        chi2 += dy * dy / sigma2;

//#ifdef DYDA
        printf("\ndipole %d: measured=%e, predicted=%e, dy=%e\n", d, measured, predicted, dy);
        fflush(stdout);
//#endif

        // compute dy/da
        for (i=0; i<P; i++) {

            // compute (+) delta predicted field
            for (j=0; j<P; j++)
                if (i == j)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) + delta[j]);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));

            // copy 'Try' vector parameters to Sensor
            Vector2OPMFIELD(Try, &Sensor);
#ifdef MTOB
            Bp = MtoB_plus(&Dipole[d], &Sensor, Current[d]);
#else
            Bp = Coil2B_plus(&Dipole[d], &Sensor, Current[d]);
#endif

            // compute (-) delta predicted field
            for (j=0; j<P; j++)
                if (i == j)
                    gsl_vector_set(Try, j, gsl_vector_get(a, j) - delta[j]);
                else
                    gsl_vector_set(Try, j, gsl_vector_get(a, j));

            // copy vector parameters to Sensor
            Vector2OPMFIELD(Try, &Sensor);
#ifdef MTOB
            Bn = MtoB_plus(&Dipole[d], &Sensor, Current[d]);
#else
            Bn = Coil2B_plus(&Dipole[d], &Sensor, Current[d]);
#endif

            // compute partial derivative
            dyda[i] = (Bp - Bn) / (2. * delta[i]);
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
    }

    // scale alpha & beta to lambda range
    gsl_matrix_scale(alpha, 1.0e-12);
    gsl_vector_scale(beta, 1.0e-12);

    gsl_vector_free(Try);
    free(dyda);

    // return reduced chi-square
    return (chi2 / (double)(D - P + 1));
}
