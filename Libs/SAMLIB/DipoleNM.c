// DipoleNM() -- Nelder-Mead downhill simplex minimization in
//  multiple dimensions. Fit dipole to multiple magnetic
//  sensor measurements.
//
//	Author:	Stephen E. Robinson
//			MEG Core Group
//			NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <samlib.h>
#include <samutil.h>

#define	Alpha           1.
#define Beta            .5
#define Gamma           2.
#define InitTolerance   1.0e-06
#define IterationLimit  100000
#define NumDimensions   6
#define NumVertices     7

double  DipoleNM(           // return chi^2 of low point
    FIELD       *Sensor,    // Sensor[S] -- OPM sensor array        (input & output)
    FIELD       *Dipole,    // Dipole[D] -- magnetic dipole array   (input)
    gsl_matrix  *Xm,        // measured field -- SxD                (input)
    gsl_matrix  *Xv,        // measurement variance -- SxD          (input)
    double      *Current,   // Current[D] -- dipole current array   (input)
    int         Dindex      // sensor index                         (input)
)
{
    double      *measured;  // measured[S] -- measurements for one magnetic dipole
    double      *sigma2;    // sigma2[S] -- measurement variance for one magnetic dipole
    double      current;    // magnetic dipole current for single dipole
    double      Px[NumVertices][NumDimensions];
    double      Yx[NumVertices];
    double      PBar[NumDimensions];
    double      PRef[NumDimensions];
    double      PRefRef[NumDimensions];
    double      YPRef;
    double      YPRefRef;
    double      Tolerance;
    double      ToleranceLimit = InitTolerance;
    double      max;
    double      min;
    int         High;
    int         NextHigh;
    int         Low;
    int         Iterations;
    int         S;
    int         s;
    int         i;
    int         j;
    double      DipoleFunction(double *Px, FIELD *Sensor, double current, double *measured, double *sigma2, int S);

    // extract vectors for measurements & variance
    S = Xm->size1;
    measured = new_arrayE(double, S, "measured[S]");
    sigma2 = new_arrayE(double, S, "sigma2[S]");
    current = Current[Dindex];
    for (s=0; s<S; s++) {
        measured[s] = gsl_matrix_get(Xm, s, Dindex);
        sigma2[s] = gsl_matrix_get(Xv, s, Dindex);
    }

    // initialize the simplex
    FIELDtoA6(&Dipole[Dindex], Px[0]);
    Yx[0] = DipoleFunction(Px[0], Sensor, current, measured, sigma2, S);
    for (i=1; i<NumVertices; i++) {
        for (j=0; j<NumDimensions; j++)
            if (j == i - 1)
                Px[i][j] = 1.25 * Px[0][j];
            else
                Px[i][j] = Px[0][j];
        Yx[i] = DipoleFunction(Px[i], Sensor, current, measured, sigma2, S);
	}

    // iterative search minimizing 'SensorFunction()'
    for (Iterations=0;;) {
        for (i=0, max=(-HUGE), min=HUGE; i<NumVertices; i++) {
            if (Yx[i] > max) {
                max = Yx[i];
                High = i;
            }
            if (Yx[i] < min) {
                min = Yx[i];
                Low = i;
            }
        }
        for (i=0, max=(-HUGE); i<NumVertices; i++)
            if (i != High && Yx[i] > max) {
                max = Yx[i];
                NextHigh = i;
            }
        if (High == NextHigh || High == Low || NextHigh == Low)
            Cleanup("degenerate simplex!");
        Tolerance = 2. * (Yx[High] - Yx[Low]) / (Yx[High] + Yx[Low]);
        if (Tolerance < ToleranceLimit) {       // normal exit upon solution
            A6toFIELD(Px[Low], &Dipole[Dindex]);
            return (Yx[Low]);
        }
        if (Iterations++ == IterationLimit) {   // if tolerance is too high
            ToleranceLimit *= 10.;              //  decrease tolerance by 10
            Iterations = 0;
        }
//            Cleanup("iteration limit reached!");
        for (i=0; i<NumDimensions; i++)
            PBar[i] = 0.;
        for (i=0; i<NumVertices; i++)
            if(i != High)
                for (j=0; j<NumDimensions; j++)
                    PBar[j] += Px[i][j];
        for (i=0; i<NumDimensions; i++) {
            PBar[i] /= (double)NumDimensions;
            PRef[i] = (1. + Alpha) * PBar[i] - Alpha * Px[High][i];
        }
        YPRef = DipoleFunction(PRef, Sensor, current, measured, sigma2, S);
        if (YPRef <= Yx[Low]) {
            for (i=0; i<NumDimensions; i++)
                PRefRef[i] = Gamma * PRef[i] + (1. - Gamma) * PBar[i];
            YPRefRef = DipoleFunction(PRefRef, Sensor, current, measured, sigma2, S);
            if (YPRefRef < Yx[Low]) {
                for (i=0; i<NumDimensions; i++)
                    Px[High][i] = PRefRef[i];
                Yx[High] = YPRefRef;
            } else {
                for (i=0; i<NumDimensions; i++)
                    Px[High][i] = PRef[i];
                Yx[High] = YPRef;
            }
        } else if (YPRef >= Yx[NextHigh]) {
            if (YPRef < Yx[High]) {
                for (i=0; i<NumDimensions; i++)
                    Px[High][i] = PRef[i];
                Yx[High] = YPRef;
            }
            for (i=0; i<NumDimensions; i++)
                PRefRef[i] = Beta * Px[High][i] + (1. - Beta) * PBar[i];
            YPRefRef = DipoleFunction(PRefRef, Sensor, current, measured, sigma2, S);
            if (YPRefRef < Yx[High]) {
                for (i=0; i<NumVertices; i++)
                    if (i != Low) {
                        for (j=0; j<NumDimensions; j++) {
                            PRef[j] = .5 * (Px[i][j] + Px[Low][j]);
                            Px[i][j] = PRef[j];
                        }
                        Yx[i] = DipoleFunction(PRef, Sensor, current, measured, sigma2, S);
                    }
            }
        } else {
            for (i=0; i<NumDimensions; i++)
                Px[High][i] = PRef[i];
            Yx[High] = YPRef;
        }
    }
}


// return reduced chi^2 for measured - predicted, for all dipoles
double  DipoleFunction(
    double      *Parms,         // sensor parameters for simplex vertex     (input)
    FIELD       *Sensor,        // Sensor[S] -- magnetic dipole array       (input)
    double      current,        // magnetic dipole current                  (input)
    double      *measured,      // measured[S] -- dipole array field        (input)
    double      *sigma2,        // sigma2[S] -- measurement variance        (input)
    int         S               // number of sensors                        (input)
) {
	FIELD       Dipole;         // test dipole
	double      predicted;      // predicted field for dipole
	double  	dy;             // measured - predicted
	double      chi2;           // goodness of fit
	int         s;              // sensor index

	A6toFIELD(Parms, &Dipole);
    for (s=0, chi2=0.; s<S; s++) {
	    predicted = Coil2B(&Dipole, &Sensor[s], current);
	    dy = measured[s] - predicted;
        chi2 += dy * dy / sigma2[s];
    }
//    printf("chi2=%e\n", chi2/(double)(S-1));fflush(stdout);
    return (chi2 /  (double)(S - 1));
}

