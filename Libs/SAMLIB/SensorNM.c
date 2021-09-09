// SensorNM() -- Nelder-Mead downhill simplex minimization in
//  multiple dimensions. Fit sensor to multiple magnetic
//  dipole measurements.
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

double  SensorNM(           // return chi^2 of low point
    FIELD       *Sensor,    // Sensor[S] -- OPM sensor array        (input & output)
    FIELD       *Dipole,    // Dipole[D] -- magnetic dipole array   (input)
    gsl_matrix  *Xm,        // measured field -- SxD                (input)
    gsl_matrix  *Xv,        // measurement variance -- SxD          (input)
    double      *current,   // current[D] -- dipole current array   (input)
    int         Sindex      // sensor index                         (input)
)
{
    double      *measured;  // measured[D] -- measurements for one OPM sensor
    double      *sigma2;    // sigma2[D] -- measurement variance for one OPM sensor
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
    int         D;
    int         d;
    int         i;
    int         j;
    double      SensorFunction(double *Px, FIELD *Dipole, double *current, double *measured, double *sigma2, int D);

    // extract vectors for measurements & variance
    D = Xm->size2;
    measured = new_arrayE(double, D, "measured[D]");
    sigma2 = new_arrayE(double, D, "sigma2[D]");
    for (d=0; d<D; d++) {
        measured[d] = gsl_matrix_get(Xm, Sindex, d);
        sigma2[d] = gsl_matrix_get(Xv, Sindex, d);
    }

    // initialize the simplex
    FIELDtoA6(&Sensor[Sindex], Px[0]);
    Yx[0] = SensorFunction(Px[0], Dipole, current, measured, sigma2, D);
    for (i=1; i<NumVertices; i++) {
        for (j=0; j<NumDimensions; j++)
            if (j == i - 1)
                Px[i][j] = 1.25 * Px[0][j];
            else
                Px[i][j] = Px[0][j];
        Yx[i] = SensorFunction(Px[i], Dipole, current, measured, sigma2, D);
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
            A6toFIELD(Px[Low], &Sensor[Sindex]);
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
        YPRef = SensorFunction(PRef, Dipole, current, measured, sigma2, D);
        if (YPRef <= Yx[Low]) {
            for (i=0; i<NumDimensions; i++)
                PRefRef[i] = Gamma * PRef[i] + (1. - Gamma) * PBar[i];
            YPRefRef = SensorFunction(PRefRef, Dipole, current, measured, sigma2, D);
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
            YPRefRef = SensorFunction(PRefRef, Dipole, current, measured, sigma2, D);
            if (YPRefRef < Yx[High]) {
                for (i=0; i<NumVertices; i++)
                    if (i != Low) {
                        for (j=0; j<NumDimensions; j++) {
                            PRef[j] = .5 * (Px[i][j] + Px[Low][j]);
                            Px[i][j] = PRef[j];
                        }
                        Yx[i] = SensorFunction(PRef, Dipole, current, measured, sigma2, D);
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
double  SensorFunction(
    double      *Parms,         // sensor parameters for simplex vertex     (input)
    FIELD       *Dipole,        // Dipole[D] -- magnetic dipole array       (input)
    double      *current,       // current[D] -- magnetic dipole current    (input)
    double      *measured,      // measured[D] -- dipole array field        (input)
    double      *sigma2,        // sigma2[D] -- measurement variance        (input)
    int         D               // number of dipoles                        (input)
) {
	FIELD       Sensor;         // test sensor
	double      predicted;      // predicted field for dipole
	double  	dy;             // measured - predicted
	double      chi2;           // goodness of fit
	int         d;              // dipole index

	A6toFIELD(Parms, &Sensor);
    for (d=0, chi2=0.; d<D; d++) {
	    predicted = Coil2B(&Dipole[d], &Sensor, current[d]);
	    dy = measured[d] - predicted;
        chi2 += dy * dy / sigma2[d];
    }
    return (chi2 /  (double)(D - 1));
}

