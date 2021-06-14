// Levy distribution function
//
//  Author: Stephen E. Robinon
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <time.h>
#include <samlib.h>

double  Levy(       // returns p-value
    double  x,      // domain
    double  mu,     // location parameter
    double  c       // scale factor (peakiness)
)

{
    double  den;    // denominator
    double  num;    // numerator
    double  xmu;    // x-mu
    double  scale;  // scale factor
 
    // sanity test
    if (x < mu)
        Cleanup("x must be >= mu");

    // compute probability
    xmu = x - mu;
    scale = sqrt(c / (2. * M_PI));
    den = pow(xmu, 1.5);
    num = exp(-c / (2. * xmu));
    return (scale * num / den);
/*
    // seed random number generaator using system time
    srand48((long)time((time_t *)0));

    // return random Levy value
    return(c / pow(erfc(1. - drand48()), 2.) + mu);
*/
}
