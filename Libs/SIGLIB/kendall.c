// Kendall -- compute Kendall's tau (nonparametric rank correlation)
//
// Given data arrays x & y, this program returns Kendall’s τ.
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>


double  kendall(
    double          *x,
    double          *y,
    unsigned long   n
)
{
    double          a1;
    double          a2;
    double          aa;
    double          tau;
    unsigned long   n1;
    unsigned long   n2;
    unsigned long   k;
    unsigned long   j;
    long            is;

    for (j=0, n1=n2=is=0; j<(n-1); j++) {   // loop over first member of pair
        for (k=(j+1); k<n; k++) {           // & second member
            a1 = x[j] - x[k];
            a2 = y[j] - y[k];
            aa = a1 * a2;
            if (aa) {                       // neither array has a tie
                ++n1;
                ++n2;
                if (aa > 0.)
                    ++is;
                else
                    --is;
            } else {                        // one or both arrays have ties
                if (a1)
                    ++n1;                   // an "extra x" event
                if (a2)
                    ++n2;                   // an "extra y" event
            }
        }
    }
    tau = (double)is / (sqrt((double)n1) * sqrt((double)n2));
    return(tau);
}
