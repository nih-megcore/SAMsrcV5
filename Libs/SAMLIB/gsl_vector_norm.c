// gsl_vector_norm() -- normalize upper 3 elements of 6-parameter vector
//  This is for a 6-parameter vector where the first three elements are
//  position & the last three are orientation
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <gsl/gsl_vector.h>


void    gsl_vector_norm(
    gsl_vector  *a
)
{
        double  r;
        double  r2;
        int     i;

        // verify that 'a' is at least 6 elements
        if (a->size < 6)
            Cleanup("can't norm orientation");

        // make orientation of 'a' unit vector
        for (i=3, r2=0.; i<6; i++)
            r2 += gsl_vector_get(a, i) * gsl_vector_get(a, i);
        r = sqrt(r2);
        for (i=3; i<6; i++)
            gsl_vector_set(a, i, gsl_vector_get(a, i) / r);
}


double  gsl_vector_norm2(
    gsl_vector  *a
)
{
        double  r;
        double  r2;
        int     i;

        // verify that 'a' is at least 6 elements
        if (a->size < 6)
            Cleanup("can't norm orientation");

        // make orientation of 'a' unit vector
        for (i=3, r2=0.; i<6; i++)
            r2 += gsl_vector_get(a, i) * gsl_vector_get(a, i);
        r = sqrt(r2);
        for (i=3; i<6; i++)
            gsl_vector_set(a, i, gsl_vector_get(a, i) / r);

    // return vector length
    return r;
}


double  gsl_vector_radius(
    gsl_vector  *a
)
{
    double  r2;
    int     i;
    
     // solve radius
     for (i=0, r2=0.; i<3; i++)
        r2 += gsl_vector_get(a, i) * gsl_vector_get(a, i);

    return sqrt(r2);
}

double  gsl_vector_angle(
    gsl_vector  *a
)
{
    double  r;
    double  r2;
    int     i;

    // verify that 'a' is 6 elements
    if (a->size != 6)
        Cleanup("can't compute angle");

    // make orientation of 'a' unit vector
    for (i=3, r2=0.; i<6; i++)
        r2 += gsl_vector_get(a, i) * gsl_vector_get(a, i);
    r = sqrt(r2);
    for (i=3; i<6; i++)
        gsl_vector_set(a, i, gsl_vector_get(a, i) / r);

    // compute phi (angle from z-axis)
    return (atan2(sqrt(gsl_vector_get(a, 0) * gsl_vector_get(a, 0) + gsl_vector_get(a, 1) * gsl_vector_get(a, 1)), gsl_vector_get(a, 2)));
}


