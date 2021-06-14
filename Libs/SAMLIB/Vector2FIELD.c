// Vector2FIELD() -- copy linear array to FIELD parameters (position, orientation)
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <geoms.h>
#include <gsl/gsl_vector.h>


void    Vector2FIELD(
    gsl_vector  *a,
    FIELD       *parms
)

{
        parms->p[X_] = gsl_vector_get(a, 0);
        parms->p[Y_] = gsl_vector_get(a, 1);
        parms->p[Z_] = gsl_vector_get(a, 2);
        parms->v[X_] = gsl_vector_get(a, 3);
        parms->v[Y_] = gsl_vector_get(a, 4);
        parms->v[Z_] = gsl_vector_get(a, 5);
}
