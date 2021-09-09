// Vector2OPMFIELD() -- copy linear array to FIELD parameters (position, orientation)
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <geoms.h>
#include <gsl/gsl_vector.h>


void    Vector2OPMFIELD(
    gsl_vector  *a,
    OPMFIELD       *parms
)

{
        parms->p[X_] = gsl_vector_get(a, 0);
        parms->p[Y_] = gsl_vector_get(a, 1);
        parms->p[Z_] = gsl_vector_get(a, 2);
        parms->v[0] = gsl_vector_get(a, 3);
        parms->v[1] = gsl_vector_get(a, 4);
        parms->g = gsl_vector_get(a, 5);
}
