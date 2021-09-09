// OPMFIELD2Vector() -- copy FIELD parameters (position, orientation) to linear array
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <geoms.h>
#include <gsl/gsl_vector.h>


void    OPMFIELD2Vector(
    OPMFIELD       *parms,
    gsl_vector  *a
)

{
        gsl_vector_set(a, 0, parms->p[X_]);
        gsl_vector_set(a, 1, parms->p[Y_]);
        gsl_vector_set(a, 2, parms->p[Z_]);
        gsl_vector_set(a, 3, parms->v[RHO_]);
        gsl_vector_set(a, 4, parms->v[THETA_]);
        gsl_vector_set(a, 5, parms->g);
}
