// Vector_FIELD_convert() -- copy linear array to FIELD parameters (position, orientation)
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <geoms.h>

void    A6toFIELD (
    double  *a6,
    FIELD   *field
)
{
    field->p[X_] = a6[0];
    field->p[Y_] = a6[1];
    field->p[Z_] = a6[2];
    field->v[X_] = a6[3];
    field->v[Y_] = a6[4];
    field->v[Z_] = a6[5];
}

void    FIELDtoA6 (
    FIELD   *field,
    double  *a6
)
{
    a6[0] = field->p[X_];
    a6[1] = field->p[Y_];
    a6[2] = field->p[Z_];
    a6[3] = field->v[X_];
    a6[4] = field->v[Y_];
    a6[5] = field->v[Z_];
}
    
