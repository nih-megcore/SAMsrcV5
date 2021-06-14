// MoveField() -- takes a transformation matrix & a vector field,
//  & moves the vector field to new reference frame
//
//  Author: Stephen E. Robinson
//          Neuromagnetism Laboratory
//          Dept.of Neurology
//          Henry Ford Hospital
//          Detroit, MI
//          Copyright (c) 2007
//

#include <geoms.h>

void MoveField(
    FIELD   *Field,     // vector field to be moved
    XFORM   *Transform  // transformation matrix
) {
    FIELD   A;
    int     i, j;

    for (i=X_; i<=Z_; i++) {
        for (j=X_, A.p[i]=A.v[i]=0.; j<=Z_; j++) {
            A.p[i] += Transform->Rotate[i][j] * Field->p[j];
            A.v[i] += Transform->Rotate[i][j] * Field->v[j];
        }
        A.p[i] += Transform->Translate[i];
    }

    // replace original vector FIELD with new one
    for (i=X_; i<=Z_; i++) {
        Field->p[i] = A.p[i];
        Field->v[i] = A.v[i];
    }
}
