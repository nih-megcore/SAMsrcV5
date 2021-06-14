#include <lfulib.h>
#include <math.h>

void OrthoPlane(
    double  *Vr,
    double  *Vx,
    double  *Vy
) {
    double      r;          // distance
    double      r2;         // r^2
    double      max;        // axis with maximum projection
    int         v;          // vector index
    int         index;      // axis index
    static int  ref = X_;   // initial axis

    // compulsively make certain that Vr is unit length
    for (v=X_, r2=0.; v<=Z_; v++)
        r2 += Vr[v] * Vr[v];
    r = sqrt(r2);
    for (v=X_; v<=Z_; v++)
        Vr[v] /= r;         // there now... I feel better already!

    // find principal coil vector
    for (v=X_, max=0.; v<=Z_; v++)
        if (fabs(Vr[v]) > max) {
            max = fabs(Vr[v]);
            index = v;
        }

    // make hysteretic decision of x, y or z reference
    if (fabs(Vr[ref]) < 0.1)
        ref = index;

    // compute x' vector perpendicular to Vr
    switch (ref) {
        case X_:        // x' = r <cross> Z_
            Vx[X_] = Vr[Y_];
            Vx[Y_] = -Vr[X_];
            Vx[Z_] = 0.;
            break;
        case Y_:        // x' = r <cross> X_
            Vx[X_] = 0.;
            Vx[Y_] = Vr[Z_];
            Vx[Z_] = -Vr[Y_];
            break;
        case Z_:        // x' = r <cross> Y_
            Vx[X_] = -Vr[Z_];
            Vx[Y_] = 0.;
            Vx[Z_] = Vr[X_];
            break;
    }

    // make x' a unit vector
    for (v=X_, r2=0.; v<=Z_; v++)
        r2 += Vx[v] * Vx[v];
    r = sqrt(r2);
    for (v=X_; v<=Z_; v++)
        Vx[v] /= r;

    // y' = v <cross> x'
    Vy[X_] = Vr[Y_] * Vx[Z_] - Vr[Z_] * Vx[Y_];
    Vy[Y_] = Vr[Z_] * Vx[X_] - Vr[X_] * Vx[Z_];
    Vy[Z_] = Vr[X_] * Vx[Y_] - Vr[Y_] * Vx[X_];

    // make y' a unit vector, too
    for (v=X_, r2=0.; v<=Z_; v++)
        r2 += Vy[v] * Vy[v];
    r = sqrt(r2);
    for (v=X_; v<=Z_; v++)
        Vy[v] /= r;
}
