// Compute integration points and coefficients for the sensor coils prior to calling FwdSensor().

#include <stdio.h>
#include <math.h>
#include <geoms.h>
#include <samlib.h>

void FwdSensIntPnt(
    SENSOR      *Sensor,    // sensor channel (magnetometer, gradiometer) structure
    int         Order       // integration order: 1, 2, 3, 4, 5, & 7 are implemented
) {
    int         i;          // integration index
    int         c;          // coil index
    int         v;          // vector index
    double      m;          // |m|
    double      m2;         // m^2
    double      Vx[3];      // unit x-vector
    double      Vy[3];      // unit y-vector
    double      Radius;     // radius of integration point
    double      XDisp;      // x position in plane of coil
    double      YDisp;      // y position in plane of coil
    double      angle;      // angle of integration point

    switch (Order) {
        case 1:
            Sensor->Coil[0].NumInt = 1;
            Sensor->Coil[0].w[0] = 1.;
            break;
        case 2:
            Sensor->Coil[0].NumInt = 3;
            Sensor->Coil[0].w[0] = 0.333333;
            Sensor->Coil[0].w[1] = 0.333333;
            Sensor->Coil[0].w[2] = 0.333333;
            break;
        case 3:
            Sensor->Coil[0].NumInt = 4;
            Sensor->Coil[0].w[0] = 0.25;
            Sensor->Coil[0].w[1] = 0.25;
            Sensor->Coil[0].w[2] = 0.25;
            Sensor->Coil[0].w[3] = 0.25;
            break;
        case 4:
            Sensor->Coil[0].NumInt = 6;
            Sensor->Coil[0].w[0] = 0.25;
            Sensor->Coil[0].w[1] = 0.25;
            Sensor->Coil[0].w[2] = 0.25;
            Sensor->Coil[0].w[3] = 0.25;
            Sensor->Coil[0].w[4] = 0.25;
            Sensor->Coil[0].w[5] = 0.15;
            break;
        case 5:
            Sensor->Coil[0].NumInt = 7;
            Sensor->Coil[0].w[0] = 0.125;
            Sensor->Coil[0].w[1] = 0.125;
            Sensor->Coil[0].w[2] = 0.125;
            Sensor->Coil[0].w[3] = 0.125;
            Sensor->Coil[0].w[4] = 0.125;
            Sensor->Coil[0].w[5] = 0.125;
            Sensor->Coil[0].w[6] = 0.25;
            break;
        case 7:     // 7th-order integration -- set 12-point integration weights
            Sensor->Coil[0].NumInt = 12;
            Sensor->Coil[0].w[0] = 0.12321;
            Sensor->Coil[0].w[1] = 0.12321;
            Sensor->Coil[0].w[2] = 0.12321;
            Sensor->Coil[0].w[3] = 0.12321;
            Sensor->Coil[0].w[4] = 0.074074;
            Sensor->Coil[0].w[5] = 0.074074;
            Sensor->Coil[0].w[6] = 0.074074;
            Sensor->Coil[0].w[7] = 0.074074;
            Sensor->Coil[0].w[8] = 0.052715;
            Sensor->Coil[0].w[9] = 0.052715;
            Sensor->Coil[0].w[10] = 0.052715;
            Sensor->Coil[0].w[11] = 0.052715;
            break;
        default:
            fatalerr("invalid Order for forward solution");
    }

    // evaluate all coil integration points
    for (c=0; c<Sensor->NumCoils; c++) {

        // compulsively make certain that the coil vector is unit
        for (v=X_, m2=0.; v<=Z_; v++)
            m2 += Sensor->Coil[c].origin.v[v] * Sensor->Coil[c].origin.v[v];
        m = sqrt(m2);
        for (v=X_; v<=Z_; v++)
            Sensor->Coil[c].origin.v[v] /= m;   // there now... I feel better already!

        // if req'd, solve for coil plane from coil vector
        if (Order > 1)
            OrthoPlane(Sensor->Coil[c].origin.v, Vx, Vy);

        // test for valid integration order & set up integration points
        switch (Order) {
            case 1:     // single point
                for (v=X_; v<=Z_; v++)
                    Sensor->Coil[c].B[0].p[v] = Sensor->Coil[c].origin.p[v];
                break;

            case 2:     // three integration points
                Radius = M_SQRT1_2 * Sensor->Coil[c].radius;
                for (i=0; i<3; i++) {
                    angle = (double)i * 2. * M_PI / 3.;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v] + Vx[v] * XDisp + Vy[v] * YDisp;
                }
                break;

            case 3:     // four integration points
                Radius = M_SQRT1_2 * Sensor->Coil[c].radius;
                for (i=0; i<4; i++) {
                    angle = (double)i * M_PI / 2.;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v] + Vx[v] * XDisp + Vy[v] * YDisp;
                }
                break;

            case 4:     // six integration points
                Radius = 0.81650 * Sensor->Coil[c].radius;  // 1st five integration points
                for (i=0; i<5; i++) {
                    angle = (double)i * 2. * M_PI / 5.;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v] + Vx[v] * XDisp + Vy[v] * YDisp;
                }
                for (v=X_; v<=Z_; v++)              // 6th integration point
                    Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v];
                break;

            case 5:     // seven integration points
                Radius = 0.81650 * Sensor->Coil[c].radius;  // 1st six integration points
                for (i=0; i<6; i++) {
                    angle = (double)i * M_PI / 3.;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v] + Vx[v] * XDisp + Vy[v] * YDisp;
                }
                for (v=X_; v<=Z_; v++)              // 7th integration point
                    Sensor->Coil[c].B[6].p[v] = Sensor->Coil[c].origin.p[v];
                break;

            case 7:     // twelve integration points
                Radius = 0.45667 * Sensor->Coil[c].radius;  // 1st four integration points
                for (i=0; i<4; i++) {
                    angle = ((double)i - 0.5) * M_PI / 2.;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v] + Vx[v] * XDisp + Vy[v] * YDisp;
                }
                Radius = 0.86603 * Sensor->Coil[c].radius;  // 2nd four integration points
                for (i=4; i<8; i++) {
                    angle =(double)i * M_PI / 2.;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v] + Vx[v] * XDisp + Vy[v] * YDisp;
                }
                Radius = 0.91100 * Sensor->Coil[c].radius;  // 3rd four integration points
                for (i=8; i<12; i++) {
                    angle = ((double)i - 0.5) * M_PI / 2.;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Sensor->Coil[c].B[i].p[v] = Sensor->Coil[c].origin.p[v] + Vx[v] * XDisp + Vy[v] * YDisp;
                }
                break;

            default:
                break;
        }
        Sensor->Coil[c].Evaluated = TRUE;
        if (c > 0) {
            Sensor->Coil[c].NumInt = Sensor->Coil[0].NumInt;
            for (i=0; i<Sensor->Coil[c].NumInt; i++)
                Sensor->Coil[c].w[i] = Sensor->Coil[0].w[i];
        }
    }
}
