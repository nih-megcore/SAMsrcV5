// EtoR() -- compute rotation matrix from Euler angles.
//  In SAMcoreg, the order of evaluation of Euler angles is psi, theta, phi
//
// Some definitions:
//  Rotation[X_][X_] = c[PSI_] * c[PHI_] - c[THETA_] * s[PHI_] * s[PSI_];
//  Rotation[X_][Y_] = c[PSI_] * s[PHI_] + c[THETA_] * c[PHI_] * s[PSI_];
//  Rotation[X_][Z_] = s[PSI_] * s[THETA_];
//  Rotation[Y_][X_] = - (s[PSI_] * c[PHI_]) - c[THETA_] * s[PHI_] * c[PSI_];
//  Rotation[Y_][Y_] = - (s[PSI_] * s[PHI_]) + c[THETA_] * c[PHI_] * c[PSI_];
//  Rotation[Y_][Z_] = c[PSI_] * s[THETA_];
//  Rotation[Z_][X_] = s[THETA_] * s[PHI_];
//  Rotation[Z_][Y_] = - (s[THETA_] * c[PHI_]);
//  Rotation[Z_][Z_] = c[THETA_];
//
//  R=
//      cosθ cosψ       sinφ sinθ cosψ − cosφ sinψ      cosφ sinθ cosψ + sinφ sinψ
//      cosθ sinψ       sinφ sinθ sinψ + cosφ cosψ      cosφ sinθ sinψ - sinφ cosψ
//      − sinθ          sinφ cosθ                       cosφ cosθ
//
//  θ = −arcsin(rzx)
//  φ = atan2(rzy, rzz)
//  ψ = atan2(ryx, rxx)
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <geoms.h>

void EtoR(
    double  Euler[3],           // Euler angles
    double  Rotation[3][3]      // rotation matrix
) {
    double          c[3];       // cosines
    double          s[3];       // sines
    register int    v;          // Euler angle index

    // compute trig functions of Euler angles
    for (v=PSI_; v<=PHI_; v++) {
        s[v] = sin(Euler[v]);
        c[v] = cos(Euler[v]);
    }

    // rotation matrix definitions:
    Rotation[X_][X_] = c[THETA_] * c[PSI_];
    Rotation[X_][Y_] = c[PHI_] * s[PSI_] + s[PHI_] * s[THETA_] * c[PSI_];
    Rotation[X_][Z_] = s[PHI_] * s[PSI_] - (c[PHI_] * s[THETA_] * c[PSI_]);
    Rotation[Y_][X_] = - (c[THETA_] * s[PSI_]);
    Rotation[Y_][Y_] = c[PHI_] * c[PSI_] - (s[PHI_] * s[THETA_] * s[PSI_]);
    Rotation[Y_][Z_] = s[PHI_] * c[PSI_] + c[PHI_] * s[THETA_] * s[PSI_];
    Rotation[Z_][X_] = s[THETA_];
    Rotation[Z_][Y_] = - (s[PHI_] * c[THETA_]);
    Rotation[Z_][Z_] = c[PHI_] * c[THETA_];
}
