// Coil2B() -- compute the magnetic field generated by the CTF head coil
//  using the sum of fields of each of ten circular loops
//  All variables are in inscaled SI units
//
//  Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <math.h>
#include <geoms.h>
#include <lfulib.h>


double  Coil2B(             // compute axial & radial magnetic fields
    FIELD   *Coil,          // position & orientation of a CTF head coil
    FIELD   *Sensor,        // position & orientation of sensor
    double  current         // current (amperes)
)
{
    // variables & constants
    double  B[3];           // vector field at sensor
    double  Pnt[3];         // sensor position
    double  radius[] = {    // table of coil turn radii for CTF "head" coil
        3.048000e-03,       // loop 1 -- outer turn radius
        2.790825e-03,       // loop 2
        2.533650e-03,       // loop 3
        2.276475e-03,       // loop 4
        2.019300e-03        // loop 5 -- inner turn radius
    };
    double  sum;            // field projected on sensor vector
    int     i;              // loop radius index
    int     v;              // vector index

    // get sensor position
    for (v=X_; v<=Z_; v++)
        Pnt[v] = Sensor->p[v];

    // sum radial & axial fields of each pair of loops (neglect spacing between pairs)
    for (i=0, sum=0.; i<5; i++) {
        LtoB(Coil, Pnt, radius[i], current, B);
        for (v=X_; v<=Z_; v++)      // sum projection to sensor axis
            sum += B[v] * Sensor->v[v];
    }
    return(2. * sum);
}