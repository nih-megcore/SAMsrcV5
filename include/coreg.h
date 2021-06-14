// coreg.h -- tuneable parameters for SAMcoreg & associated subroutines
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#ifndef COR_H
#define COR_H

// coregistration parameters
#define DEFAULT_ORDER   10          // default order for Nolte spherical harmonics
#define INT_ORDER       3           // sensor coil integration order (3 points)
#define MIN_Z           0.04        // z=4 cm minimum for inclusive mesh elements
#define MIN_SPAN        0.5         // fraction of sorted vertices for minimum eigenvalue ratio for solving moment vector
#define MIN_SNR         .25         // minimum S/N
#define MAX_ANGLE       (M_PI / 3.) // maximum angular error
#define MAX_VERTEX      10000       // sample 10000 random mesh vertices from Freesurfer segmentation
#define MAX_DOT         0.20        // maximum dot product for finding sulcal vertices
#define SMIDGE          1.0e-09     // add 1 nanometre to stopping value to succeed loop test
#define DELTA_TRANS     0.001       // line-search & partial derivative translation delta (m)
#define DELTA_ROTAT     0.002       // line-search & partial derivative rotation delta (mrad)
#define INIT_TRANS      0.002       // initial translation step-size
#define END_TRANS       0.0001      // stopping translations step-size (m)
#define DAMPING         0.5         // damping factpr fpr sequenial line-search
#define NSTEP           8           // line-search steps -8 to +8

// Marquardt-Levenberg tuneable fit parameters
#define DELTA           1.0e-06     // translation & rotation delta
#define ILAMBDA         .001        // initial lambda
#define DLAMBDA         10.         // change in lambda on success of failure

// indices of 6-parameter transform vector
#define aX              0
#define aY              1
#define aZ              2
#define aPSI            3
#define aTHETA          4
#define aPHI            5

#endif
