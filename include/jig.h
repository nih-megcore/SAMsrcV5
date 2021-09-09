#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <geoms.h>

#define P           6                           // number of parameters in fit
#define N_OPM       19                          // 16 primary sensors + 3 reference magnetometers
#define N_DIP       37                          // 37 CTF magnetic dipole coils
#define L           6.957e-07                   // CTF magnetic dipole coil inductance (H)
#define X_COND_GAIN 2.87e-04                    // 287 µA/Vin -- transconductance gain (Acoil/Vin)
#define X_IMP_GAIN  2.0e+03                     // 2 mVout/µA -- transimpedance gain (Vout/Acoil)
#define BOX_GAIN    (X_COND_GAIN * X_IMP_GAIN)  // Vout/Vin -- 0.574
#define V2A         (1. / X_IMP_GAIN)           // 500 µA/Vout
#define ADGAIN      2.980232416e-07             // ADC gain factor (V/bit)
#define AM2         2.0583e-04                  // Ampere-meters^2 per Ampere coil current for CTF magnetic dipole coil (sum of loop areas)

// OPMstats parameters (not req'd)
#define FREQ        21.0                        // nominal dipole frequency (tbd, based on OPM spectrum)
#define SRATE       1000.0                      // nominal sample rate

// Levenberg-Marquardt parameters

#define INIT_LAMBDA 1.0e-03     // initial lambda
#define END_LAMBDA  1.0e-10
#define MAX_LAMBDA  1.0e+10
//#define LAMBDA_UP   11.
//#define LAMBDA_DN   9.

//#define INIT_LAMBDA 1.0e+08
//#define MAX_LAMBDA  1.0e+15
//#define END_LAMBDA  1.0e-09
#define LAMBDA_UP   10.
#define LAMBDA_DN   10.

#define MAXITS      100

#define MAXTEST     3           // convergence test
#define DCHSQ       .01         // minimum change in chi-square
#define EPSILON     1.0e-06     // change in chi-square that is significant
#define MIN_OPM_RAD 0.068
#define MAX_OPM_RAD 0.150
#define MAX_OPM_ORI (M_PI/5.)

#define MAXROT      0.087  	// maximum allowable angular deviation of sensor orientation vector, in radians
#define MAXTRANS    0.005       // maximum allowable translation of sensor location variable, in meters 
#define MINGAIN     0.5        // minimum gain 
#define MAXGAIN     1.5        // minimum gain 

#define MTOB    TRUE            // switch from Coil2B to MTOB (must recompile SAMLIB)




