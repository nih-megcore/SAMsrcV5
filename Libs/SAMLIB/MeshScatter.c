// MeshScatter() -- compute Spearman rank correlation and goodness of fit to
//  a linear function between anatomic moment vector & that predicted from
//  the MEG signal.
//
// This subroutine examines all mesh vertices and computes the z^2 values &
//  angular error between the beamformer moment vector Vu & the anatomic moment
//  vector Vc. If a mesh element's unconstrained moment vector is greater than
//  zero (signifying success of the eigensystem solution from 'ECDSolveFwd'),
//  we solve for z2. Compute angular error and save along with z2. The
//  correlation between the angular error and z^2 is output.
//
// Iterative calls to MeshScatter are used by SAMcoreg to refine the coordinates
//  of the MEG primary and reference sensors, relative to an MRI of the brain,
//  under the hypothesis that, as z2 inceases the angular error decreases.
//  Therefore, the best 'fit' for the sensor coordinates should improve correlation
//  & minimize chi^2
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <coreg.h>
#include <mesh.h>


void MeshScatter(
    HeaderInfo  *Header,        // MEG data header
    ChannelInfo *NewChannel,    // transformed MEG channel info
    gsl_matrix  *C,             // MxM covariance matrix
    gsl_matrix  *Cinv,          // MxM inverse covariance matrix
    double      Noise,          // mean sensor noise power
    HULL        *Hull,          // conductive model with normals
    MESHINFO    *Mesh,          // Mesh[N] -- cortical mesh vertices and normals
    PARMINFO    *Params,        // coregistration parameters
    int         *nv,            // number of vertices used
    double      *slope,         // slope of the regression
    double      *x2,            // goodness of fit (chi^2)
    double      *corr,          // correlation of fit
    double      *rms            // rms z^2 value
) {
    VOXELINFO   Voxel;          // voxel format for ECDsolve
    COEFFS      *coeff;         // coefficient structure for sensor frame
    gsl_vector  *Bc;            // Bc - M constrained beamformer forward solution
    gsl_vector  *Bu;            // Bu - M unconstrained beamformer forward solution
    gsl_vector  *W;             // W - M beamformer coefficients
    double      *Snr;           // Snr[MAX_VERTEX] -- array of z2 values
    double      *Rad;           // Rad[MAX_VERTEX] -- array of angular error values
    double      *Wgt;           // Wgt[MAX_VERTEX] -- array of weights
    double      *work;          // work[2*MAX_VERTEX] -- workspace for Spearman rank correlation
    double      s2;             // source power
    double      n2;             // noise power
    double      z2;             // S/N (power)
    double      dot;            // dot product of vu.vc
    double      span;           // span of eigensystem solution
    double      max;            // maximum z^2
    double      tmp;            // 'ol reliable...
    double      cov00, cov01, cov11;
    double      intercept;
    double      tilt;           // slope, of course
    double      chisq;
    double      ratio;          // proportion for weight
    double      sx;             // sum of x
    double      sx2;            // sum of x^2
    static double   sigma2 = -1.;   // variance at ~pi/2
    int         M;              // number of primary sensors
    int         n;              // mesh vertex index
    int         nn;             // alternate vertex index
    int         v;              // vector index
    int         Nused;          // number of vertices used

    // allocations
    M = C->size1;
    Snr = new_arrayE(double, Params->MaxVertex, "Snr[]");
    Rad = new_arrayE(double, Params->MaxVertex, "Rad[]");
    Wgt = new_arrayE(double, Params->MaxVertex, "Wgt[]");
    work = new_arrayE(double, 2 * Params->MaxVertex, "work[]");
    coeff = new(COEFFS);

    // (1) reevaluate basis coefficients & sensor geometry in Nolte solution
    ECDIntPnt(Header, NewChannel, Hull, coeff);

    // (2) step through selected vertices & solve for unconstrained moment vector, constrained moment vector, & constrained z^2
#pragma omp parallel private(Bc, Bu, W)
{
    // allocate vectors for forward solutions & weights
    Bc = gsl_vector_alloc(M);
    Bu = gsl_vector_alloc(M);
    W = gsl_vector_alloc(M);

    max = 0.;
//    min = HUGE;
#pragma omp for private(Voxel, v, span, tmp, s2, n2, z2, dot)
    for (n=0; n<Params->MaxVertex; n++) {

        // (2a) set target position
        for (v=X_; v<=Z_; v++) {
            Voxel.p[v] = Mesh[n].Vp[v];         // set vertex position vector
            Voxel.v[v] = Mesh[n].Vc[v];         // set constrained moment vector in vc
            Voxel.Solve = TRUE;                 // solve unconstrained moment vector
        }

        // (2b) compute weakly-constrained (optimized) moment vectors & forward solution vectors
        ECDSolveFwd0(&Voxel, Bc, Bu, C, Cinv, Header, NewChannel, Hull, &span, coeff);
        if (span >= Params->MinSpan) {           // span tests reliability of unconstrained moment vector
            for (v=X_, tmp=0.; v<=Z_; v++) {    // save beamformer (unconstrained) moment vector in Vu
                Mesh[n].Vu[v] = Voxel.v[v];     // Voxel.v is modified on output to contain unconstrained moment vector
                tmp += Voxel.v[v] * Voxel.v[v];
            }
        } else {
            tmp = 0.;
        }

        // (2c) solve [constrained] z^2 & angular error, if unconstrained vertex normal vector is non-zero
        if (tmp > 0.) {

            // (2d) solve constrained beamformer output S/N (z2), sum mean square S/N, & save in 'Mesh'
            SAMsolve(C, Cinv, Bc, W, Noise, &s2, &n2);
            z2 = (n2 > 0.)? s2 / n2: 0.;
#pragma omp critical
{
            if (z2 > max)
                max = z2;
}
            Mesh[n].z2 = z2;

            // (2e) compute angular error (acos|Vu<dot>Vc|), & save in 'Mesh'
            for (v=X_, dot=0.; v<=Z_; v++)
                dot += Mesh[n].Vu[v] * Mesh[n].Vc[v];
            Mesh[n].angle = acos(fabs(dot));

        } else {
            Mesh[n].z2 = 0.;
            Mesh[n].angle = -1.;
        }
    }

    // free GSL vectors
    gsl_vector_free(Bc);
    gsl_vector_free(Bu);
    gsl_vector_free(W);

}   // end parallel

    // compute initial variance estimate
    if (sigma2 < 0.) {
        for (n=nn=0, sx=sx2=0.; n<Params->MaxVertex; n++)
            if (Mesh[n].angle > 0.9 * M_PI_2 && Mesh[n].z2 > 0.) {
                sx += Mesh[n].angle;
                sx2 += Mesh[n].angle * Mesh[n].angle;
                nn++;
            }
        sigma2 = (sx2 - sx * sx / (double)nn) / (double)nn;
    }

    // (3) form Rad & Snr vectors and accumulate sx2
    for (n=nn=0, sx2=0.; n<Params->MaxVertex; n++)
        if (Mesh[n].z2 > max * Params->MinSNR) {
            Snr[nn] = Mesh[n].z2;
            Rad[nn] = Mesh[n].angle;
            sx2 += Snr[nn] * Snr[nn];
            nn++;
        }
    Nused = nn;

    // form Wgt vector
    for (n=1; n<=Nused; n++) {
        ratio = (double)n / (double)Nused;
        Wgt[n-1] = 1. / (ratio * sigma2);
    }
    
    // (4) test for Nused >= 100
    if (Nused >= 100) {

        // (4a) output rms z^2 values
        *rms = sqrt(sx2 / (double)Nused);
        gsl_fit_wlinear(Rad, 1, Wgt, 1, Snr, 1, Nused, &intercept, &tilt, &cov00, &cov01, &cov11, &chisq);
        *slope = -tilt;
        *x2 = chisq;
        *corr = -gsl_stats_spearman(Rad, 1, Snr, 1, Nused, work);
    } else {    // this should never happen!
        *slope = 0.;
        *rms = 0.;
        *x2 = -1.;
        *corr = 0.;
    }
    *nv = Nused;

    free(Snr);
    free(Rad);
    free(Wgt);
    free(work);
    ECDIntPntFree(coeff);
    free(coeff);
}
