// MeshScatterMS() -- compute Spearman rank correlation and goodness of fit to
//  a linear function between anatomic moment vector & that predicted from
//  the MEG signal.
//
// This subroutine examines all mesh vertices and computes the z^2 values &
//  angular error between the beamformer moment vector Vu & the anatomic moment
//  vector Vc. If a mesh element's unconstrained moment vector is greater than
//  zero (signifying success of the eigensystem solution from 'DIPSolveFwd'),
//  we solve for z2. Compute angular error and save along with z2. The
//  correlation between the angular error and z^2 is output.
//
// Iterative calls to MeshScatter are used by SAMcoreg to refine the coordinates
//  of the MEG primary and reference sensors, relative to an MRI of the brain,
//  under the hypothesis that, as z2 inceases the angular error decreases.
//  Therefore, the best 'fit' for the sensor coordinates should improve correlation
//  & minimize chi^2
//
// This version uses the Multisphere forward solution
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_fit.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <coreg.h>
#include <mesh.h>

#define Order   3               // 3rd-order integration


void MeshScatterMS(
    HeaderInfo  *Header,        // MEG data header
    ChannelInfo *Channel,       // transformed MEG channel info
    gsl_matrix  *C,             // MxM covariance matrix
    gsl_matrix  *Cinv,          // MxM inverse covariance matrix
    double      Noise,          // mean sensor noise power
    MESHINFO    *Mesh,          // Mesh[N] -- cortical mesh vertices and normals
    PARMINFO    *Params,        // coregistration parameters
    int         *nv,            // number of vertices used
    double      *slope,         // slope of the regression
    double      *x2,            // goodness of fit (chi^2)
    double      *corr,          // correlation of fit
    double      *rms            // rms z^2 value
) {
    VOXELINFO   Voxel;          // voxel format for ECDsolve
    gsl_vector  *Bc;            // Bc - M constrained beamformer forward solution
    gsl_vector  *Bu;            // Bu - M unconstrained beamformer forward solution
    gsl_vector  *W;             // W - M beamformer coefficients
    double      *Snr;           // Snr[MAX_VERTEX] -- array of z2 values
    double      *Rad;           // Rad[MAX_VERTEX] -- array of angular error values
    double      *work;          // work[2*MAX_VERTEX] -- workspace for Spearman rank correlation
    double      *Span;          // Span[MAX_VERTEX] -- array of spans for finding median value
    double      s2;             // source power
    double      n2;             // noise power
    double      z2;             // S/N (power)
    double      dot;            // dot product of vu.vc
    double      sx2;            // sum of x^2
    double      span;           // span of eigensystem solution
    double      minSpan;        // minimum span
    double      max;            // maximum z^2
    double      cov00, cov01, cov11;
    double      intercept;
    double      tilt;           // slope, of course
    double      chisq;          // chi^2
    int         M;              // number of primary sensors
    int         m;              // channel index
    int         n;              // mesh vertex index
    int         nn;             // alternate vertex index
    int         v;              // vector index
    int         Nused;          // number of vertices used

    // allocations
    M = C->size1;
    Snr = new_arrayE(double, Params->MaxVertex, "Snr[]");
    Rad = new_arrayE(double, Params->MaxVertex, "Rad[]");
    work = new_arrayE(double, 2 * Params->MaxVertex, "work[]");
    Span = new_arrayE(double, Params->MaxVertex, "Span[]");

    // initialize integration points for all primary & reference sensors w/o changing local spheres
    for (m=0; m<Header->NumChannels; m++)
        if (Channel[m].ChannelType == TYPE_REF_MAG
         || Channel[m].ChannelType == TYPE_REF_GRAD
         || Channel[m].ChannelType == TYPE_MEG
         && Channel[m].Flag == TRUE) {
            FwdSensIntPnt(&Channel[m].Geom.MEGSensor, Order);
        }

    // step through selected vertices & solve for unconstrained moment vector, constrained moment vector, & constrained z^2
#pragma omp parallel private(Bc, Bu, W)
{
    // allocate vectors for forward solutions & weights
    Bc = gsl_vector_alloc(M);
    Bu = gsl_vector_alloc(M);
    W = gsl_vector_alloc(M);

    // for each vertex...
    max = 0.;
#pragma omp for private(Voxel, v, span, s2, n2, z2, dot)
    for (n=0; n<Params->MaxVertex; n++) {

        // set target position
        for (v=X_; v<=Z_; v++) {
            Voxel.p[v] = Mesh[n].Vp[v];         // set vertex position vector
            Voxel.v[v] = Mesh[n].Vc[v];         // set constrained moment vector in vc
            Voxel.Solve = TRUE;                 // solve unconstrained moment vector
        }

        // compute weakly-constrained (optimized) moment vectors & forward solution vectors
        DIPSolveFwd(&Voxel, Bc, Bu, C, Cinv, Header, Channel, &span);
        Mesh[n].span = span;
        for (v=X_; v<=Z_; v++)                  // save beamformer (unconstrained) moment vector in Vu
            Mesh[n].Vu[v] = Voxel.v[v];         // Voxel.v is modified on output to contain unconstrained moment vector
        
        // solve [constrained] z^2 & angular error
        SAMsolve(C, Cinv, Bc, W, Noise, &s2, &n2);
        z2 = (n2 > 0.)? s2 / n2: 0.;
#pragma omp critical
{
        if (z2 > max)
            max = z2;
}
        Mesh[n].z2 = z2;

        // compute angular error (acos|Vu<dot>Vc|), & save in 'Mesh'
        for (v=X_, dot=0.; v<=Z_; v++)
            dot += Mesh[n].Vu[v] * Mesh[n].Vc[v];
        Mesh[n].angle = acos(fabs(dot));
    }

    // free GSL vectors
    gsl_vector_free(Bc);
    gsl_vector_free(Bu);
    gsl_vector_free(W);

}   // end parallel

    // generate Span vector
    for (n=0; n<Params->MaxVertex; n++)
        Span[n] = Mesh[n].span;
    gsl_sort(Span, 1, Params->MaxVertex);   // sort Span in ascending order
    n = (int)(Params->MinSpan * (double)Params->MaxVertex);
    minSpan = Span[n];
#ifdef EBUG
    printf("\ntheshold span=%f\n", minSpan);
    fflush(stdout);
#endif

    // finally, limit mesh solutions to constraints & form Rad & Snr vectors and accumulate sx2
    for (n=nn=0, sx2=0.; n<Params->MaxVertex; n++)
        if (Mesh[n].span >= minSpan && Mesh[n].z2 > max * Params->MinSNR && Mesh[n].angle <= Params->MaxAngle) {
            Snr[nn] = Mesh[n].z2;
            Rad[nn] = Mesh[n].angle;
            sx2 += Snr[nn] * Snr[nn];
            nn++;
        } else {
            Mesh[n].z2 = -1.;
            Mesh[n].angle = -1.;
        }
    Nused = nn;

    // test for Nused >= 100
    if (Nused >= 100) {

        // output rms z^2 values
        *rms = sqrt(sx2 / (double)Nused);
        gsl_fit_linear(Rad, 1, Snr, 1, Nused, &intercept, &tilt, &cov00, &cov01, &cov11, &chisq);
        *slope = -tilt;
        *x2 = chisq / (double)Nused;
        *corr = -kendall(Rad, Snr, Nused);
//      *corr = -gsl_stats_spearman(Rad, 1, Snr, 1, Nused, work);
    } else {    // this should never happen!
        *slope = 0.;
        *rms = 0.;
        *x2 = -1.;
        *corr = 0.;
    }
    *nv = Nused;

    // free allocations
    free(Snr);
    free(Rad);
    free(work);
    free(Span);
}
