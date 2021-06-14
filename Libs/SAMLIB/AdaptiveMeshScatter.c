// AdaptiveMeshScatter() -- compute correlayion & chi^2 for z^2 vs angular error using
//  anatomic moment vector & that predicted from the MEG signal.
//
// This subroutine examines all mesh vertices and computes the z^2 values &
//  angular error between the beamformer moment vector Vu & the anatomic moment
//  vector Vc. If a mesh element's unconstrained moment vector is greater than
//  zero (signifying success of the eigensystem solution from 'ECDSolveFwd'),
//  we solve for z2. Compute angular error and save along with z2. The
//  correlation between the angular error and z^2 is output.
//
// Iterative calls to AdaptiveMeshScatter are used by SAMcoreg to refine the coordinates
//  of the MEG primary and reference sensors, relative to an MRI of the brain,
//  under the hypothesis that, as z2 inceases the angular error decreases.
//  Therefore, the best 'fit' for the sensor coordinates should minimize chi^2
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

void AdaptiveMeshScatter(
    HeaderInfo  *Header,        // MEG data header
    ChannelInfo *NewChannel,    // transformed MEG channel info
    gsl_matrix  *C,             // MxM covariance matrix
    gsl_matrix  *Cinv,          // MxM inverse covariance matrix
    double      Noise,          // mean sensor noise power
    HULL        *Hull,          // conductive model with normals
    MESHINFO    *Mesh,          // Mesh[N] -- cortical mesh vertices and normals
    PARMINFO    *Parms,         // input vertex count & output MaxAngle & MinSNR
    int         *nv,            // number of vertices used
	double      *x2,			// goodness of fit (chi^2)
	double		*corr,			// Spearman correlation coefficient
    double      *rms            // rms z^2 value
) {
    VOXELINFO   Voxel;          // voxel format for ECDsolve
    COEFFS      *coeff;         // coefficient structure for sensor frame
    gsl_vector  *Bc;            // Bc - M constrained beamformer forward solution
    gsl_vector  *Bu;            // Bu - M unconstrained beamformer forward solution
    gsl_vector  *W;             // W - M beamformer coefficients
    double      *Snr;           // Snr[MAX_VERTEX] -- array of z2 values
    double      *Rad;           // Rad[MAX_VERTEX] -- array of angular error values
    double      *work;          // work[2*MAX_VERTEX] -- workspace for Spearman rank correlation
    double      s2;             // source power
    double      n2;             // noise power
    double      z2;             // S/N (power)
    double      dot;            // dot product of vu.vc
    double      sx2;            // sum of squares
    double      fract;          // fraction
    double      max;            // fraction of maximum z2 for entire mesh
    double      z2max;          // maximum z2 for entire mesh
    double      angle;          // maximum angular error
    double      r;              // correlation coefficient
    double      rbest;          // largest correlation
    double      bestAngle;
    double      bestSnr;
    double      span;           // span of eigensystem solution
    double      tmp;            // 'ol reliable...
	double		cov00, cov01, cov11;
	double		c0, c1;
	double		chisq;
    int         M;              // number of primary sensors
    int         n;              // mesh vertex index
    int         nn;             // alternate vertex index
    int         v;              // vector index
    int         Nused;          // number of vertices used

    // allocations
    M = C->size1;
    Snr = new_arrayE(double, Parms->MaxVertex, "Snr[]");
    Rad = new_arrayE(double, Parms->MaxVertex, "Rad[]");
    work = new_arrayE(double, 2 * Parms->MaxVertex, "work[]");

    coeff = new(COEFFS);

    // (1) reevaluate basis coefficients & sensor geometry in Nolte solution
    ECDIntPnt(Header, NewChannel, Hull, coeff);

    // (2) step through selected vertices & solve for unconstrained moment vector, constrained moment vector, & constrained z^2
    z2max = 0.;
#pragma omp parallel private(Bc, Bu, W)
{
    Bc = gsl_vector_alloc(M);
    Bu = gsl_vector_alloc(M);
    W = gsl_vector_alloc(M);

#pragma omp for private(Voxel, v, span, tmp, s2, n2, z2, dot)
    for (n=0; n<Parms->MaxVertex; n++) {

        // (2a) set target position
        for (v=X_; v<=Z_; v++) {
            Voxel.p[v] = Mesh[n].Vp[v];     // set vertex position vector
            Voxel.v[v] = Mesh[n].Vc[v];     // set constrained moment vector in vc
            Voxel.Solve = TRUE;
        }

        // (2b) compute weakly-constrained (optimized) moment vectors & forward solution vectors
        ECDSolveFwd0(&Voxel, Bc, Bu, C, Cinv, Header, NewChannel, Hull, &span, coeff);
        for (v=X_, tmp=0.; v<=Z_; v++) {    // save beamformer (unconstrained) moment vector in vu
            Mesh[n].Vu[v] = Voxel.v[v];     // Voxel.v is modified on output to contain unconstrained moment vector
            tmp += Voxel.v[v] * Voxel.v[v];
        }

        // (2c) test for non-zero unconstrained vertex normal vector
        if (tmp > 0.) {

            // (2d) solve constrained beamformer output S/N (z2) & angular error -- save in 'Mesh'
            SAMsolve(C, Cinv, Bc, W, Noise, &s2, &n2);
            z2 = (n2 > 0.)? s2 / n2: 0.;
#pragma omp critical
            if (z2 > z2max)
                z2max = z2;
            Mesh[n].z2 = z2;

            // (2e) compute angular error (acos|Vu<dot>Vc|), & save in 'Mesh'
            for (v=X_, dot=0.; v<=Z_; v++)
                dot += Mesh[n].Vu[v] * Mesh[n].Vc[v];
            Mesh[n].angle = acos(fabs(dot));

        } else {
            Mesh[n].z2 = 0.;
            Mesh[n].angle = 0.;
        }
    }

    gsl_vector_free(Bc);
    gsl_vector_free(Bu);
    gsl_vector_free(W);

} // end parallel

    // (3) adaptive test of limits (best correlation)
//  printf("\nadaptive limits search:\n");
//  fflush(stdout);
    for(fract=0.1, rbest=0.; fract<=0.60; fract+=0.01) {
        max = fract * z2max;
        for(angle=0.25; angle<=1.70; angle+=0.01) {
            for(n=nn=0; n<Parms->MaxVertex; n++)
                if(Mesh[n].z2 >= max && Mesh[n].angle <= angle) {
                    Snr[nn] = Mesh[n].z2;
                    Rad[nn] = Mesh[n].angle;
                    nn++;
                }
            r = fabs(gsl_stats_correlation(Snr, 1, Rad, 1, nn));
//          printf("r=%f, fract=%f, angle=%f, n=%d\n", r, fract, angle, nn);
//          fflush(stdout);
            if(r > rbest) {
                bestAngle = angle;
                bestSnr = fract * z2max;
                rbest = r;
            }
        }
    }
//  printf("r=%f, maxAngle=%f, fract=%f, n=%d\n", rbest, bestAngle, bestSnr, nn);
//  fflush(stdout);
    
    // (3a) set limits
    for (n=nn=0, sx2=0.; n<Parms->MaxVertex; n++)
        if (Mesh[n].z2 > bestSnr && Mesh[n].angle < bestAngle) {
            Snr[nn] = Mesh[n].z2;
            Rad[nn] = Mesh[n].angle;
            sx2 += Snr[nn] * Snr[nn];
            nn++;
        }
    Nused = nn;

    // (3b) output limits in Parms
    Parms->MaxAngle = bestAngle;
    Parms->MinSNR = bestSnr;

    // (4) test for Nused >= 12
    if (Nused >= 12) {

        // (4a) output rms z^2 values
        *rms = sqrt(sx2 / (double)Nused);
		gsl_fit_linear(Rad, 1, Snr, 1, Nused, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
//  printf("\t\tc0=%e, c1=%e, chisq=%e\n", c0, c1, chisq); fflush(stdout);
		*x2 = chisq / (Nused - 2);        // reduced chi^2
		*corr = -gsl_stats_spearman(Rad, 1, Snr, 1, Nused, work);
    } else {    // this should never happen!
        *rms = 0.;
		*x2 = -1.;
		*corr = HUGE;
    }
    *nv = Nused;

    free(Snr);
    free(Rad);
    free(work);
    ECDIntPntFree(coeff);
    free(coeff);
}
