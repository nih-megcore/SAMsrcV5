// ROISolveWts -- compute beamformer weight for extended source defined by cortical
//	patch. Simple algorithm reconstructs a low-order forward solution for a patch of
//	vertices, using an SVD of the MxN forward solutions for each vertex.
//
//			B[i] = Sum W[k] * U[i][k] * V[k][j]
//
//	where W[k] is set to zero at the rank limit for MAX_FRACT. This version
//	uses either multiple local spheres (if Hull == NULL) or the realistic head model
//  (Nolte) forward solution! For the latter, initial startup is slow due to
//  computation of the spherical harmonics for the cortical hull.
//
//	Author: Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <fcntl.h>
#include <model.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <DataFiles.h>
#include <SAMfiles.h>

#define MAX_FRACT	0.95


void    ROISolveWts(
    VOXELINFO   *ROI,           // list of ROI vertices
    int         NumVertices,    // number of ROI vertices (N)
    gsl_vector  *W,             // W - Mx1 beamformer coefficients
    gsl_matrix  *Cs,            // Cs - MxM solution covariance matrix
    gsl_matrix  *Co,            // Co - MxM orientation covariance matrix
    HeaderInfo  *Header,        // MEG header structure
    ChannelInfo *Channel,       // MEG channel structure
    HULL        *Hull,          // cortical hull (if Hull.nv == 1, then use multisphere)
    int         *Rank,          // rank -- number of bases
    int         Solve           // solve moment vector -- TRUE | FALSE
)

{
	gsl_matrix			*A;			// A SVD matrix
	gsl_matrix			*V;			// V SVD matrix
	gsl_vector			*S;			// S SVD vector -- singular values
	gsl_matrix			*X;			// X SVD matrix
	gsl_vector			*work;		// SVD workspace
	gsl_matrix			*Sz;		// zero padded S
	gsl_matrix			*K;			// intermediate matrix
	gsl_matrix			*L;			// L - Mx3 primary sensor lead field
	static gsl_matrix	*Csinv;		// Csinv - MxM inverse covariance matrix
	static gsl_matrix	*Coinv;		// Coinv - MxM inverse orientation covariance matrix
	static gsl_vector	*B;			// rank-reduced approximation to forward solution
	double				Vu[3];		// unconstrained source vector
	double				gain;		// gain (B' C^-1 B)
	double				sum;		// sum of singular values
	double				ratio;		// ratio of running sum of singular values / sum
	double				span;		// span of source moment vector solution
	double				tmp;		// the ubiquitous variable
	int					M;			// number rows
	int					N;			// number columns
	int					i;			// i-index
	int					j;			// j-index
	int					k;			// k-index
	int					v;			// vector index
	int					R;			// rank of bases
	int					transpose;	// transposition flag
	static int			NumSensors = 0;	// saved number of primary sensors (M)

	// memory allocation of vectors & matrices that depend on NumSensors
	if(NumSensors == 0) {
		NumSensors = Cs->size1;
		Csinv = gsl_matrix_alloc(NumSensors, NumSensors);
		pinv(Cs, Csinv);
		if(Solve == TRUE) {
			Coinv = gsl_matrix_alloc(NumSensors, NumSensors);
			pinv(Co, Coinv);
		}
		B = gsl_vector_alloc(NumSensors);
	}

	// initialize dimensions & rank
	if(NumVertices > NumSensors) {
		M = NumVertices;
		N = NumSensors;
		transpose = FALSE;
	} else {
		M = NumSensors;
		N = NumVertices;
		transpose = TRUE;
	}
	R = 0;

	// allocate gsl matrices
	L = gsl_matrix_alloc(NumSensors, 3);
	A = gsl_matrix_alloc(M, N);
	S = gsl_vector_alloc(N);
	V = gsl_matrix_alloc(N, N);
	X = gsl_matrix_alloc(N, N);
	work = gsl_vector_alloc(N);

	// compute H for all ROI vertices
	for(i=0; i<NumVertices; i++) {				// loop on NumVertices

		// compute lead-field
		if (Hull->nv == 1)
		    DIPLeadField(ROI[i].p, L, Header, Channel, 3);      // multisphere with 3rd-order coil integration
		else
		    ECDLeadField(ROI[i].p, L, Header, Channel, Hull);   // Nolte solution

		// if req'd, solve unconstrained moment vector
		if(Solve == TRUE) {
			SolveMoment(Coinv, L, &span, Vu);
			for(j=0; j<NumSensors; j++)
				for(v=X_; v<=Z_; v++)
					ROI[i].v[v] = Vu[v];		// replace vertex normal vector with moment solution
		}

		// compute forward solution
		for(j=0; j<NumSensors; j++) {
			for(v=X_, tmp=0.; v<=Z_; v++)
				tmp += ROI[i].v[v] * gsl_matrix_get(L, j, v);
			if(transpose == FALSE)
				gsl_matrix_set(A, i, j, tmp);
			else
				gsl_matrix_set(A, j, i, tmp);
		}
	}

	// compute SVD of A -- note that V is returned, not its transpose!
	if(gsl_linalg_SV_decomp_mod(A, X, V, S, work) == -1)
		cleanup("gsl_linalg_SV_decomp_mod() failed");
	gsl_matrix_free(X);
	gsl_vector_free(work);
	gsl_matrix_free(L);

	// find rank for max_fract threshold -- assumes SVD orders singular values from largest to smallest
    Sz = gsl_matrix_calloc(N, N);           // allocate Sz with zero elements
    for(i=0, sum=0.; i<N; i++)              // sum all singular values
        sum += gsl_vector_get(S, i);
    for(i=R=0, tmp=0.; i<N; i++) {          // keep running sum of singular values until threshold is met
        tmp += gsl_vector_get(S, i);
        ratio = tmp / sum;
        if(ratio < MAX_FRACT || i == 0) {   // solution must be at least rank 1!
            gsl_matrix_set(Sz, i, i, gsl_vector_get(S, i));
            R++;
		} else {
			gsl_matrix_set(Sz, i, i, 0.);		// zero remaining singular values
		}
	}
	*Rank = R;
	gsl_vector_free(S);

	// reconstruct rank-reduced approximation to forward solution
	K = gsl_matrix_calloc(N, N);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Sz, 0., K);
	gsl_matrix_free(Sz);
	gsl_matrix_free(V);
	for(i=0; i<NumSensors; i++)
		for(j=0, tmp=0.; j<R; j++) {
			for(k=0; k<R; k++)
				tmp += gsl_matrix_get(A, i, k) * gsl_matrix_get(K, k, j);
			gsl_vector_set(B, i, tmp);
		}
	gsl_matrix_free(A);
	gsl_matrix_free(K);

	// compute beamformer coefficients: W = B' C^-1 B
	gsl_blas_dgemv(CblasNoTrans, 1., Csinv, B, 0., W);		// W = C^(-1) B
	gsl_blas_ddot(B, W, &tmp);
	gain = 1. / tmp;
	gsl_vector_scale(W, gain);
}
