// NSolveFwd() - solve for beamformer optimal forward solution
//	for any N current dipole in realistic head model. Call
//  NSolveMoments() to set moment vectors.
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <lfulib.h>
#include <math.h>
#include <voxel.h>
#include <siglib.h>

void    NSolveFwd(
	VOXELINFO	*Voxel,			// Voxel[N] -- (source) voxel structure -- normal vector req'd for solving orientation (input)
	gsl_matrix  *B,             // B - MxN forward solution matrix (externally allocated)
	N_LEADFIELD *Lfld,          // Lfld[N] - array of N-source Mx3 lead-field matrices
	gsl_matrix	*C,				// C - MxM orientation covariance matrix
	int			N				// number of sources
)
{
    gsl_matrix          *Vec;   // Vec - Nx3 moment orientation vectors
    gsl_matrix          *Cinv;  // Cinv - MxM inverse orientation covariance matrix
	double              tmp;    // what else?
	int                 solve;	// TRUE | FALSE -- flag for whether to solve for orientation
	int					m;		// sensor index
	int                 n;      // target index
	int					v;		// vector index
    int                 M;      // number of primary sensors

	// allocate memory & invert covariance
	M = C->size1;
    Cinv = gsl_matrix_alloc(M, M);
    pinv(C, Cinv);
    Vec = gsl_matrix_alloc(N, 3);

    // test voxel flags as to whether to solve dipole orientation
    for(n=0, solve=TRUE; n<N; n++)
         if(Voxel[n].Solve == FALSE) {   // if any one flag is false, don't solve for orientation
            solve = FALSE;
            break;
		}

    // if req'd, solve for N moment vectors & copy to voxel/targets
    if(solve == TRUE) {                 // linear solutions for N moment vectors
        NSolveMoments(Vec, Cinv, Lfld, N, &tmp);
        for(n=0; n<N; n++)
            for(v=X_; v<=Z_; v++)
                Voxel[n].v[v] = gsl_matrix_get(Vec, n, v);
    }

    // compute forward solutions for N sources
    for(n=0; n<N; n++)
        for(m=0; m<M; m++) {
            for(v=X_, tmp=0.; v<=Z_; v++)
                tmp += Voxel[n].v[v] * gsl_matrix_get(Lfld[n].L, m, v);
            gsl_matrix_set(B, m, n, tmp);
        }

    // free allocations
    gsl_matrix_free(Cinv);
    gsl_matrix_free(Vec);
}
