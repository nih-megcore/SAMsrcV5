// NLeadField() - compute the M x 3N lead field matrix
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <stdlib.h>
#include <lfulib.h>
#include "sam_param.h"

void	NLeadField(
	VOXELINFO	*Voxel,			// Voxel[N] -- array of N voxel structures (input)
	N_LEADFIELD *Lfld,          // Lfld[N] -- array of Mx3 lead-field matrices for N sources (output)
	HeaderInfo	*Header,		// header information with sensor list (input)
	ChannelInfo	*Channel,		// sensor information (input)
	HULL		*Hull,			// cortical hull model (input)
	int			N,				// number of sources (input)
	int         Model           // SSPHERE | MSPHERE | NOLTE (input)
)

{
    gsl_matrix  *Ls;            // Ls - Mx3 single source lead-field matrix
    gsl_vector  *xyz;           // single sensor lead-field vector
    double      Vc[3];          // cartesian position vector
    int         m;              // sensor index
    int         n;              // target index
	int         v;              // vector index
	int         M;              // number of sensors

    // get dimensions of Lfld
    M = Lfld[0].L->size1;
    if(Lfld[0].L->size2 != 3)
        cleanup("column dimension of Lfld != 3");

    // allocate gsl matrix
    Ls = gsl_matrix_alloc(M, 3);
    xyz = gsl_vector_alloc(3);

	// compute ECD lead field matrix - compute forward solution matrix for unit dipole terms
	for(n=0; n<N; n++) {
		for(v=X_; v<=Z_; v++)
			Vc[v] = Voxel[n].p[v];
        switch(Model) {
            case SSPHERE:
            case MSPHERE:
                DIPLeadField(Vc, Ls, Header, Channel, 2);   // note: 2nd-order integration in sensor
                break;
            case NOLTE:
#if ELEKTA
                ECDLeadFieldElekta(Vc, Ls, Header, Channel, Hull);
#else
                ECDLeadField(Vc, Ls, Header, Channel, Hull);
#endif
                break;
            default:
                break;
        }

        // save Ls in Lfld
        for(m=0; m<M; m++) {
            gsl_matrix_get_row(xyz, Ls, m);
            gsl_matrix_set_row(Lfld[n].L, m, xyz);
        }
	}

    // free allocations
    gsl_matrix_free(Ls);
    gsl_vector_free(xyz);
}
