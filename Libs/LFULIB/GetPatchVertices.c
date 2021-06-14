// GetPatchVertices -- routine to reduce number of vertices for patch solution

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <model.h>
#include <samlib.h>
#include <lfulib.h>
#include <voxel.h>

#define PATCH_SIZE	10000				// sufficient number of vertices for 'VertexList' patch
#define MIN_VERTEX	10					// minimum number of vertices for solution


void	GetPatchVertices (
	VOXELINFO	*Vertex,				// Vertex[N] -- external array				                (input)
	int			NumVert,				// number of vertices						                (input)
	double		*Vcentroid,				// Vcentroid[3] -- patch origin				                (input)
	double		Extent,					// radial extent from patch origin			                (input)
	int			*NumUsed,				// number of vertices used					                (output)
	int         Solve                   // TRUE: solve moment vector; FALSE: use cortical normal    (input)
)
{
	static VOXELINFO	*VertexList;	// intermediate list of all vertices within target ROI
	double				r2;				// r^2
	double				r;				// distance
	double				tmp;			// an oft-used variable name
	int					Nroi;			// number of ROI vertices
	int					i;				// all-purpose index
	int					n;				// vertex index
	int					v;				// vector index
	static int			Npatch = 0;

	// allocate 'VertexList' sufficient to accomodate the patch of your choice
	if(Npatch == 0) {
		if((VertexList = (VOXELINFO *)malloc(PATCH_SIZE * sizeof(VOXELINFO))) == NULL)
			allocfailed("VertexList[]");
		Npatch = PATCH_SIZE;
	}

	// generate 'VertexList' array for _all_ vertices within target ROI
	for(n=i=0; n<NumVert; n++) {

		// compute distance between ROI centroid & vertex
		for(v=X_, r2=0.; v<=Z_; v++) {
			tmp = Vcentroid[v] - Vertex[n].p[v];
			r2 += tmp * tmp;
		}
		r = sqrt(r2);

		// populate list of vertices within 'Extent'
		if(r <= Extent) {
			for(v=X_; v<=Z_; v++) {
				VertexList[i].p[v] = Vertex[n].p[v];
				VertexList[i].v[v] = Vertex[n].v[v];
			}
			i++;
		}
	}
	Nroi = i;

	// sanity check
	if(Nroi > PATCH_SIZE)
		cleanup("too many vertices -- decrease 'Extent'");
	if(Nroi >= MIN_VERTEX) {		// sufficient number of vertices for this voxel/target?
		for(n=0; n<Nroi; n++)		// copy vertices marked true to 'VertexUsed' array
			for(v=X_; v<=Z_; v++) {
				Vertex[n].p[v] = VertexList[n].p[v];
				Vertex[n].v[v] = VertexList[n].v[v];
			}
			Vertex[i].ROI = TRUE;
			Vertex[i].Solve = Solve;
		*NumUsed = Nroi;
	} else {
		*NumUsed = 0;
	}
}
