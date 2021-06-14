// GetVoxels() -- populate an array of voxels for the define ROI
//      & set limits according to whether or not a hull model is present
//
//      Author: Stephen E. Robinson
//              MEG Core Facility
//              NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <geoms.h>
#include <voxel.h>
#include <lfulib.h>
#include <samlib.h>

#define MAX_RAD		0.12			// 120 mm maximum voxel distance from local sphere origin
#define MIN_RAD		0.005			// 5 mm minimum voxel distance from local sphere origin
#define ETA			1.0e-12			// a small fraction of a step

void GetVoxels(
    VOXELINFO	**OutVoxel,
    int			*NumVoxels,
    HULL		*Hull,				// if Hull->nv == 0, no hull present
    double		Ts[3],
    double		Te[3],
    double		step)
{
	static VOXELINFO	*Voxel;		// array of voxels
	double				Vc[3];		// cartesian vector
	double				u1[3];		// unit vector from voxel to vertex
	double				r;			// radius
	double				r2;			// radius^2
	double				u2;			// vertex normal vector^2
	double				dot;		// dot product of u1 & u2
	double				min;		// minimum distance
	double				tmp;		// the variable who's real name we don't dare speak!
	int					V;			// number of voxels
	int					c;			// voxel index
	int					h;			// hull vertex index
	int					index = 0;	// nearest neigbor index
	int					i;
	int					j;
	int					k;
	int					v;			// vector index

	// count number of voxels in ROI
	for(Vc[X_]=Ts[X_], i=0; Vc[X_]<=(Te[X_]+ETA); Vc[X_]+=step, i++);
	for(Vc[Y_]=Ts[Y_], j=0; Vc[Y_]<=(Te[Y_]+ETA); Vc[Y_]+=step, j++);
	for(Vc[Z_]=Ts[Z_], k=0; Vc[Z_]<=(Te[Z_]+ETA); Vc[Z_]+=step, k++);
	V = i * j * k;

	// allocate array of voxels & set coordinates
	if((Voxel = (VOXELINFO *)malloc((size_t)V * sizeof(VOXELINFO))) == NULL)
		cleanup("'Voxel[]' allocation failed");

	// set voxel coordinates with ROI flag=FALSE
	for(Vc[X_]=Ts[X_], c=0; Vc[X_]<=(Te[X_]+ETA); Vc[X_]+=step)
		for(Vc[Y_]=Ts[Y_]; Vc[Y_]<=(Te[Y_]+ETA); Vc[Y_]+=step)
			for(Vc[Z_]=Ts[Z_]; Vc[Z_]<=(Te[Z_]+ETA); Vc[Z_]+=step, c++) {
				for(v=X_; v<=Z_; v++)
					Voxel[c].p[v] = Vc[v];
				Voxel[c].Solve = TRUE;
				Voxel[c].ROI = FALSE;
			}

	// set voxel flags
	for(c=0; c<V; c++) {

		if(Hull->nv != 1) {		// if Hull present, test if voxel is inside or outside the hull

			// find hull vertex nearest to this voxel -- with non-zero normal!!!
			for(h=0, min=HUGE; h<Hull->nv; h++) {
				for(v=X_, r2=u2=0.; v<=Z_; v++) {
					tmp = Voxel[c].p[v] - Hull->vertex[h].p[v];
					r2 += tmp * tmp;
					u2 += Hull->vertex[h].v[v] * Hull->vertex[h].v[v];
				}
				if(r2 < min && u2 > 0.) {
					min = r2;
					index = h;
				}
			}

			// compute unit vector for voxel to vertex
			for(v=X_, r2=0.; v<=Z_; v++) {
				u1[v] = Voxel[c].p[v] - Hull->vertex[index].p[v];
				r2 += u1[v] * u1[v];
			}
			r = sqrt(r2);
			for(v=X_; v<=Z_; v++)
				u1[v] /= r;

			// compute dot u1 dot u2
			for(v=X_, dot=0.; v<=Z_; v++)
				dot += u1[v] * Hull->vertex[index].v[v];
			if(dot < 0.) {
				Voxel[c].Solve = TRUE;
				Voxel[c].ROI = TRUE;
			} else {
				Voxel[c].Solve = FALSE;
				Voxel[c].ROI = FALSE;
			}

		} else {						// there is no hull, but first element has local sphere origin

			// compute voxel radial distance from local sphere origin
			for(v=X_, r2=0.; v<=Z_; v++) {
				tmp = Voxel[c].p[v] - Hull->Vo[v];
				r2 += tmp * tmp;
			}
			r = sqrt(r2);
			if(r > MIN_RAD && r < MAX_RAD) {
				Voxel[c].Solve = TRUE;
				Voxel[c].ROI = TRUE;
			} else {
				Voxel[c].Solve = FALSE;
				Voxel[c].ROI = FALSE;
			}
		}
	}
	*NumVoxels = V;
	*OutVoxel = Voxel;
}
