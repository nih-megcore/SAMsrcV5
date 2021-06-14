// RandomMesh-- select vertices randomly from N input vertices.
//  Tuneable parameters are found in 'coreg.h'. The number of vertices sampled is
//  designated by 'MAX_VERTEX'.
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdlib.h>
#include <math.h>
#include <lfulib.h>
#include <mesh.h>
#include <voxel.h>
#include <coreg.h>

// pick MaxVertex vertices of the Vertex array and store them in Mesh
// |Vertex| must be >= |Mesh|

void RandomMesh(
    VOXELINFO   *Vertex,    // Vertex[N] -- raw vertex list
    PARMINFO    *Params,    // coregistration parameters
    MESHINFO    *Mesh,      // Mesh[MAX_VERTEX] -- output pointer to refined Mesh list for SAMcoreg
    int         N           // number of vertices
) {
    int i, j, n, v;
    int *idx;

    // randomize the order of the Vertex array
    idx = new_array(int, N);
    for (i = 0; i < N; i++) {
        idx[i] = i;
    }
    shuffle(idx, N);

    // choose the first MaxVertex of them to fill the Mesh array
    n = j = 0;
    while (n < Params->MaxVertex && j < N) {
        i = idx[j++];
        if (Vertex[i].ROI) {
            for (v=X_; v<=Z_; v++) {
                Mesh[n].Vp[v] = Vertex[i].p[v];
                Mesh[n].Vc[v] = Vertex[i].v[v];     // use anatomically-constrained vertex normal
            }
            n++;
        }
    }
    if (n < Params->MaxVertex) {
        fatalerr("n < MaxVertex in RandomMesh!");   // shouldn't happen
    }

    free(idx);
}
