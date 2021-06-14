// Atlas2RandomMesh-- select vertices randomly from N input vertices.
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

void Atlas2RandomMesh(
    VOXELINFO   *Atlas,     // Atlas[Natlas] - cortical mesh vertices and normals (input)
    int         Natlas,     // number of Atlas structs (input)
    MESHINFO    *Mesh,      // Mesh[Nmesh] -- cortical mesh for solution info (output)
    PARMINFO    *Params     // parameter structure
) {
    double  Vo[3];          // origin (centroid) of atlas
    double  Vr[3];          // radial vector from origin to vertex
    double  r;              // length of radial vector
    double  r2;             // r^2
    double  dot;            // dot product
    int     i;
    int     j;
    int     n;              // atlas index
    int     v;              // vector index
    int     *idx;           // array of random vertices

    // compute centroid of atlas
    Vo[X_] = Vo[Y_] = Vo[Z_] = 0.;
    for (n=0; n<Natlas; n++)
        for (v=X_; v<=Z_; v++)
            Vo[v] += Atlas[n].p[v];
    for (v=X_; v<=Z_; v++)
        Vo[v] /= (double)Natlas;

    // randomize the order of the Vertex array
    idx = new_array(int, Natlas);
    for (n=0; n<Natlas; n++)
        idx[n] = n;
    shuffle(idx, Natlas);

    // choose the first MaxVertex of them to fill the Mesh array
    n = j = 0;
    while (n < Params->MaxVertex && j < Natlas) {
        i = idx[j++];

        // compute unit vector from origin to vertex
        for(v=X_, r2=0.; v<=Z_; v++) {
            Vr[v] = Atlas[i].p[v] - Vo[v];  // radial vector from origin to vertex
            r2 += Vr[v] * Vr[v];
        }
        r = sqrt(r2);
        for(v=X_; v<=Z_; v++)               // normalize
            Vr[v] /= r;

        // compute dot product of Vr & atlas normal vector
        for(v=X_, dot=0.; v<=Z_; v++)
            dot += Atlas[i].v[v] * Vr[v];   // dot product with atlas normal
        if(fabs(dot) < Params->MaxDot && Atlas[n].p[Z_] > Params->MinZ) {
            for(v=X_; v<=Z_; v++) {
                Mesh[n].Vp[v] = Atlas[i].p[v];
                Mesh[n].Vc[v] = Atlas[i].v[v];
            }
            n++;
        }
    }
    if (n < Params->MaxVertex) {
        fatalerr("n < MaxVertex in Atlas2RandomMesh!"); // shouldn't happen
    }

    free(idx);
}
