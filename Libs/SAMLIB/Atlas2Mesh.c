// Atlas2Mesh() -- create mesh from atlas, using MaxDot, MinZ, & origin
//  Mesh must be allocated using the number of vertices in the entire atlas
//
//	Author:
//		Stephen E. Robinson, PhD
//		MEG Core Group
//		NIMH
//

#include <stdlib.h>
#include <math.h>
#include <lfulib.h>
#include <mesh.h>
#include <voxel.h>


void    Atlas2Mesh(
    VOXELINFO   *Atlas,     // Atlas[Natlas] - cortical mesh vertices and normals (input)
    int         Natlas,     // number of Atlas structs (input)
    MESHINFO    *Mesh,      // Mesh[Nmesh] -- cortical mesh for solution info (output)
    PARMINFO    *Params     // parameter structure
)
{
    double  Vo[3];          // origin (centroid) of atlas
    double  Vr[3];          // radial vector from origin to vertex
    double  r;              // length of radial vector
    double  r2;             // r^2
    double  dot;            // dot product
    int     n;              // atlas index
    int     nn;             // mesh index
    int     v;              // vector index

    // compute centroid of atlas
    Vo[X_] = Vo[Y_] = Vo[Z_] = 0.;
    for (n=0; n<Natlas; n++)
        for (v=X_; v<=Z_; v++)
            Vo[v] += Atlas[n].p[v];
    for (v=X_; v<=Z_; v++)
        Vo[v] /= (double)Natlas;

    // for each vertex in the atlas...
    for(n=nn=0; n<Natlas; n++) {

        // compute unit vector from origin to vertex
        for(v=X_, r2=0.; v<=Z_; v++) {
            Vr[v] = Atlas[n].p[v] - Vo[v];  // radial vector from origin to vertex
            r2 += Vr[v] * Vr[v];
        }
        r = sqrt(r2);
        for(v=X_; v<=Z_; v++)                // normalize
            Vr[v] /= r;

        // compute dot product of Vr & atlas normal vector
        for(v=X_, dot=0.; v<=Z_; v++)
            dot += Atlas[n].v[v] * Vr[v];   // dot product with atlas normal
        if(fabs(dot) < Params->MaxDot && Atlas[n].p[Z_] > Params->MinZ) {
            for(v=X_; v<=Z_; v++) {
                Mesh[nn].Vp[v] = Atlas[n].p[v];
                Mesh[nn].Vc[v] = Atlas[n].v[v];
            }
            nn++;
        }
    }
    Params->MaxVertex = nn;      // MaxVertex has the actual number of vertices
}
