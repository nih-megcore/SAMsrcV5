// GetAtlas() -- read cortical atlas of vertices & normals
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <geoms.h>
#include <voxel.h>
#include "samutil.h"

void GetAtlas(
    char        *PathName,      // full path of named atlas (output of FSnormals.py)
    VOXELINFO   **VertexList,   // list of vertices & their unit normals, with flags
    int         *N              // number of vertices
) {
    VOXELINFO   *Vertex;        // list of vertices with normals & flags
    double      Vc[3];          // vertex position
    double      Vr[3];          // vertex normal
    double      r;              // distance
    double      r2;             // r^2
    int         code;           // roi code
    int         A;              // total atlas lines
    int         i;              // vertex index
    int         v;              // vector index
    char        Line[256];      // line buffer
    FILE        *fp;

    // open atlas file & count the lines
    fp = fileopen(PathName, "r");
    A = 0;
    while (fgets(Line, sizeof(Line), fp) != NULL)   // count lines
        A++;
    if (A == 0)
        cleanup("empty atlas");
    rewind(fp);

    // allocate vertices
    Vertex = new_array(VOXELINFO, A);
    *VertexList = Vertex;

    // read lines while ignoring those beginning with non-numeric characters
    i = 0;
    while (fgets(Line, sizeof(Line), fp) != NULL) {
        if (isalpha(Line[0]) == 0) {
            if (sscanf(Line, "%lf%lf%lf%lf%lf%lf%d", &Vc[X_], &Vc[Y_], &Vc[Z_], &Vr[X_], &Vr[Y_], &Vr[Z_], &code) != 7)
                cleanup("vertex read failed");
            for (v=X_, r2=0.; v<=Z_; v++)
                r2 += Vr[v] * Vr[v];
            r = sqrt(r2);
            for (v=X_; v<=Z_; v++) {
                Vertex[i].p[v] = Vc[v];
                Vertex[i].v[v] = Vr[v] / r;
            }
            Vertex[i].ROI = FALSE;
            i++;
        }
    }
    *N = i;
    fclose(fp);
}
