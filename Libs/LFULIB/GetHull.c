// GetHull() -- read cortical hull model from text file in dataset directory
//      note that vertices are indexed base 1, for now...
//
//      Author: Stephen E. Robinson
//              MEG Core Facility
//              NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lfulib.h>

void GetHull(
    char    *hullName,          // path of hull.shape file
    HULL    *Hull               // cortical hull model with position & normals of vertices, including origin & scale
) {
    double  Vp[3];              // cartesian position vector
    double  Vo[3];              // cartesian orientation vector
    double  r;                  // normal vector magnitude
    double  r2;                 // r^2
    double  sum;                // simple sum
    double  tmp;                // what else?
    int     i;
    int     j;
    int     v;
    FILE    *fp;

    // read number of vertices & allocate memory
    fp = fileopen(hullName, "r");
    if (fscanf(fp, "%d", &Hull->nv) != 1)
        cleanup("can't read number of vertices");
    Hull->vertex = new_arrayE(FIELD, Hull->nv, "hull vertex[]");

    // read vertices & normals
    for (i=j=0; i<Hull->nv; i++) {
        if (fscanf(fp, "%lf%lf%lf", &Vp[X_], &Vp[Y_], &Vp[Z_]) != 3)
            cleanup("can't read vertex position");
        for (v=X_; v<=Z_; v++)
            if (Vp[v] > 0.2)                    // if position > 20 cm, something is very wrong!
                cleanup("vertex > 20 cm from origin -- corrupted 'hull.shape'?");
        if (fscanf(fp, "%lf%lf%lf", &Vo[X_], &Vo[Y_], &Vo[Z_]) != 3)
            cleanup("can't read vertex orientation");
        for (v=X_, r2=0.; v<=Z_; v++)           // enforce unit normal vector
            r2 += Vo[v] * Vo[v];
        if (r2 > 0.) {
            r = sqrt(r2);
            for (v=X_; v<=Z_; v++) {
                Hull->vertex[i].p[v] = Vp[v];   // hull.shape position is already metres
                Hull->vertex[i].v[v] = Vo[v] / r;
            }
            j++;                                // count number of vertices with non-zero normals
        }
    }
    fclose(fp);

    // test good vertices
    if (j < (Hull->nv-100)) // missing more than 100 vertices!
        cleanup("more than 100 zero-length vertices -- corrupted 'hull.shape'?");
    if (j != Hull->nv) {
        fprintf(stderr, " - note: %d hull vertices have zero-length normals", Hull->nv - j);
        Hull->nv = j;       // revise number of vertices with non-zero normals
    }

    // compute hull geometric center
    for (i=0, Vp[X_]=Vp[Y_]=Vp[Z_]=0.; i<Hull->nv; i++)
        for (v=X_; v<=Z_; v++)
            Vp[v] += Hull->vertex[i].p[v];
    for (v=X_; v<=Z_; v++)
        Hull->Vo[v] = Vp[v] / (double)Hull->nv;

    // compute hull scale
    for (i=0, sum=0.; i<Hull->nv; i++) {
        for (v=X_, r2=0.; v<=Z_; v++) {
            tmp = Hull->vertex[i].p[v] - Hull->Vo[v];
            r2 += tmp * tmp;
        }
        sum += sqrt(r2);
    }
    Hull->scale = sum / (double)Hull->nv;

    // make inflated hull for Nolte solution

    Hull->infvertex = new_arrayE(FIELD, Hull->nv, "Hull->infvertex[]");
    for (v = 0; v < Hull->nv; v++) {
        for (i = X_; i <= Z_; i++) {
            // assumes that normals are all outward-pointing
            Hull->infvertex[v].p[i] = Hull->vertex[v].p[i] + INFLATE * Hull->vertex[v].v[i];
            Hull->infvertex[v].v[i] = Hull->vertex[v].v[i];
        }
    }
}
