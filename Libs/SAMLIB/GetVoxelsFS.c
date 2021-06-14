// GetVoxelsFS() -- populate an array of voxels using a
//  FreeSurfer atlas (from FSnormals.py)
//
//      Author: Stephen E. Robinson
//              MEG Core Facility
//              NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <geoms.h>
#include <voxel.h>
#include <lfulib.h>
#include <samlib.h>
#include "samutil.h"

#define MAX_RAD     0.12        // 120 mm maximum voxel distance from local sphere origin
#define MIN_RAD     0.005       // 5 mm minimum voxel distance from local sphere origin
#define ETA         1.0e-12     // a small fraction of a step

/* detect a string of decimal digits */

static int isint(char *s)
{
    do {
        if (!isdigit(*s)) {
            return FALSE;
        }
    } while (*++s);
    return TRUE;
}

void GetVoxelsFS(
    char        *AtlasName,
    VOXELINFO   **OutVoxel,
    int         *NumVoxels,
    double      Vmin[3],
    double      Vmax[3],
    double      *stepsize)
{
    VOXELINFO           *Voxel;         // array of voxels
    double              delta[3];       // step size
    double              Vc[3];          // position vector
    double              Vn[3];          // normal vector
    double              Vlast[3];       // last voxel
    double              r;              // normal length
    double              r2;             // r^2
    double              tmp;            // you know what this means, don't you?
    int                 V;              // number of voxels
    int                 Vstep[3];       // voxel steps
    int                 N;              // number of vertices
    int                 i;              // array index
    int                 j;              // j-index
    int                 k;              // k-index
    int                 n;              // voxel index
    int                 v;              // vector index
    int                 flag;           // voxel found flag
    int                 nflag;          // file begins with integer
    char                buf[256];       // line buffer
    FILE                *fp;            // file pointer for headshape

    // read Atlas file
    fp = fileopen(AtlasName, "r");

    // first read -- determine N, ROI, & step size
    N = 0;
    for (v=X_; v<=Z_; v++) {
        Vmax[v] = 0.;
        Vmin[v] = HUGE;
        Vlast[v] = HUGE;
        delta[v] = HUGE;
    }

    /* For backwards compatibility, if the first line is a
    single integer, ignore it. */

    nflag = FALSE;
    fgetline(buf, sizeof(buf), fp);
    if (isint(buf)) {
        nflag = TRUE;
        fgetline(buf, sizeof(buf), fp);
    }

    do {
        N++;

        // parse one voxel, position and normal
        if (sscanf(buf, "%lf%lf%lf%lf%lf%lf", &Vc[X_], &Vc[Y_], &Vc[Z_], &Vn[X_], &Vn[Y_], &Vn[Z_]) != 6)
            cleanup("error reading vertex");

        // find ROI maxima
        for (v=X_; v<=Z_; v++)
            if (Vc[v] > Vmax[v])
                Vmax[v] = Vc[v];

        // find ROI minima
        for (v=X_; v<=Z_; v++)
            if (Vc[v] < Vmin[v])
                Vmin[v] = Vc[v];

        // find smallest non-zero separation (step size) in each direction
        for (v=X_; v<=Z_; v++) {
            tmp = fabs(Vc[v] - Vlast[v]);
            if (tmp < delta[v] && tmp > ETA)
                delta[v] = tmp;
        }
        for (v=X_; v<=Z_; v++)
            Vlast[v] = Vc[v];

    } while (fgetline(buf, sizeof(buf), fp) != NULL);

#if DEBUG
    printf("\nread %d voxels from '%s'\n", N, AtlasName);
    printf("x: start=%f, end=%f (cm)\n", 100.*Vmin[X_], 100.*Vmax[X_]);
    printf("y: start=%f, end=%f (cm)\n", 100.*Vmin[Y_], 100.*Vmax[Y_]);
    printf("z: start=%f, end=%f (cm)\n", 100.*Vmin[Z_], 100.*Vmax[Z_]);
    printf("steps = %f, %f, %f\n", 100.*delta[X_], 100.*delta[Y_], 100.*delta[Z_]);
#endif

    // test for equal steps!
    for (v=X_; v<=Z_; v++)
        Vstep[v] = (int)rint(10000.*delta[v]);  // convert delta to nearest integer 0.1 mm
    if (Vstep[X_] != Vstep[Y_] && Vstep[X_] != Vstep[Z_])
        cleanup("unequal voxel step size!");
    *stepsize = (double)Vstep[X_] / 10000.;

    // determine number of voxels in ROI (this is greater than N!)
    for (r=Vmin[X_], Vstep[X_]=0; r<=(Vmax[X_]+ETA); r+=delta[X_], Vstep[X_]++);
    for (r=Vmin[Y_], Vstep[Y_]=0; r<=(Vmax[Y_]+ETA); r+=delta[Y_], Vstep[Y_]++);
    for (r=Vmin[Z_], Vstep[Z_]=0; r<=(Vmax[Z_]+ETA); r+=delta[Z_], Vstep[Z_]++);
    V = Vstep[X_] * Vstep[Y_] * Vstep[Z_];

#if DEBUG
    printf("allocating %d voxel array\n", V);
#endif
    Voxel = new_array(VOXELINFO, V);

    // populate voxel grid with 3D coordinates
    for (n=i=0, Vc[X_]=Vmin[X_]; i<Vstep[X_]; i++, Vc[X_]+=delta[X_])
        for (j=0, Vc[Y_]=Vmin[Y_]; j<Vstep[Y_]; j++, Vc[Y_]+=delta[Y_])
            for (k=0, Vc[Z_]=Vmin[Z_]; k<Vstep[Z_]; k++, Vc[Z_]+=delta[Z_]) {
                for (v=X_; v<=Z_; v++) {
                    Voxel[n].p[v] = Vc[v];
                    Voxel[n].v[v] = 0.;
                }
                Voxel[n].Solve = FALSE;
                Voxel[n].ROI = FALSE;
                n++;
            }

    assert(V == n);

    // re-read atlas list & put vertex normals in correct indices
    rewind(fp);

    // for each voxel...

    if (nflag) {
        fgetline(buf, sizeof(buf), fp);
    }
    while (fgetline(buf, sizeof(buf), fp) != NULL) {

        // parse one voxel
        if (sscanf(buf, "%lf%lf%lf%lf%lf%lf", &Vc[X_], &Vc[Y_], &Vc[Z_], &Vn[X_], &Vn[Y_], &Vn[Z_]) != 6)
            cleanup("error reading vertex");

        // enforce unit normal vector
        for (v=X_, r2=0.; v<=Z_; v++)
            r2 += Vn[v] * Vn[v];
        if (r2 == 0.)
            cleanup("found zero-length normal vector in atlas");
        r = sqrt(r2);
        for (v=X_; v<=Z_; v++)
            Vn[v] = Vn[v] / r;

        // determine index for this coordinate
        for (i=0, flag=FALSE; i<V; i++) {

            // compute distance from ith voxel to this coordinate
            for (v=X_, r2=0.; v<=Z_; v++) {
                tmp = Vc[v] - Voxel[i].p[v];
                r2 += tmp * tmp;
            }
            r = sqrt(r2);

            // is this our voxel?
            if (r < ETA) {
                for (v=X_; v<=Z_; v++)
                    Voxel[i].v[v] = Vn[v];
                Voxel[i].ROI = TRUE;
#if DEBUG
    printf("found index %d\n", i);
#endif
                flag = TRUE;
                break;
            }
        }
        if (!flag)
            cleanup("couldn't insert atlas voxel coordinate");
    }
    fclose(fp);

    // output
    *OutVoxel = Voxel;
    *NumVoxels = V;
}
