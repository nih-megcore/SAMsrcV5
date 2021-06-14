// GetHDM() -- reads *.hdm file into HEAD_MODEL structure
//
//  Author: Stephen E. Robinson
//          Neuromagnetism Laboratory
//          Dept.of Neurology
//          Henry Ford Hospital
//          Detroit, MI
//          Copyright (c) 2007
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>      // get axes definitions
#include <samlib.h>
#include <lfulib.h>


void    GetHDM(
    char        *HdmName,
    HEAD_MODEL  *Model
) {
    double  x;
    double  y;
    double  z;
    double  r;
    int     M;
    int     v;
    int     count;
    char    Label[256];
    char    Line[256];
    char    Name[256];
    char    *token;
    FILE    *fp;

    // open head model file
    if((fp = fopen(HdmName, "r")) == NULL)
        cleanup("'default.hdm' open failed");

    // read mri file (path)
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'MRI_File:' label");
    } while(strncmp(Label, "MRI_FILE:", 9) != 0);
    if(fscanf(fp, "%s", Line) != 1)
        cleanup("can't read MRI file name");

    // parse mri file name from full path string
    token = strtok(Line, "/");
    for(;;) {
        token = strtok(NULL, "/");
        if(token != NULL)
            strcpy(Model->MriName, token);
        else
            break;
    }

    // read single (or mean) local sphere values
    do {                                                    // synch up to line containing 'ORIGIN_X'
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'ORIGIN_X' label");
    } while(strncmp(Label, "ORIGIN_X:", 9) != 0);
    if(fscanf(fp, "%lf", &Model->SphereOrigin[X_]) != 1)    // read sphere x-origin
        cleanup("can't read sphere x-origin");
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'ORIGIN_Y' label");
    } while(strncmp(Label, "ORIGIN_Y:", 9) != 0);
    if(fscanf(fp, "%lf", &Model->SphereOrigin[Y_]) != 1)
        cleanup("can't read sphere y-origin");
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'ORIGIN_Z' label");
    } while(strncmp(Label, "ORIGIN_Z:", 9) != 0);
    if(fscanf(fp, "%lf", &Model->SphereOrigin[Z_]) != 1)
        cleanup("can't read sphere z-origin");
    for(v=X_; v<=Z_; v++) {                 // convert to metres
        Model->SphereOrigin[v] *= 0.01;
        if(Model->SphereOrigin[v] > 0.1)    // origin must be < 10 cm
            cleanup("origin >= 10 cm!");
    }
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'RADIUS:' label");
    } while(strncmp(Label, "RADIUS:", 7) != 0);
    if(fscanf(fp, "%lf", &Model->SphereRadius) != 1)
        cleanup("can't read sphere radius");
    Model->SphereRadius *= 0.01;            // convert to metres
    if(Model->SphereRadius > 0.12)          // radius must be < 12 cm
        cleanup("radius > 12 cm!");

    // read voxel resolution
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'SAGITTAL:' label");
    } while(strncmp(Label, "SAGITTAL:", 9) != 0);
    if(fscanf(fp, "%lf", &Model->MriRes[S_]) != 1)      // mm resolution
        cleanup("can't read MRI sagittal resolution");
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'CORONAL:' label");
    } while(strncmp(Label, "CORONAL:", 8) != 0);
    if(fscanf(fp, "%lf", &Model->MriRes[C_]) != 1)      // mm resolution
        cleanup("can't read MRI coronal resolution");
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'AXIAL:' label");
    } while(strncmp(Label, "AXIAL:", 6) != 0);
    if(fscanf(fp, "%lf", &Model->MriRes[A_]) != 1)      // mm resolution
        cleanup("can't read MRI axial resolution");

    // read fiduciary coordinates into 3x3 matrix (in voxels)
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'NASION:' label");
    } while(strncmp(Label, "NASION:", 7) != 0);
    if(fscanf(fp, "%d%d%d", &Model->Fiducial[NA_][S_], &Model->Fiducial[NA_][C_], &Model->Fiducial[NA_][A_]) !=3)
        cleanup("can't read nasion voxel coordinates");
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'LEFT_EAR:' label");
    } while(strncmp(Label, "LEFT_EAR:", 9) != 0);
    if(fscanf(fp, "%d%d%d", &Model->Fiducial[LE_][S_], &Model->Fiducial[LE_][C_], &Model->Fiducial[LE_][A_]) !=3)
        cleanup("can't read left ear voxel coordinates");
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'RIGHT_EAR:' label");
    } while(strncmp(Label, "RIGHT_EAR:", 10) != 0);
    if(fscanf(fp, "%d%d%d", &Model->Fiducial[RE_][S_], &Model->Fiducial[RE_][C_], &Model->Fiducial[RE_][A_]) !=3)
        cleanup("can't read right ear voxel coordinates");

    // read multisphere fields
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'SEARCH_RADIUS:' label");
    } while(strncmp(Label, "SEARCH_RADIUS:", 14) != 0);
    if(fscanf(fp, "%lf", &Model->SearchRadius) !=1)
        cleanup("can't read search radius value");
    Model->SearchRadius *= 0.01;            // convert to metres
    do {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'HEADSHAPE_FILE:' label");
    } while(strncmp(Label, "HEADSHAPE_FILE:", 15) != 0);
    if(fscanf(fp, "%s", Name) != 1)
        cleanup("can't read head-shape name");
#if 0
    token = strtok(Name, "/");      // parse headshape file name from full path string
    for(;;) {
        token = strtok(NULL, "/");
        if(token != NULL)
            strcpy(Model->HeadShapeName, token);
        else
            break;
    }

    if(fscanf(fp, "%s", Label) != 1)        // read next line for optional indicator
        cleanup("can't read 'SURFACE_TYPE:' label");
    if(strncmp(Label, "SURFACE_TYPE:", 13) != 0) {
        Model->SurfaceType = SCALP;         // assume scalp if this field not present
    } else {
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read surface type label");
        if(!strncmp(Label, "SCALP", 5))
            Model->SurfaceType = SCALP;
        else if(!strncmp(Label, "CORTICAL", 8))
            Model->SurfaceType = CORTICAL;
        else if(!strncmp(Label, "INNERSKULL", 10))
            Model->SurfaceType = INNERSKULL;
        else
            cleanup("unknown surface type");
    }
#endif

    strcpy(Model->HeadShapeName, "hull.shape");
    Model->SurfaceType = INNERSKULL;

    // read values for each channel in multisphere list
    do {                                            // synch up to start of multisphere list
        if(fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'Radius' label");
    } while(strncmp(Label, "Radius", 6) != 0);      // this is the last word in header line
    M = 0;
    do {
        count = fscanf(fp, "%s%lf%lf%lf%lf", Name, &x, &y, &z, &r);
        if(count == 5) {
            strncpy(Model->MultiSphere[M].ChanName, Name, sizeof(Model->MultiSphere[M].ChanName));
            Model->MultiSphere[M].Origin[X_] = 0.01 * x;
            Model->MultiSphere[M].Origin[Y_] = 0.01 * y;
            Model->MultiSphere[M].Origin[Z_] = 0.01 * z;
            Model->MultiSphere[M].Radius = 0.01 * r;
            M++;
        }
    } while(count == 5 && M < 500);
    Model->NumSensors = M;

    // exit stage left...
    fclose(fp);
}
