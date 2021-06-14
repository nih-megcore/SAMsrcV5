// GetVoxelsGII() -- populate an array of voxels using a
//  FreeSurfer surface in GIFTI format, exported with mkGiiAtlas.py

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

#ifndef HAVE_GIFTI

GIIATLAS *GetGiiAtlas(
    char        *AtlasName
) {
    fatalerr("GIFTI format is unavailable! Install GIFTI support and recompile.");
}

#else

#include "gifti_io.h"

static void get_centroid(giiMetaData *meta, double *cent);
static int get_hemi(giiMetaData *meta);
#define RH 0
#define LH 1

GIIATLAS *GetGiiAtlas(
    char        *AtlasName
) {
    GIIATLAS            *atlas;         // returned atlas info
    VOXELINFO           *Voxel;         // array of voxels
    int                 *ridx;          // right hemisphere index
    int                 *lidx;          // left hemisphere index
    int                 nrh, nlh;       // number of right and left hemisphere voxels
    int                 V;              // total number of voxels
    int                 i, j, k;        // array indices
    int                 n;              // voxel index
    int                 v;              // vector index
    gifti_image         *gim;           // gifti image pointer
    giiDataArray        *da;            // gifti data array pointer
    int                 numDA;          // number of data arrays
    int                 hemi;           // hemisphere, RH or LH
    int                 *ip, *idx;      // pointer into node index arrays
    float               *fp;            // pointer into data arrays

    // Read the Atlas file.

    gifti_set_verb(0);                  // suppress coordsys warning
    gim = gifti_read_image(AtlasName, 1);
    if (gim == NULL) {
        fatalerr("Failed to read gifti file '%s'.", AtlasName);
    }
    numDA = gim->numDA;

    // Figure out how many voxels there are, sanity check the file,
    // and read the R and L indices.

    nrh = nlh = 0;
    for (i = 0; i < numDA; i++) {
        da = gim->darray[i];
        if (da->ind_ord != GIFTI_IND_ORD_ROW_MAJOR) {
            fatalerr("unsupported Atlas data array format");
        }
        if (da->intent == NIFTI_INTENT_NODE_INDEX) {
            hemi = get_hemi(&da->meta);
            n = (int)da->dims[0];
            idx = new_array(int, n);
            if (hemi == RH) {
                nrh = n;
                ridx = idx;
            } else {
                nlh = n;
                lidx = idx;
            }
            ip = (int *)da->data;
            for (j = 0; j < n; j++) {
                *idx++ = *ip++;
            }
        }
    }
    V = nrh + nlh;
    if (V == 0) {
        fatalerr("No data in Atlas file!");
    }

    // Allocate the voxel array and initialize the flags.

    Voxel = new_array(VOXELINFO, V);
    for (i = 0; i < V; i++) {
        Voxel[i].Solve = FALSE;
        Voxel[i].ROI = TRUE;
    }

    // Fill in the positions, and normals.
    // We create a single Voxel array from the atlas,
    // by convention mkGiiAtlas puts the RH first.
    // Do the POINTSET and VECTOR arrays separately
    // so the voxels are in the right order.

    k = 0;
    for (i = 0; i < numDA; i++) {
        da = gim->darray[i];
        if (da->intent == NIFTI_INTENT_POINTSET) {
            n = (int)da->dims[0];
            assert(da->dims[1] == 3);
            fp = (float *)da->data;
            for (j = 0; j < n; j++) {
                for (v = 0; v < 3; v++) {
                    Voxel[k].p[v] = *fp++;
                }
                k++;
            }
        }
    }

    k = 0;
    for (i = 0; i < numDA; i++) {
        da = gim->darray[i];
        if (da->intent == NIFTI_INTENT_VECTOR) {
            n = (int)da->dims[0];
            assert(da->dims[1] == 3);
            fp = (float *)da->data;
            for (j = 0; j < n; j++) {
                for (v = 0; v < 3; v++) {
                    Voxel[k].v[v] = *fp++;
                }
                k++;
            }
        }
    }

#if DEBUG
    printf("\nread %d voxels from '%s'\n", V, AtlasName);
#endif

    // create gifti atlas structure

    atlas = new(GIIATLAS);
    atlas->Voxel = Voxel;
    atlas->ridx = ridx;
    atlas->lidx = lidx;
    atlas->nrh = nrh;
    atlas->nlh = nlh;
    atlas->meta = new(giiMetaData);
    gifti_clear_nvpairs(atlas->meta);
    gifti_copy_nvpairs(atlas->meta, &gim->meta);
    get_centroid(atlas->meta, atlas->cent);

    gifti_free_image(gim);

    return atlas;
}

/* Scan a DA's metadata, get the hemisphere code. */

static int get_hemi(giiMetaData *meta)
{
    int i;

    for (i = 0; i < meta->length; i++) {
        if (strcmp(meta->name[i], "AnatomicalStructurePrimary") == 0) {
            if (strcmp(meta->value[i], "CortexRight") == 0) {
                return RH;
            } else if (strcmp(meta->value[i], "CortexLeft") == 0) {
                return LH;
            }
        }
    }
    fatalerr("Improper Atlas file, no metadata.");
}

/* Copy the centroid metadata. */

static void get_centroid(giiMetaData *meta, double *cent)
{
    int i;

    for (i = 0; i < meta->length; i++) {
        if (strcmp(meta->name[i], "centroid") == 0) {
            if (sscanf(meta->value[i], "%lf%lf%lf", cent, cent + 1, cent + 2) != 3) {
                fatalerr("improper centroid in atlas metadata");
            }
            return;
        }
    }
    fatalerr("Improper Atlas file, no centroid.");
}

#endif /* HAVE_GIFTI */
