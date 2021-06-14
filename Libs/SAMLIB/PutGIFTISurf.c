// PutGIFTISurf() -- Write GIFTI format surface files that SUMA can read.

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

void PutGIFTISurf(char *surfname, float **d, int nd, int start, int *idx, int n, void *meta, float TimeStep)
{
    fatalerr("GIFTI format is unavailable! Install GIFTI support and recompile.");
}

#else

#include "gifti_io.h"

extern gifti_image *new_gifti_image();
extern giiDataArray *new_gifti_DA(int intent, int datatype);

/* This routine writes GIFTI surface files. It can write timeseries data,
which is stored as (index, da1, da2, ..., daN), or a regular image with
just one data array. TimeStep is ignored in that case, otherwise it is
added to the file's metadata. The meta data is copied from the original
Atlas file to preserve the original surface information, which is
useful for converting back to a non-indexed version (for smoothing). */

void PutGIFTISurf(
    char        *surfname,              // pathname for writing GIFTI file
    float       **d,                    // array(s) of floats to write from
    int         nd,                     // number of arrays in d
    int         start,                  // index of the first voxel to write
    int         *idx,                   // node index
    int         n,                      // number of voxels to write
    void        *meta,                  // metadata from original atlas
    float       TimeStep                // time step if nd > 1
) {
    int                 i, j;           // array indices
    int                 *ip;            // node index pointer
    float               *fp;            // data pointer
    char                buf[100];       // temp string buffer
    gifti_image         *gim;           // gifti image pointer
    giiDataArray        *da;            // data array

    // Output surface files have a node index and at least one float array.

    gim = new_gifti_image();
    gim->numDA = nd + 1;
    gim->darray = new_array(giiDataArray *, gim->numDA);
    gim->darray[0] = new_gifti_DA(NIFTI_INTENT_NODE_INDEX, NIFTI_TYPE_INT32);
    i = NIFTI_INTENT_NONE;              // intent for a single image
    if (nd > 1) {
        i = NIFTI_INTENT_TIME_SERIES;   // ... and for more than one
    }
    for (j = 0; j < nd; j++) {
        gim->darray[j + 1] = new_gifti_DA(i, NIFTI_TYPE_FLOAT32);
    }

    // node index

    da = gim->darray[0];
    da->num_dim = 1;
    da->dims[0] = n;
    da->nvals = n;
    da->nbyper = 4;
    ip = new_array(int, n);
    da->data = ip;
    for (i = 0; i < n; i++) {
        *ip++ = idx[i];
    }

    // float array(s)

    for (j = 0; j < nd; j++) {
        da = gim->darray[j + 1];
        da->num_dim = 1;
        da->dims[0] = n;
        da->nvals = n;
        da->nbyper = 4;
        fp = new_array(float, n);
        da->data = fp;
        for (i = 0; i < n; i++) {
            *fp++ = d[j][start + i];
        }
    }

    // Copy the passed image metadata. Update the metadata to
    // include TimeStep, if needed.

    gifti_copy_nvpairs(&gim->meta, (giiMetaData *)meta);
    if (nd > 1) {
        snprintf(buf, sizeof(buf), "%f", TimeStep);
        gifti_add_to_nvpairs(&gim->meta, "TimeStep", buf);
    }

    // Write it.

    if (gifti_write_image(gim, surfname, 1) == 1) {
        fatalerr("there was a problem writing '%s'", surfname);
    }

    gifti_free_image(gim);
}

#endif /* HAVE_GIFTI */
