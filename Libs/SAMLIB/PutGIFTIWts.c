// PutGIFTIWts() -- write beamformer weights to a GIFTI file.

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

void PutGIFTIWts(
    char        *WgtFilePath,
    gsl_matrix  *Wgt,
    float       Noise
) {
    fatalerr("GIFTI format is unavailable! Install GIFTI support and recompile.");
}

#else

#include "gifti_io.h"

extern gifti_image *new_gifti_image();
extern giiDataArray *new_gifti_DA(int intent, int datatype);
static void set_meta(giiMetaData *meta, char *name, char *value);
static void make_vector_da(giiDataArray *da, gsl_matrix *w);

void PutGIFTIWts(
    char        *WgtFilePath,           // full path for writing GIFTI file
    gsl_matrix  *Wgt,                   // Wgt - VxM SAM coefficients
    float       Noise                   // sensor mean noise floor (unused?)
) {
    int                 i;              // array index
    gifti_image         *gim;           // gifti image pointer
    char                buf[100];

    //gifti_set_verb(10);

    // Fix the filename. .nii => .gii
    i = strlen(WgtFilePath);
    if (i > 4 && WgtFilePath[i-3] == 'n') {
        WgtFilePath[i-3] = 'g';
    }

    // We'll write everything out as one data array.
    // Later the node indices from the atlas will be
    // used to split the results into two files.

    gim = new_gifti_image();
    gim->numDA = 1;
    gim->darray = new_array(giiDataArray *, gim->numDA);
    gim->darray[0] = new_gifti_DA(NIFTI_INTENT_VECTOR, NIFTI_TYPE_FLOAT32);

    make_vector_da(gim->darray[0], Wgt);

    // Store the noise in the metadata for the image.
    // @@@ It would be better to store this as a 1-element data array.

    snprintf(buf, sizeof(buf), "%g", Noise);
    set_meta(&gim->meta, "noise", buf);

    // Write it.

    if (gifti_write_image(gim, WgtFilePath, 1) == 1) {
        fatalerr("there was a problem writing '%s'", WgtFilePath);
    }
    gifti_free_image(gim);
}

/* Set a meta field. Only 1 name/value pair here! */

static void set_meta(giiMetaData *meta, char *name, char *value)
{
    meta->length = 1;
    meta->name = new(char *);
    meta->value = new(char *);
    meta->name[0] = copy_string(name);
    meta->value[0] = copy_string(value);
}

/* Fill in a vector array. */

static void make_vector_da(giiDataArray *da, gsl_matrix *w)
{
    int v, m, V, M;
    float *fp;

    da->num_dim = 2;
    da->dims[0] = V = w->size1;
    da->dims[1] = M = w->size2;
    da->nvals = V * M;
    da->nbyper = 4;
    fp = new_array(float, V * M);
    da->data = fp;
    for (v = 0; v < V; v++) {
        for (m = 0; m < M; m++) {
            *fp++ = (float)gsl_matrix_get(w, v, m);
        }
    }
}

#endif /* HAVE_GIFTI */
