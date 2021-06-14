// GetGIFTIWts() -- read beamformer weights from a GIFTI file.

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

gsl_matrix *GetGIFTIWts(
    char        *WgtFilePath
) {
    fatalerr("GIFTI format is unavailable! Install GIFTI support and recompile.");
}

#else

#include "gifti_io.h"

gsl_matrix *GetGIFTIWts(
    char        *WgtFilePath            // full path for GIFTI file
) {
    int                 V;              // number of voxels
    int                 M;              // number of channels
    int                 i, v, m;        // array indices
    gifti_image         *gim;           // gifti image pointer
    giiDataArray        *da;            // gifti data array pointer
    float               *fp;            // pointer into data array
    gsl_matrix          *Wgt;           // output weights

    // Fix the filename. .nii => .gii
    i = strlen(WgtFilePath);
    if (i > 4 && WgtFilePath[i-3] == 'n') {
        WgtFilePath[i-3] = 'g';
    }

    gifti_set_verb(0);                  // suppress coordsys warning
    gim = gifti_read_image(WgtFilePath, 1);
    if (gim == NULL) {
        fatalerr("Failed to read gifti file '%s'.", WgtFilePath);
    }
    assert(gim->numDA == 1);
    da = gim->darray[0];
    assert(da->intent == NIFTI_INTENT_VECTOR);
    assert(da->datatype == NIFTI_TYPE_FLOAT32);
    assert(da->num_dim == 2);

    V = da->dims[0];
    M = da->dims[1];
    Wgt = gsl_matrix_alloc(V, M);

    fp = (float *)da->data;
    for (v = 0; v < V; v++) {
        for (m = 0; m < M; m++) {
            gsl_matrix_set(Wgt, v, m, *fp++);
        }
    }

    gifti_free_image(gim);
    return Wgt;
}

#endif /* HAVE_GIFTI */
