// Some utility routines to make GIFTI file handling easier.

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "samutil.h"

#ifdef HAVE_GIFTI

#include "gifti_io.h"

gifti_image *new_gifti_image()
{
    gifti_image *gim;

    gim = new(gifti_image);
    gifti_clear_gifti_image(gim);
    gim->version = copy_string(GIFTI_XML_VERSION);

    return gim;
}

giiDataArray *new_gifti_DA(int intent, int datatype)
{
    giiDataArray *da;

    da = new(giiDataArray);
    gifti_clear_DataArray(da);
    da->intent = intent;
    da->datatype = datatype;
    da->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    da->encoding = GIFTI_ENCODING_B64GZ;
    da->endian = gifti_get_this_endian();
    return da;
}

#endif /* HAVE_GIFTI */
