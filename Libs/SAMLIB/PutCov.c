// PutCov() -- write covariance file (with header) to disk
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <SAMfiles.h>
#include <samlib.h>
#include <version.h>

void PutCov(
    char    *CovName,       // fully-qualified path name
    COV_HDR *CovHeader,     // the filled-in header
    int     *ChanIndex,     // channel index (base-0)
    double  **Covariance    // Covariance[M][M]
) {
    int     i, j;           // channel indices
    int     M;              // number of primary sensors used
    int32_t *WriteIndex;    // copy of channel index array to write
    char    Identity[16] =  // covariance identity string
                "SAMCOVAR";
    FILE    *fout;          // output file pointer

    // allocate array for indexing channels
    M = CovHeader->NumChans;
    WriteIndex = new_arrayE(int32_t, M, "WriteIndex[]");    // make sure it's int32_t

    // open file for write
    fout = fileopen(CovName, "w");
    if (fwrite((void *)Identity, 8, 1, fout) != 1)
        cleanup("covariance ID string write failed");
    if (fwrite((void *)CovHeader, sizeof(COV_HDR), 1, fout) != 1)
        cleanup("covariance header write failed");
    for (i=0; i<M; i++)
        WriteIndex[i] = ChanIndex[i];
    if (fwrite((void *)WriteIndex, sizeof(int32_t), M, fout) != M)
        cleanup("covariance index write failed");
    for (i=0; i<M; i++)
        for (j=0; j<M; j++)
            if (fwrite((void *)&Covariance[i][j], sizeof(double), 1, fout) != 1)
                cleanup("covariance matrix write failed");
    fclose(fout);

    free(WriteIndex);
}
