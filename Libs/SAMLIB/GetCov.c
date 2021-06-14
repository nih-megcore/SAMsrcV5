// GetCov()
//
//      Author: Stephen E. Robinson
//              MEG Core Facility
//              NIMH
//

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <fcntl.h>
#include <gsl/gsl_matrix.h>
#include <SAMfiles.h>
#include <samlib.h>
#include "samutil.h"
#include <version.h>

void GetCov(
    char        *CovName,           // fully-qualified path name
    COV_HDR     *CovHeader,         // the covariance header
    int         **ChanOut,          // ChanIndex[M] -- out
    gsl_matrix  *CovOut             // Covariance - MxM
) {
    double          Ctmp;           // covariance element
    int             *ChanIndex;     // ChanIndex[M]
    int             i;              // channel index
    int             j;              // a useful, but entirely different index
    int             M;              // number of channels
    char            Identity[16];   // file identification string
    char            CovID[16] =     // covariance identity string
        "SAMCOVAR";
    off_t           offset;         // bytes offset from start of file
    struct stat     sbuf;           // stat buffer
    FILE            *fin;           // input data pointer

    // open file for read
    fin = fileopen(CovName, "r");
    stat(CovName, &sbuf);
    if (fread((void *)Identity, 8, 1, fin) != 1)
        cleanup("can't read covariance ID");
    if (strncmp(Identity, CovID, 8) != 0)
        cleanup("incorrect covariance ID string");
    if (fread((void *)CovHeader, sizeof(COV_HDR), 1, fin) != 1)
        cleanup("can't read covariance header");

    M = CovHeader->NumChans;
    ChanIndex = new_array(int, M);

    if (CovHeader->Version != SAM_REV)
        cleanup("covariance file has incorrect version number");
    offset = sbuf.st_size - (M * 4 + M * M * 8);
    if (fseek(fin, offset, SEEK_SET) == -1)
        cleanup("can't seek to channel index");
    if (fread((void *)ChanIndex, sizeof(int), M, fin) != M)
        cleanup("can't read channel index");
    offset = sbuf.st_size - M * M * 8;
    if (fseek(fin, offset, SEEK_SET) == -1)
        cleanup("can't seek to covariance matrix");
    for (i=0; i<M; i++)
        for (j=0; j<M; j++) {
            if (fread((void *)&Ctmp, sizeof(double), 1, fin) != 1)
                cleanup("can't read covariance matrix");
            gsl_matrix_set(CovOut, i, j, Ctmp);
        }
    *ChanOut = ChanIndex;
    fclose(fin);
}
