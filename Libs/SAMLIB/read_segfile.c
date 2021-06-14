/* Read a "segment" file, which contains start and end times for segments to
use in the covariance. */

#include "samutil.h"
#include <DataFiles.h>
#include <TimeSegs.h>

/* Read a file of segments, trial t0 t1, and return a pointer to a COV_SEG
array with nsegs elements. nsegs is set to the number of entries. */

COV_SEG *read_segfile(char *name, int *nsegs, HeaderInfo *header)
{
    int n, tr;
    char *s, *e;
    double t0, t1;
    FILE *f;
    COV_SEG *segs, *seg;
    char buf[256];

    /* Count the lines. */

    n = 0;
    f = fileopen(name, "r");
    while (fgetline(buf, sizeof(buf), f)) {
        n++;
    }
    rewind(f);

    /* Read them into a COV_SEG array. */

    segs = new_array(COV_SEG, n);
    *nsegs = n;
    seg = segs;
    while (fgetline(buf, sizeof(buf), f)) {
        s = buf;

        /* trial t0 t1 */

        tr = (int)strtol(s, &e, 10);
        if (e == s) {
            fatalerr("badly formed trial # in segment file %s: '%s'", name, buf);
        }
        s = e;

        t0 = strtod(s, &e);
        if (e == s) {
            fatalerr("badly formed t0 in file %s: '%s'", name, buf);
        }
        s = e;

        t1 = strtod(s, &e);
        if (e == s) {
            fatalerr("badly formed t1 in file %s: '%s'", name, buf);
        }

        seg->Epoch = tr;
        seg->TS = (int)rint(header->SampleRate * t0) + header->PreTrigSamples;
        seg->TE = (int)rint(header->SampleRate * t1) + header->PreTrigSamples;
        seg->TT = seg->TE - seg->TS;
        seg++;
    }

    fclose(f);
    return segs;
}
