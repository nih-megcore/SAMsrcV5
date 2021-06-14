// Get the _Noise file value.

#include <stdio.h>
#include "samutil.h"

void GetNoise(char *dir, char *name, double *r)
{
    int i;
    FILE *np;
    char *fpath, *s;

    // space to work

    i = strlen(dir) + strlen(name) + 9;     /* %s/%s_MuNoise */
    fpath = new_string(i);

    // Look for _MuNoise first. If that fails, try _Noise

    s = fpath;
    s = strecpy(s, dir);
    *s++ = '/';
    s = strecpy(s, name);
    s = strecpy(s, "_MuNoise");
    np = fopen(fpath, "r");
    if (np == NULL) {
        s = fpath;
        s = strecpy(s, dir);
        *s++ = '/';
        s = strecpy(s, name);
        s = strecpy(s, "_Noise");
        np = fopen(fpath, "r");
    }

    if (np == NULL)
        Cleanup("can't open noise file for %s (did you run sam_cov?)", name);
    if (fscanf(np, "%le", r) != 1)
        Cleanup("can't read noise file for %s", name);

    fclose(np);
    free(fpath);
}
