// Construct various pathnames using patterns.

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include "samlib.h"
#include "samutil.h"
#include "sam_param.h"

// Copy pattern, with replacements, into path.

// For example, a pattern of "%M/%P/%s" is the default for files
// in the MRI directory.

// ~ is the value of $HOME
// %M is p->MRIDirectory
// %P is setname[:p->NumPrefix]
// %s is the passed name
// %d is the dataset name
// %H is the hash code (first _ delimited field of setname)
// %S is the study (second _ field)
// %D is the date (third)
// %R is the run (fourth)

// path is output, using the given pattern. name may override the pattern.
// If exist is set, the resulting path must exist.
// len is the size of the path buffer, and must be large enough.

void GetMRIPath(char *path, int len, PARMINFO *p, char *name, int exist)
{
    GetFilePath(p->MRIPattern, path, len, p, name, exist);
}

void GetFilePath(char *pattern, char *path, int len, PARMINFO *p, char *name, int exist)
{
    int i;
    char *s, *t, *setname, *home;
    char *f[4];     // H, S, D, and R pointers

    // If name is absolute or contains metacharacters, it overrides the pattern.

    i = 0;
    for (s = name; *s; s++) {
        if (*s == '%') {
            i = 1;
        }
    }
    if (i || name[0] == '/' || name[0] == '~') {
        pattern = name;
        name = "";          // so %s does nothing
    }

    // Parse the setname, for %H, etc. Up to 4 underscore delimited fields.

    s = p->DataSetName;
    t = rindex(s, '/');     // just the setname
    if (t != NULL) {
        s = t + 1;
    }
    setname = s;

    t = new_string(strlen(s) + 1);
    f[0] = t;
    i = 1;
    while (*s) {
        if (*s == '_' && i < 4) {
            *t++ = '\0';
            f[i] = t;
            i++;
        } else {
            *t++ = *s;
        }
        s++;
    }
    *t = '\0';

    // Parse the pattern, output what it says.

    s = pattern;
    t = path;

    if (*s == '~') {
        home = getenv("HOME");
        if (home) {
            t = strecpy(t, home);
        }
        s++;
    }

    while (*s) {
        if (*s != '%') {
            *t++ = *s++;
        } else {
            s++;
            if (*s == 'M') {
                t = strecpy(t, p->MRIDirectory);
            } else if (*s == 'P') {
                strncpy(t, setname, p->NumPrefix);
                t += p->NumPrefix;
                *t = '\0';
            } else if (*s == 's') {
                t = strecpy(t, name);
            } else if (*s == 'd') {
                t = strecpy(t, p->DataSetName);
            } else if (*s == 'H') {
                t = strecpy(t, f[0]);
            } else if (*s == 'S') {
                t = strecpy(t, f[1]);
            } else if (*s == 'D') {
                t = strecpy(t, f[2]);
            } else if (*s == 'R') {
                t = strecpy(t, f[3]);
            }
            s++;
        }
    }
    *t = '\0';
    if (t >= path + len) {  // we already borked it
        fatalerr("filename too long: %s", path);
    }
    free(f[0]);

    if (exist) {
        if (access(path, F_OK) != 0) {
            fatalerr("can't access %s", path);
        }
    }
}
