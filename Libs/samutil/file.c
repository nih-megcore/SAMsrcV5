/* File utilities. */

#include "samutil.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* Open a file or exit on failure. */

FILE *fileopen(char *name, char *mode)
{
    FILE *f;

    if ((f = fopen(name, mode)) == NULL) {
        fatalerr("can't %s '%s'", *mode == 'r' ? "open" : "write", name);
    }
    return f;
}

/* Check to see if a file exists before we open it. */

int fileexists(char *name)
{
    struct stat buf;

    if (stat(name, &buf) < 0 || (buf.st_mode & S_IFMT) != S_IFREG) {
        return FALSE;
    }
    return TRUE;
}

/* Check to see if a directory exists. */

int direxists(char *name)
{
    struct stat buf;

    if (stat(name, &buf) < 0 || (buf.st_mode & S_IFMT) != S_IFDIR) {
        return FALSE;
    }
    return TRUE;
}

/* Like fgets() but remove the trailing newline. Return a pointer to
the nul at the end of the string or NULL on error. */

char *fgetline(char *buf, int maxlen, FILE *infile)
{
    int i, c;
    char *s;

    s = buf;
    i = maxlen - 1;
    if (i < 0) {
        return NULL;
    }
    if (i == 0) {
        *s = '\0';
        return s;
    }
    c = getc(infile);
    if (c == EOF) {
        return NULL;
    }
    while (c != '\n' && c != EOF) {
        *s++ = c;
        if (--i <= 0) {
            break;
        }
        c = getc(infile);
    }
    *s = '\0';
    return s;
}
