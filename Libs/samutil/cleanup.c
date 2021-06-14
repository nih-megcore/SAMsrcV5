// Cleanup routine

#include <lfulib.h>
#include <stdarg.h>

char *Progname = "library";     // default for backwards compat.

/* Fatal error, print a message and exit. */

void Cleanup(char *str, ...)
{
    va_list ap;

    fprintf(stderr, "\n\07%s: ", Progname);
    va_start(ap, str);
    vfprintf(stderr, str, ap);
    va_end(ap);
    fputc('\n', stderr);

    exit(-1);
}

/* Fatal allocation error, print a message and exit. */

void allocfailed(char *str, ...)
{
    va_list ap;

    fprintf(stderr, "\n\07%s: allocation failure -- ", Progname);
    va_start(ap, str);
    vfprintf(stderr, str, ap);
    va_end(ap);
    fputc('\n', stderr);

    exit(-1);
}

/* Just print an information or warning message to stderr. */

void msg(char *str, ...)
{
    va_list ap;

    va_start(ap, str);
    vfprintf(stderr, str, ap);
    va_end(ap);
    fflush(stderr);
}
