/* Error and warning messages. */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include "defs.h"

char *Progname = NULL;

/* Print a message and exit. */

void fatalerr(const char *str, ...)
{
	va_list ap;

	if (Progname) {
		fprintf(stderr, "%s: ", Progname);
	}
	va_start(ap, str);
	vfprintf(stderr, str, ap);
	va_end(ap);
	fputc('\n', stderr);

	exit(1);
}

/* Just print a message. */

void msg(const char *str, ...)
{
	va_list ap;

	va_start(ap, str);
	vfprintf(stderr, str, ap);
	va_end(ap);
}
