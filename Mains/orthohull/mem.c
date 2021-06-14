/* Memory utilities. */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include "defs.h"

/* Allocate some space and return a pointer to it. See new() in defs.h. */

void *new_(size_t s)
{
	void *x;

	if ((x = (void *)malloc(s)) == NULL) {
		fatalerr("can't malloc(%d)", s);
	}

	return x;
}

/* Like strdup() but exits on error. */

char *copy_string(const char *s)
{
	char *t;

	t = new_string(strlen(s));
	strcpy(t, s);

	return t;
}

/* Allocate an n x m data matrix. */

float **new_matrix(int n, int m)
{
	int i;
	float *d, **a;

	/* Allocate the array in contiguous storage. This allows some
	routines to work more efficiently by avoiding a pointer dereference
	on each access. */

	d = new_array(float, n * m);

	/* Also return an array of pointers to the rows, for convenience
	where speed isn't an issue. The pointer to the whole array is
	stored in a[0] (or *a). Note that on some architectures, the
	pointer dereference may be cached and take very little time. */

	a = new_array(float *, n);
	for (i = 0; i < n; i++) {
		a[i] = d;
		d += m;
	}

	return a;
}

/* Free a matrix. */

void free_matrix(float **a)
{
	free(*a);
	free(a);
}

/* Allocate or reallocate an array using realloc(). */

void *new_block_(void *ptr, size_t s)
{
	void *x;

	if ((x = (void *)realloc(ptr, s)) == NULL) {
		fatalerr("can't realloc(%d)", s);
	}

	return x;
}
