/* Memory utilities. */

#include "samutil.h"

/* Allocate some space and return a pointer to it. See new() in samutil.h. */

void *new_(size_t s, char *err)
{
    void *x;

    if ((x = (void *)malloc(s)) == NULL) {
        fatalerr("[%s] can't malloc(%d)", err, s);
    }

    return x;
}

/* Like strdup() but exits on error. */

char *copy_string(char *s)
{
    char *t;

    t = new_string(strlen(s));
    strcpy(t, s);

    return t;
}

/* Like strcpy(), but return a pointer to the new nul at the end of t. */

char *strecpy(char *t, char *s)
{
    while (*s) {
        *t++ = *s++;
    }
    *t = '\0';
    return t;
}

/* Allocate an n x m data matrix. */

double **new_matrixE(int n, int m, char *err)
{
    int i;
    double *d, **a;

    /* Allocate the array in contiguous storage. This allows some
    routines to work more efficiently by avoiding a pointer dereference
    on each access. */

    d = new_arrayE(double, n * m, err);

    /* Also return an array of pointers to the rows, for convenience
    where speed isn't an issue. The pointer to the whole array is
    stored in a[0] (or *a). Note that on some architectures, the
    pointer dereference may be cached and take very little time. */

    a = new_arrayE(double *, n, err);
    for (i = 0; i < n; i++) {
        a[i] = d;
        d += m;
    }

    return a;
}

/* Free a matrix. */

void free_matrix(double **a)
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
