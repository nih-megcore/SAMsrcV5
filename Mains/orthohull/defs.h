/* Standard definitions. */

#ifndef _DEFS_H
#define _DEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef EPS
#define EPS 1e-7
#endif

#ifndef __COMPAR_FN_T
#define __COMPAR_FN_T
typedef int (*__compar_fn_t)(const void *, const void *);
#endif

#ifdef  __cplusplus
extern "C" {
#endif

/* A function and some macros for allocating memory. [mem.c] */

extern void *new_(size_t s);
extern void *new_block_(void *ptr, size_t s);

#define new(type) (type *)new_(sizeof(type))
#define new_string(size) (char *)new_((size) + 1)
#define new_array(type, size) (type *)new_((size) * sizeof(type))
#define new_block(ptr, type, size) (type *)new_block_((void *)(ptr), \
	(size) * sizeof(type))

/* General memory allocation. [mem.c] */

extern char *copy_string(const char *s);
extern float **new_matrix(int n, int m);
extern void free_matrix(float **);

/* Error messages. [err.c] */

extern char *Progname;
extern void fatalerr(const char *, ...);
extern void msg(const char *, ...);

#ifdef  __cplusplus
}
#endif

/* Handy typedefs for declaring storage of 2- & 3-vectors. */

typedef float XY[2];
typedef float XYZ[3];

#endif  /* _DEFS_H */
