#ifndef _SAMUTIL_H
#define _SAMUTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef HUGE
#define HUGE HUGE_VAL
#endif

// cleanup.c -- error and warning messages

extern char *Progname;
extern void Cleanup(char *, ...) __attribute__ ((noreturn));
extern void allocfailed(char *, ...) __attribute__ ((noreturn));
extern void msg(char *, ...);

#define cleanup Cleanup // backwards compat.
#define fatalerr Cleanup

// file.c -- file utilities

extern FILE *fileopen(char *name, char *mode);
extern int fileexists(char *name);
extern int direxists(char *name);

extern char *fgetline(char *buf, int maxlen, FILE *infile);

// mem.c -- functions and macros for allocating memory

extern void *new_(size_t s, char *err);
extern void *new_block_(void *ptr, size_t s);

#define new(type) (type *)new_(sizeof(type), "new")
#define newE(type, err) (type *)new_(sizeof(type), err)
#define new_string(size) (char *)new_((size) + 1, "string")
#define new_array(type, size) (type *)new_((size) * sizeof(type), "array")
#define new_arrayE(type, size, err) (type *)new_((size) * sizeof(type), err)
#define new_block(ptr, type, size) \
    (type *)new_block_((void *)(ptr), (size) * sizeof(type))

extern char *copy_string(char *s);
extern char *strecpy(char *t, char *s);
extern double **new_matrixE(int n, int m, char *err);
extern void free_matrix(double **);

// misc

void shuffle(int *array, int n);

// Environment variable access.

extern char **environ;

#endif
