#ifndef SAM_PLATFORM_H
#define SAM_PLATFORM_H

/* Compatibility shims for the native Windows (MinGW-w64) wheel build. */
#ifdef _WIN32
#include <direct.h>
#include <io.h>
#include <math.h>
#include <process.h>
#include <stdlib.h>
#include <string.h>

#ifndef F_OK
#define F_OK 0
#endif
#ifndef R_OK
#define R_OK 4
#endif
#ifndef W_OK
#define W_OK 2
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#define access _access
#define chdir _chdir
#define getcwd _getcwd
#define getpid _getpid
#define mkdir(path, mode) _mkdir(path)
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#define strdup _strdup
#define unlink _unlink
#define index strchr
#define rindex strrchr

long sam_random(void);
void sam_srandom(unsigned int seed);
double sam_drand48(void);
void sam_srand48(long seed);

#define random sam_random
#define srandom sam_srandom
#define drand48 sam_drand48
#define srand48 sam_srand48
#endif

#endif
