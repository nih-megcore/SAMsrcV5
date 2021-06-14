/* Calculate some common stats from a file of numbers. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <math.h>

char *Progname;
char Usage[] = "usage: %s [-q] [-p <precision>] [file]\n";
#define USAGE() msg(Usage, Progname)

int Quiet = 0;
int Prec = 8;

void doit(FILE *);
void getstats(double *fp, int n, double *mp, double *vp, double *adp);

FILE *fileopen(char *name, char *mode);
char *fgetline(char *buf, int maxlen, FILE *infile);

void *new_block_(void *ptr, size_t s);
#define new_block(ptr, type, size) (type *)new_block_((void *)(ptr), \
	(size) * sizeof(type))

void fatalerr(char *, ...);
void msg(char *, ...);

int main(int argc, char **argv)
{
	int c;
	FILE *infile;

	Progname = *argv;

	/* Parse option arguments. */

	while ((c = getopt(argc, argv, "qp:")) != EOF) {
		switch (c) {

		case 'q':
			Quiet = 1;
			break;

		case 'p':
			Prec = atoi(optarg);
			if (Prec < 0) {
				Prec = 0;
			}
			if (Prec > 15) {
				Prec = 15;
			}
			break;

		default:
			USAGE();
			exit(1);
		}
	}
	argc -= optind;
	argv += optind;

	infile = stdin;
	if (argc == 1) {
		infile = fileopen(argv[0], "r");
	}

	doit(infile);

	return 0;
}

#define BLOCK 5000

void doit(FILE *infile)
{
	int n, remain;
	char buf[100];
	double *fp0, *fp, x, sum, min, max, mean, variance, adev;

	fp0 = NULL;
	n = remain = 0;
	min = 1.0e30;
	max = -1.0e30;
	sum = 0.;

	while (fgetline(buf, sizeof(buf), infile)) {
		x = atof(buf);
		if (remain == 0) {
			remain = BLOCK;
			fp0 = new_block(fp0, double, n + remain);
			fp = fp0 + n;
		}
		*fp++ = x;
		sum += x;
		n++;
		remain--;
		if (x > max) max = x;
		if (x < min) min = x;
	}

	/* Get the stats. */

	getstats(fp0, n, &mean, &variance, &adev);

	if (Quiet) {
		printf("%d %.*10$g %.*10$g %.*10$g %.*10$g %.*10$g %.*10$g %.*10$g %.*10$g\n",
			n, mean, variance, sum, min, max, adev,
			sqrt(variance), sqrt(variance) / sqrt((double)n), Prec);
	} else {
		printf("n = %d\nmean = %g\nvariance = %g\n", n, mean, variance);
		printf("sum is %g\n", sum);
		printf("range is [%g, %g]\n", min, max);
		printf("avg. dev. = %g\n", adev);
		printf("std. dev. = %g\nstd. error = %g\n", sqrt(variance),
			sqrt(variance) / sqrt((double)n));
	}

	free(fp0);
}

/* Return the mean, variance, and average deviation. */

void getstats(double *fp0, int n, double *mp, double *vp, double *adp)
{
	int i;
	double *fp, sum, sumsq, t;

	/* First, get the mean. */

	fp = fp0;
	sum = 0.;
	for (i = 0; i < n; i++) {
		sum += *fp++;
	}
	*mp = sum / n;

	/* Now the average deviation and the variance. */

	fp = fp0;
	sum = sumsq = 0.;
	for (i = 0; i < n; i++) {
		t = *fp++ - *mp;
		sum += fabs(t);
		sumsq += t * t;
	}
	*adp = sum / n;
	if (n == 1) n = 2;
	*vp = sumsq / (n - 1);
}

/* Open a file or exit on failure. */

FILE *fileopen(char *name, char *mode)
{
	FILE *f;

	if ((f = fopen(name, mode)) == NULL) {
		fatalerr("can't %s '%s'", *mode == 'r' ? "open" : "write", name);
	}
	return (f);
}

/* Like fgets() but remove the trailing newline.  Return a pointer to
the nul at the end of the string or NULL on error. */

char *fgetline(char *buf, int maxlen, FILE *infile)
{
	int i, c;
	char *s;

	s = buf;
	i = maxlen - 1;
	if (i < 0) {
		return (NULL);
	}
	if (i == 0) {
		*s = '\0';
		return (s);
	}
	c = getc(infile);
	if (c == EOF) {
		return (NULL);
	}
	while (c != '\n' && c != EOF) {
		*s++ = c;
		if (--i <= 0) {
			break;
		}
		c = getc(infile);
	}
	*s = '\0';
	return (s);
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

/* Print a message and exit. */

void fatalerr(char *str, ...)
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

void msg(char *str, ...)
{
	va_list ap;

	va_start(ap, str);
	vfprintf(stderr, str, ap);
	va_end(ap);
}
