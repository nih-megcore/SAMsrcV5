/* Read a triangle mesh and make normals for each vertex. */

#include "defs.h"
#include <unistd.h>

typedef struct {
	XYZ v;                          /* vertex */
	XYZ n;                          /* normal */
	struct vlist_struct *vlist;     /* and a list of adjacent vertices */
} VERTEX;

typedef struct vlist_struct {
	VERTEX *vp;
	struct vlist_struct *next;
} VLIST;

VERTEX *Verts;
int Nvert;

typedef struct {                        /* triangle */
	VERTEX *vp[3];                  /* 3 vertices */
} TRIANGLE;

TRIANGLE *Triangles;
int Ntri;

void doit(FILE *infile);
void read_mesh(FILE *infile);
void inflate_mesh(float frac);
void build_normals(void);
void make_edge(VERTEX *a, VERTEX *b);
void write_mesh(void);
static inline float triareaH(float *v1, float *v2, float *v3);

int iflag = FALSE;
float frac;

char Usage[] = "usage: %s [-i<inflate_frac>] [file.off]\n";
#define USAGE() msg(Usage, Progname)

int main(int argc, char **argv)
{
	int c;
	FILE *infile;

	Progname = *argv;

	/* Parse option arguments. */

	while ((c = getopt(argc, argv, "i:")) != EOF) {
		switch (c) {

		case 'i':
			frac = atof(optarg);
			iflag = TRUE;
			break;

		default:
			USAGE();
			exit(1);
		}
	}
	argc -= optind;
	argv += optind;

	/* Process the named file, or stdin if no file given.
	The name '-' also specifies stdin. */

	infile = stdin;
	if (argc && **argv != '-') {
		infile = fopen(*argv, "r");
		if (infile == NULL) {
			fatalerr("can't open '%s'\n", *argv);
		}
	}

	doit(infile);

	return 0;
}

void doit(FILE *infile)
{
	read_mesh(infile);

	build_normals();

	if (iflag) {
		inflate_mesh(frac);
	}

	write_mesh();
}

/* Read a triangle mesh from a file. */

void read_mesh(FILE *infile)
{
	int i, a, b, c, n;
	float *v, x;
	char buf[100];
	VERTEX *vp;
	TRIANGLE *tp;
	XYZ norm;

	/* This reads .off files. */

	fgets(buf, 100, infile);
	if (strcmp(buf, "OFF\n") != 0) {
		fatalerr("not an OFF file\n");
	}

	fgets(buf, 100, infile);
	sscanf(buf, "%d %d %*d", &Nvert, &Ntri);

	/* Read the verts. */

	Verts = new_array(VERTEX, Nvert);
	for (i = 0, vp = Verts; i < Nvert; i++, vp++) {
		v = vp->v;
		vp->vlist = NULL;
		vzero(vp->n);                   /* init the normal */
		fgets(buf, 100, infile);
		sscanf(buf, "%f %f %f", v, v + 1, v + 2);

#if 0
		/* Convert from RAI mm to PRI m. */

		x = v[0];
		v[0] = -v[1];
		v[1] = x;
		vscale(v, .001);
#endif
	}

	/* Next is the triangles. Keep track of the edges. */

	Triangles = new_array(TRIANGLE, Ntri);
	for (i = 0, tp = Triangles; i < Ntri; i++, tp++) {
		fgets(buf, 100, infile);
		sscanf(buf, "%*d %d %d %d", &c, &b, &a);
		tp->vp[0] = &Verts[a];
		tp->vp[1] = &Verts[b];
		tp->vp[2] = &Verts[c];
		make_edge(&Verts[a], &Verts[b]);
		make_edge(&Verts[b], &Verts[c]);
		make_edge(&Verts[c], &Verts[a]);
	}
}

void make_edge(VERTEX *a, VERTEX *b)
{
	VLIST *vlp;

	/* See if a already points to b. */

	for (vlp = a->vlist; vlp; vlp = vlp->next) {
		if (vlp->vp == b) {
			return;
		}
	}

	/* Make two new VLIST structures, and insert each at the head of
	the other's list.  */

	vlp = new(VLIST);
	vlp->vp = b;
	vlp->next = a->vlist;
	a->vlist = vlp;

	vlp = new(VLIST);
	vlp->vp = a;
	vlp->next = b->vlist;
	b->vlist = vlp;
}

void build_normals(void)
{
	int i, j;
	float x;
	VERTEX *a, *b, *c;
	TRIANGLE *tp;
	XYZ norm;

	/* Build up normals incrementally for each vertex of each triangle. */

	for (i = 0, tp = Triangles; i < Ntri; i++, tp++) {
		a = tp->vp[0];
		b = tp->vp[1];
		c = tp->vp[2];

		x = triareaH(a->v, b->v, c->v);
		if (1e8 * x < 50.) {
			continue;
		}

		/* Calculate the normal for this triangle. */

		calc_normal(norm, a->v, b->v, c->v);

		/* Add it to the normals for each of the vertices. */

		vadd(norm, a->n, a->n);
		vadd(norm, b->n, b->n);
		vadd(norm, c->n, c->n);
	}

	/* Normalize all the normals. */

	j = 0;
	for (i = 0, a = Verts; i < Nvert; i++, a++) {
		if (vlength(a->n) == 0.) {
			j++;
		} else {
			vnormalize(a->n);
		}
	}
}

void inflate_mesh(float frac)
{
	int i, j;
	float *v, *n;
	VERTEX *vp;
	XYZ d;

	/* Move each vertex a distance frac along its normal. */

	for (i = 0, vp = Verts; i < Nvert; i++, vp++) {
		v = vp->v;
		n = vp->n;
		vcopy(n, d);
		vscale(d, frac);
		vadd(v, d, v);
	}
}

void write_mesh(void)
{
	int i, j, id;
	float *v, *n;
	VERTEX *vp;
	TRIANGLE *tp;
	XYZ min, max;

	/* Write the verts and normals. */

	printf("%d\n", Nvert);
	for (i = 0, vp = Verts; i < Nvert; i++, vp++) {
		v = vp->v;
		n = vp->n;
		printf("%f %f %f %f %f %f\n", v[0], v[1], v[2], n[0], n[1], n[2]);
	}
	printf("%d\n", Ntri);
	for (i = 0, tp = Triangles; i < Ntri; i++, tp++) {
		for (j = 0; j < 3; j++) {
			vp = tp->vp[j];
			id = vp - Verts;
			printf("%d%c", id, j == 2 ? '\n' : ' ');
		}
	}

	/* Report the bounding box. */

	vset(min, 1e30, 1e30, 1e30);
	vset(max, -1e30, -1e30, -1e30);
	for (i = 0, vp = Verts; i < Nvert; i++, vp++) {
		v = vp->v;
		if (v[0] > max[0]) max[0] = v[0];
		if (v[1] > max[1]) max[1] = v[1];
		if (v[2] > max[2]) max[2] = v[2];
		if (v[0] < min[0]) min[0] = v[0];
		if (v[1] < min[1]) min[1] = v[1];
		if (v[2] < min[2]) min[2] = v[2];
	}
	msg("bb: %g %g, %g %g, %g %g\n", min[0], max[0], min[1], max[1], min[2], max[2]);
}

static inline float triareaH(float *v1, float *v2, float *v3)
{
	float a, b, c, p;
	XYZ d;

	/* Heron's formula. */

	vsub(v1, v2, d);
	a = vlength(d);
	vsub(v2, v3, d);
	b = vlength(d);
	vsub(v3, v1, d);
	c = vlength(d);

	p = (a + b + c) / 2.;
	return sqrt(p * (p - a) * (p - b) * (p - c));
}
