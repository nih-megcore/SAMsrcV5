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
static float tvol(VERTEX *p1, VERTEX *p2, VERTEX *p3);
static float meshvol(void);
void run_qhull(void);

int iflag = FALSE;      /* inflate the mesh by frac */
float frac;

int qflag = FALSE;      /* convexify the hull with qhull */

char Usage[] = "usage: %s [-i<inflate_frac>] [-q] [file.asc]\n";
#define USAGE() msg(Usage, Progname)

int main(int argc, char **argv)
{
	int c;
	FILE *infile;

	Progname = *argv;

	/* Parse option arguments. */

	while ((c = getopt(argc, argv, "i:q")) != EOF) {
		switch (c) {

		case 'i':
			frac = atof(optarg);
			iflag = TRUE;
			break;

		case 'q':
			qflag = TRUE;
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

	if (qflag) {
		run_qhull();
	}

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

	/* This reads .asc files. Parse the header. */

	do {
		fgets(buf, 100, infile);
	} while (buf[0] == '#');
	sscanf(buf, "%d %d", &Nvert, &Ntri);

	/* Read the verts. */

	Verts = new_array(VERTEX, Nvert);
	for (i = 0, vp = Verts; i < Nvert; i++, vp++) {
		v = vp->v;
		vp->vlist = NULL;
		vzero(vp->n);                   /* init the normal */
		fgets(buf, 100, infile);
		sscanf(buf, "%f %f %f %*f", v, v + 1, v + 2);

		/* Convert from LPI mm to PRI m. */

		x = v[0];
		v[0] = v[1];
		v[1] = -x;
		vscale(v, .001);
	}

	/* Next is the triangles. Keep track of the edges. */

	Triangles = new_array(TRIANGLE, Ntri);
	for (i = 0, tp = Triangles; i < Ntri; i++, tp++) {
		fgets(buf, 100, infile);
		sscanf(buf, "%d %d %d %*d", &a, &b, &c);
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

	/* Now write the triangles. */

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
	msg("bb: %g %g, %g %g, %g %g cm\n",
		min[0] * 100., max[0] * 100.,
		min[1] * 100., max[1] * 100.,
		min[2] * 100., max[2] * 100.);

	/* Report the volume. */

	msg("vol: %g ml\n", meshvol() * 1e6);
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

static float tvol(VERTEX *p1, VERTEX *p2, VERTEX *p3)
{
	float v321 = p3->v[0] * p2->v[1] * p1->v[2];
	float v231 = p2->v[0] * p3->v[1] * p1->v[2];
	float v312 = p3->v[0] * p1->v[1] * p2->v[2];
	float v132 = p1->v[0] * p3->v[1] * p2->v[2];
	float v213 = p2->v[0] * p1->v[1] * p3->v[2];
	float v123 = p1->v[0] * p2->v[1] * p3->v[2];

	return (v123 + v231 + v312 - v132 - v213 - v321) / 6.;
}

static float meshvol(void)
{
	int i;
	float s;
	TRIANGLE *tp;

	s = 0.;
	for (i = 0, tp = Triangles; i < Ntri; i++, tp++) {
		s += tvol(tp->vp[0], tp->vp[1], tp->vp[2]);
	}

	return s;
}

void run_qhull(void)
{
	int i, j, ntri, *triangles, *t, *map, *imap, a, b, c;
	char *s, *cmd, buf[100], pointname[100], hullname[100];
	FILE *infile, *outfile, *hullfile;
	float *v;
	VERTEX *vp;
	TRIANGLE *tp;

	/* Output the vertices in the format qhull likes. */

	sprintf(pointname, "/tmp/points.%d", getpid());
	outfile = fopen(pointname, "w");
	if (outfile == NULL) {
		fatalerr("can't write %s", pointname);
	}
	fprintf(outfile, "3\n%d\n", Nvert);
	for (i = 0, vp = Verts; i < Nvert; i++, vp++) {
		v = vp->v;
		fprintf(outfile, "%g %g %g\n", v[0], v[1], v[2]);
	}
	fclose(outfile);

	/* Run the qhull command with the right arguments. */

	sprintf(hullname, "/tmp/hull.%d", getpid());
	sprintf(buf, "qhull QJ i < %s > %s", pointname, hullname);
	msg("computing hull\n");
	system(buf);

	/* Now read the hull in. */

	hullfile = fopen(hullname, "r");
	if (hullfile == NULL) {
		fatalerr("can't read %s", hullname);
	}
	fgets(buf, 100, hullfile);
	ntri = atoi(buf);
	triangles = t = new_array(int, ntri * 3);
	for (i = 0; i < ntri; i++, t += 3) {
		fgets(buf, 100, hullfile);
		sscanf(buf, "%d %d %d", t, t + 1, t + 2);
	}
	fclose(hullfile);

	/* The hull has fewer points than the original vertex set.
	We want to renumber the vertices that are actually used so
	we can output minimal geometry. Make a map from old vertex
	numbers to new. We'll run through the hull, and when we find
	a new vertex we'll map it to a new number. */

	map = new_array(int, Nvert);
	imap = new_array(int, Nvert);
	for (i = 0; i < Nvert; i++) {
		map[i] = -1;
	}
	j = 0;
	for (i = ntri * 3, t = triangles; i > 0; i--, t++) {
		if (map[*t] < 0) {
			map[*t] = j;
			imap[j++] = *t;
		}
	}

	/* Now we can recreate the geometry for the hull, containing
	only the points that are actually used. Invert the order so the
	normals point outwards. */

	Nvert = j;
	vp = new_array(VERTEX, Nvert);
	for (i = 0; i < Nvert; i++) {
		vcopy(Verts[imap[i]].v, vp[i].v); /* V -> v */
		vp[i].vlist = NULL;
		vzero(vp->n);                   /* init the normal */
	}

	free(Verts);
	free(Triangles);
	Verts = vp;

	/* Next is the triangles. Keep track of the edges. */

	Ntri = ntri;
	Triangles = new_array(TRIANGLE, Ntri);
	for (i = 0, tp = Triangles, t = triangles; i < Ntri; i++, tp++, t += 3) {
		a = map[t[2]];
		b = map[t[1]];
		c = map[t[0]];
		tp->vp[0] = &Verts[a];
		tp->vp[1] = &Verts[b];
		tp->vp[2] = &Verts[c];
		make_edge(&Verts[a], &Verts[b]);
		make_edge(&Verts[b], &Verts[c]);
		make_edge(&Verts[c], &Verts[a]);
	}

	/* Clean up. */

	unlink(pointname);
	unlink(hullname);
	free(map);
	free(imap);
}

