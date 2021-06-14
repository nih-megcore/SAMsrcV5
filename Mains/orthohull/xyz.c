/* Routines for the manipulation of 3-vectors. */

#ifndef _XYZ_C
#define _XYZ_C

#define INLINE inline

INLINE float *vnew(float x, float y, float z)
{
	float *v;

	v = new_array(float, 3);
	v[0] = x;
	v[1] = y;
	v[2] = z;

	return v;
}

INLINE void vset(float *v, float x, float y, float z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

INLINE void vcopy(const float *vsrc, float *vdst)
{
	vdst[0] = vsrc[0];
	vdst[1] = vsrc[1];
	vdst[2] = vsrc[2];
}

INLINE void vprint(const float *v)
{
	printf("%g %g %g\n", v[0], v[1], v[2]);
}

INLINE void vfprint(FILE *f, const float *v)
{
	fprintf(f, "%g %g %g", v[0], v[1], v[2]);
}

INLINE void vzero(float *v)
{
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;
}

INLINE float vlength(const float *v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

INLINE void vscale(float *v, float f)
{
	v[0] *= f;
	v[1] *= f;
	v[2] *= f;
}

INLINE void vnormalize(float *v)
{
#if defined(_IEEE) || defined(_IEEE_FP)
	vscale(v, 1.0 / vlength(v));
#else
	float l;

	l = vlength(v);
	if (l > EPS) {
		vscale(v, 1.0 / l);
	}
#endif
}

INLINE void vmult(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] * src2[0];
	dst[1] = src1[1] * src2[1];
	dst[2] = src1[2] * src2[2];
}

INLINE void vadd(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] + src2[0];
	dst[1] = src1[1] + src2[1];
	dst[2] = src1[2] + src2[2];
}

INLINE void vsub(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] - src2[0];
	dst[1] = src1[1] - src2[1];
	dst[2] = src1[2] - src2[2];
}

INLINE float vdot(const float *v1, const float *v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

INLINE void vcross(const float *v1, const float *v2, float *cross)
{
	XYZ temp;

	  /* x        y       z         z       y */
	temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
	temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
	temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
	vcopy(temp, cross);
}

INLINE void vreflect(const float *in, const float *mirror, float *out)
{
	XYZ temp;

	vcopy(mirror, temp);
	vscale(temp, vdot(mirror, in));
	vsub(temp, in, out);
	vadd(temp, out, out);
}

/* Transform a vector (treated as a column vector) by post-multiplication
with a (3x3) matrix. */

INLINE void vtransform(/*const*/ float * const *m, float *v)
{
	int i;
	float s;
	XYZ temp;

	for (i = 0; i < 3; i++) {
		s  = m[i][0] * v[0];
		s += m[i][1] * v[1];
		s += m[i][2] * v[2];
		temp[i] = s;
	}
	vcopy(temp, v);
}

/* Transform a vector (treated as a column vector) by post-multiplication
with a (4x4) matrix, where the fourth element of the vector is assumed to be
equal to 1. */

INLINE void vtransform4(/*const*/ float * const *m, float *v)
{
	int i;
	float s;
	XYZ temp;

	for (i = 0; i < 3; i++) {
		s  = m[i][0] * v[0];
		s += m[i][1] * v[1];
		s += m[i][2] * v[2];
		s += m[i][3];
		temp[i] = s;
	}
	vcopy(temp, v);
}

INLINE void calc_normal(float *norm,
	const float *a, const float *b, const float *c)
{
	XYZ s, t;

	/* Generate a normal for a triangle.

		   c
		  / \           n = b-a x c-a
		 /   \
		a-----b
	*/

	vsub(b, a, s);
	vsub(c, a, t);
	vcross(s, t, norm);
	vnormalize(norm);
}

#undef INLINE
#endif /* XYZ_C */
