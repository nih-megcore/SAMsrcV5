// debase() - remove mean to quartic fit from time-series
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <math.h>
#include <stdlib.h>
#include <samlib.h>
#include <lfulib.h>
#include <siglib.h>

void	debase(
	double	*y,				// x time-series
	int	T,					// number of samples
	int	Order				// baseline order
)

{
	static double	**a;	// a[5][5] -- matrix for inverse solution
	static double	*b;		// b[5] -- solution coefficient vector
	double			nt;		// total number of samples
	double			t1;		// t
	double			t2;		// t^2
	double			t3;		// t^3
	double			t4;		// t^4
	double			t5;		// t^5
	double			t6;		// t^6
	double			t7;		// t^7
	double			t8;		// t^8
	double			sx;		// sum of x
	double			sx2;	// sum of x^2
	double			sx3;	// sum of x^3
	double			sx4;	// sum of x^4
	double			sx5;	// sum of x^5
	double			sx6;	// sum of x^6
	double			sx7;	// sum of x^7
	double			sx8;	// sum of x^8
	double			sy;		// sum of y
	double			syx;	// sum of y*x
	double			syx2;	// sum of y*x^2
	double			syx3;	// sum of y*x^3
	double			syx4;	// sum of y*x^4
	int				i;		// matrix index
	register int	t;		// time index;
	static int		once = FALSE;	// matrix allocation flag

	// one-time allocation for inverse solutions
	if(once == FALSE) {
		if((a = (double **)malloc(5 * sizeof(double *))) == NULL)
			allocfailed("a[]");
		for(i=0; i<5; i++)
			if((a[i] = (double *)malloc(5 * sizeof(double))) == NULL)
				allocfailed("a[][]");
		if((b = (double *)malloc(5 * sizeof(double))) == NULL)
			allocfailed("b[]");
		once = TRUE;
	}
	nt = (double)T;

	// solutions for each order
	switch(Order) {
		case 0:		// remove mean

			// accumulate statistics
			for(t=0, sy=0.; t<T; t++)
				sy += y[t];

			// solve:
			//	sy = b0 * nt
			b[0] = sy / nt;

			// remove mean baseline
			for(t=0; t<T; t++)
				y[t] -= b[0];
			break;

		case 1:	// remove trend

			// accumulate statistics & compute sums of squares & product
			for(t=0, sx=sx2=sy=syx=0.; t<T; t++) {
				t1 = (double)t;
				t2 = t1 * (double)t;
				sx += t1;
				sx2 += t2;
				sy += y[t];
				syx += y[t] * t1;
			}

			// complete solution matrix, solving:
			//	sy  = b0 * nt + b1 * sx
			//	syx = b0 * sx + b1 * sx2
			a[0][0] = nt;	a[0][1] = sx;
			a[1][0] = sx;	a[1][1] = sx2;
			b[0] = sy;
			b[1] = syx;
			LUsolve(a, b, 2);

			// remove trend
			for(t=0; t<T; t++) {
				t1 = (double)t;
				y[t] -= b[0] + b[1] * t1;
			}
			break;

		case 2:		// remove quadratic

			// accumulate statistics
			for(t=0, sx=sx2=sx3=sx4=sy=syx=syx2=0.; t<T; t++) {
				t1 = (double)t;
				t2 = t1 * (double)t;
				t3 = t2 * (double)t;
				t4 = t3 * (double)t;
				sx += t1;
				sx2 += t2;
				sx3 += t3;
				sx4 += t4;
				sy += y[t];
				syx += y[t] * t1;
				syx2 += y[t] * t2;
			}

			// complete solution matrix, solving
			//	sy   = b[0] * nt  + b[1] * sx  + b[2] * sx2
			//	syx  = b[0] * sx  + b[1] * sx2 + b[2] * sx3
			//	syx2 = b[0] * sx2 + b[1] * sx3 + b[2] * sx4
			a[0][0] = nt;	a[0][1] = sx;	a[0][2] = sx2;
			a[1][0] = sx;	a[1][1] = sx2;	a[1][2] = sx3;
			a[2][0] = sx2;	a[2][1] = sx3;	a[2][1] = sx4;
			b[0] = sy;
			b[1] = syx;
			b[2] = syx2;
			LUsolve(a, b, 3);

			// remove quadratic curve
			for(t=0; t<T; t++) {
				t1 = (double)t;
				t2 = t1 * (double)t;
				y[t] -= b[0] + b[1] * t1 + b[2] * t2;
			}
			break;

		case 3:		// remove cubic trend

			// accumulate statistics
			for(t=0, sx=sx2=sx3=sx4=sx5=sx6=sy=syx=syx2=syx3=0.; t<T; t++) {
				t1 = (double)t;
				t2 = t1 * (double)t;
				t3 = t2 * (double)t;
				t4 = t3 * (double)t;
				t5 = t4 * (double)t;
				t6 = t5 * (double)t;
				sx += t1;
				sx2 += t2;
				sx3 += t3;
				sx4 += t4;
				sx5 += t5;
				sx6 += t6;
				sy += y[t];
				syx += y[t] * t1;
				syx2 += y[t] * t2;
				syx3 += y[t] * t3;
			}

			// complete solution matrix, solving:
			//	sy   = b[0] * nt  + b[1] * sx  + b[2] * sx2 + b[3] * sx3
			//	syx  = b[0] * sx  + b[1] * sx2 + b[2] * sx3 + b[3] * sx4
			//	syx2 = b[0] * sx2 + b[1] * sx3 + b[2] * sx4 + b[3] * sx5
			//	syx3 = b[0] * sx3 + b[1] * sx4 + b[2] * sx5 + b[3] * sx6
			a[0][0] = nt;	a[0][1] = sx;	a[0][2] = sx2;	a[0][3] = sx3;
			a[1][0] = sx;	a[1][1] = sx2;	a[1][2] = sx3;	a[1][3] = sx4;
			a[2][0] = sx2;	a[2][1] = sx3;	a[2][2] = sx4;	a[2][3] = sx5;
			a[3][0] = sx3;	a[3][1] = sx4;	a[3][2] = sx5;	a[3][3] = sx6;
			b[0] = sy;
			b[1] = syx;
			b[2] = syx2;
			b[3] = syx3;
			LUsolve(a, b, 4);

			// remove cubic curve
			for(t=0; t<T; t++) {
				t1 = (double)t;
				t2 = t1 * (double)t;
				t3 = t2 * (double)t;
				y[t] -= b[0] + b[1] * t1 + b[2] * t2 + b[3] * t3;
			}
			break;

		case 4:		// remove quatric curve

			// accumulate statistics
			for(t=0, sx=sx2=sx3=sx4=sx5=sx6=sx7=sx8=sy=syx=syx2=syx3=syx4=0.; t<T; t++) {
				t1 = (double)t;
				t2 = t1 * (double)t;
				t3 = t2 * (double)t;
				t4 = t3 * (double)t;
				t5 = t4 * (double)t;
				t6 = t5 * (double)t;
				t7 = t6 * (double)t;
				t8 = t7 * (double)t;
				sx += t1;
				sx2 += t2;
				sx3 += t3;
				sx4 += t4;
				sx5 += t5;
				sx6 += t6;
				sx7 += t7;
				sx8 += t8;
				sy += y[t];
				syx += y[t] * t1;
				syx2 += y[t] * t2;
				syx3 += y[t] * t3;
				syx4 += y[t] * t4;
			}

			// complete solution matrix, solving:
			//	sy   = b[0] * nt  + b[1] * sx  + b[2] * sx2 + b[3] * sx3 + b[4] * sx4
			//	syx  = b[0] * sx  + b[1] * sx2 + b[2] * sx3 + b[3] * sx4 + b[4] * sx5
			//	syx2 = b[0] * sx2 + b[1] * sx3 + b[2] * sx4 + b[3] * sx5 + b[4] * sx6
			//	syx3 = b[0] * sx3 + b[1] * sx4 + b[2] * sx5 + b[3] * sx6 + b[4] * sx7
			//	syxr = b[0] * sx4 + b[1] * sx5 + b[2] * sx6 + b[3] * sx7 + b[4] * sx8
			a[0][0] = nt;	a[0][1] = sx;	a[0][2] = sx2;	a[0][3] = sx3;	a[0][4] = sx4;
			a[1][0] = sx;	a[1][1] = sx2;	a[1][2] = sx3;	a[1][3] = sx4;	a[1][4] = sx5;
			a[2][0] = sx2;	a[2][1] = sx3;	a[2][2] = sx4;	a[2][3] = sx5;	a[2][4] = sx6;
			a[3][0] = sx3;	a[3][1] = sx4;	a[3][2] = sx5;	a[3][3] = sx6;	a[3][4] = sx7;
			a[4][0] = sx4;	a[4][1] = sx5;	a[4][2] = sx6;	a[4][3] = sx7;	a[4][4] = sx8;
			b[0] = sy;
			b[1] = syx;
			b[2] = syx2;
			b[3] = syx3;
			b[4] = syx4;
			LUsolve(a, b, 5);

			// remove quatric curve
			for(t=0; t<T; t++) {
				t1 = (double)t;
				t2 = t1 * (double)t;
				t3 = t2 * (double)t;
				t4 = t3 * (double)t;
				y[t] -= b[0] + b[1] * t1 + b[2] * t2 + b[3] * t3 + b[4] * t4;
			}
			break;

		default:
			cleanup("no such baseline removal order");
	}
}
