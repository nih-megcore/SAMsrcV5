// detrend() - remove trend (regression) line from time-series
//
//	Author:	Stephen E. Robinson
//			Neuromagnetism Laboratory
//			Dept. of Neurology
//			Henry Ford Hospital
//			Detroit, MI
//			Copyright (c) 2007
//

#include <stdlib.h>

void	detrend(
	double	*y,			// x time-series
	int		T			// number of samples
)

{
	double			sx;			// sum of x
	double			sx2;		// sum of x^2
	double			sy;			// sum of y
	double			sxy;		// sum of x*y
	double			Sxx;		// sum-of-squares of x
	double			Sxy;		// sum-of-products
	double			mx;			// mean x
	double			my;			// mean y
	double			slope;		// trend slope
	double			intercept;	// trend y-intercept
	double			nt;			// number of samples
	register int	t;			// time index

	// accumulate statistics & compute sums of squares & product
	nt = (double)T;
	for(t=0, sx=sx2=sy=sxy=0.; t<T; t++) {
		sx += (double)t;
		sx2 += (double)(t * t);
		sy += y[t];
		sxy += (double)t * y[t];
	}
	Sxx = sx2 - sx * sx / nt;	// sum of squares of x
	Sxy = sxy - sx * sy / nt;	// sum of products
	mx = sx / nt;
	my = sy / nt;

	// compute slope & intercept
	slope = Sxy / Sxx;			// slope
	intercept = my - slope * mx;		// intercept

	/* remove trend line */
	for(t=0; t<T; t++)
		y[t] -= slope * (double)t + intercept;
}
