// power() -- compute power of time series
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <math.h>
#include <lfulib.h>

double	power(
	double	*x,				// x[t] -- time-series
	int		TS,				// starting time index
	int		TE				// ending time index
) {
	double			ts;		// total number of samples
	double			sx;		// sum of x
	double			sx2;	// sum of x^2
	double			s2;		// power
	int				t;		// sample index

	// check size of time segment
	ts = (double)(TE - TS + 1);
	if(ts < 12.)
		cleanup("power() has too few samples");

	// accumulate sums
	for(t=TS, sx=sx2=0.; t<=TE; t++) {
		sx += x[t];
		sx2 += x[t] * x[t];
	}

	// compute segment mean-square power
	s2 = (sx2 - ((sx * sx) / ts)) / ts;
	
	return s2;
}

