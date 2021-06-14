// demean() - trivial routine to remove mean value from time-series
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdlib.h>

void	demean(
	double	*x,			// x time-series
	int	T				// number of samples
)

{
	double	nt;		// number of points
	double	sx;		// sum of x
	double	mx;		// mean x
	int		t;		// time index

	// compute mean of x(t)
	nt = (double)T;
	for(t=0, sx=0.; t<T; t++)
		sx += x[t];
	mx = sx / nt;

	// remove mean from x(t)
	for(t=0; t<T; t++)
		x[t] -= mx;
}
