// Hanning() -- generate a Hanning window function
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <math.h>
#include <lfulib.h>

void	Hanning(
	double	*Window,	// Window[N] -- the window function
	int		N			// the number of points in the window
)
	
{
	register int	n;	// n-index

	if(N > 2) {
		for(n=0; n<N; n++)
			Window[n] = 0.5 * (1. - cos((2. * M_PI * (double)n) / ((double)N - 1.)));
	} else {
		cleanup("Hanning window has less than 3 samples");
	}
}
