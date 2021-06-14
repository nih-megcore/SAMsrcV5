// bdiir() - zero phase-shift infinite impulse response digital filter
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdlib.h>
#include "filters.h"
#include <lfulib.h>

void	bdiir(
	double	*In,				// input data array
	double	*Out,				// output data array
	int		T,					// number of points in data array
	IIRSPEC	*Filter				// iir filter structure pointer
)

{
	static double	*Tmp;		// temporary time-series
	register double x;			// input sample value;
	register double y;			// output sample value;
	register int	t;			// time-index
	register int	i;			// coefficient-index
	register int	j;			// offset index
	static int		mem = 0;

	// if req'd, allocate or modify intermediate array
	if(mem != T) {
		if(mem != 0)
			free((void *)Tmp);
		if((Tmp = (double *)malloc((size_t)(T * sizeof(double)))) == NULL)
			allocfailed("Tmp[]");
		mem = T;
	}

	// filter is not enabled, move data to output
	if(Filter->enable == 0) {
		for(t=0; t<T; t++)
			Out[t] = In[t];
		return;
	}
	Filter->den[0] = 0.;		// set 1st term of denominator to zero
	for(t=0; t<T; t++)
		for(i=0, Tmp[t]=0.; i<Filter->NC; i++) {
			j = t - i;
			if(j < 0) {
				x = y = 0.;
			} else {
				x = In[j];
				y = Tmp[j];
			}
			Tmp[t] += Filter->num[i] * x - Filter->den[i] * y;
		}
	for(t=(T-1); t>=0; t--)
		for(i=0, Out[t]=0.; i<Filter->NC; i++) {
			j = t + i;
			if(j >= T) {
				x = y = 0.;
			} else {
				x = Tmp[j];
				y = Out[j];
			}
			Out[t] += Filter->num[i] * x - Filter->den[i] * y;
		}
}
