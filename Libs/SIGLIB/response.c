// response() - evaluate response of IIR coefficient
//
//	Author:	Stephen E. Robinson
//			Neuromagnetism Laboratory
//			Dept. of Neurology
//			Henry Ford Hospital
//			Detroit, MI
//			Copyright (c) 2007
//

#include <math.h>
#include <complex.h>
#include <filters.h>

#define TWOPI	(2.* M_PI)
#define TRUE	1
#define FALSE	0


double	response(
	IIRSPEC	*Filter,		// filter coefficients
	double	SampleRate,		// sample rate
	double	fc				// evaluation frequency
)

{
	double _Complex		wp;
	double _Complex		numerator;
	double _Complex		denominator;
	double				omega;
	double				rads;
	int					i;

	if(Filter->enable == TRUE) {
		omega = TWOPI * fc / SampleRate;
		Filter->den[0] = 0.;
		for(i=0, numerator=0., denominator=1.; i<Filter->NC; i++) {
			rads = omega * (double)i;
			wp = -1. * (cos(rads) + I * sin(rads));
			numerator += Filter->num[i] * wp;
			denominator -= Filter->den[i] * wp;
		}
		return(cabs(numerator / denominator));
	} else {
		return(1.);
	}
}
