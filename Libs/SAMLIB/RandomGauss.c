// RandomGauss() -- returns a Gaussian random number with
//	mean 'm' and standard deviation 's'.  Routine is
//	initialized with the time as a seed.
//
//	Author:	Stephen E. Robinson (after Sarvas & Hamalainen)
//			Neuromagnetism Laboratory
//			Dept.of Neurology
//			Henry Ford Hospital
//			Detroit, MI
//			Copyright (c) 2007

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

double	RandomGauss(
	double	m,			// mean
	double	s			// sigma
) {
	double			r;			// vector magnitude
	static double	v1;			// uniform number 1 (-1.0 to +1.0)
	static double	v2;			// uniform number 2 (-1.0 to +1.0)
	static double	BoxMuller;	// Box-Muller factor
	static int		mem = 0;	// initializer memory
	static int		pair = 0;	// odd-even flag for pairs of Gaussians

	// on first call, initialize seed using system time
	if(mem == 0) {
		srand48((long)time((time_t *)0));
		mem++;
	}

	// now compute Gaussian distributed noise
	if(pair == 0) {
		pair = 1;
		do {
			v1 = 2. * drand48() - 1.;
			v2 = 2. * drand48() - 1.;
			r = v1 * v1 + v2 * v2;
		} while(r >= 1.);
		BoxMuller = sqrt(-2. * log(r) / r);
		return v2 * BoxMuller * s + m;
	} else {
		pair = 0;
		return v1 * BoxMuller * s + m;
	}
}
