// mode() -- compute mode of designated time segment
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <math.h>
#include <stdlib.h>


#define FALSE	0
#define TRUE	1


typedef struct {
	double	hi;				// histogram bin upper limit
	double	lo;				// histogram bin lower limit
	double	sum;			// sum of bin member sample values
	int		n;				// number of members
} HISTOBIN;

double	mode(
	double	*x,				// x[TT] -- signal from time segment
	int		TS,				// starting sample
	int		TE				// ending sample
)

{
	static HISTOBIN	*hist;			// hist[NB] -- histogram structure
	double			max;			// maximum value
	double			min;			// minimum value
	double			step;			// bin step-size
	int				mode;			// mode count
	int				mtest;			// mode test count
	int				B;				// current number of bins
	int				b;				// bin index
	int				bm;				// modal bin
	int				t;				// time sample index
	int				done;			// single mode flag
	static int		NB = 0;			// number of bins in histogram
	static int		TT = 0;			// number of samples in time segment

	// if req'd, initialize histogram
	b = TE - TS + 1;
	if(b < 10)
		return HUGE_VAL;
	if(b != TT) {
		if(TT != 0)
			free((void *)hist);
		TT = b;
		NB = TT / 3;
		if((hist = (HISTOBIN *)malloc((size_t)(NB * sizeof(HISTOBIN)))) == NULL)
			return HUGE_VAL;
	}

	// compute peak-to-peak of time segment
	for(t=TS, max=(-1.), min=1.; t<=TE; t++) {
		if(x[t] > max)
			max = x[t];
		if(x[t] < min)
			min = x[t];
	}

	// initialize bin count, & do until there is only one mode (or mean)
	B = NB;
	do {

		// divide histogram into steps, & initialize
		step = fabs(max - min) / (double)(B + 1);
		for(b=0; b<B; b++) {
			hist[b].lo = min + (double)b * step;
			hist[b].hi = min + (double)(b + 1) * step;
			hist[b].sum = 0.;
			hist[b].n = 0;
		}
		hist[B-1].hi += 1.0e-12;	// make the last bin larger, so that sort does not miss maximum sample

		// sort segment amplitude values into bins
		for(t=TS; t<=TE; t++)
			for(b=0; b<B; b++)
				if(x[t] >= hist[b].lo && x[t] < hist[b].hi) {
					hist[b].sum += x[t];
					hist[b].n++;
				}

		// search for modal bin
		for(b=mode=0; b<B; b++)
			if(hist[b].n > mode) {
				mode = hist[b].n;
				mtest = mode - 3;
				bm = b;
			}

		// test to see if the mode is at least 3 counts larger than next mode
		for(b=0, done=TRUE; b<B; b++)
			if(b != bm)				// don't check the mode!
				if(hist[b].n > mtest)
					done = FALSE;	// set flag for next pass

		// set up cond'x for next pass
		if(done == FALSE || B > 1)
			B--;
		else
			done = TRUE;

	} while(done == FALSE);

	// compute average of modal bin, & return value
	if(hist[bm].n > 1)
		return hist[bm].sum / (double)hist[bm].n;
	else
		return 0.;
}
