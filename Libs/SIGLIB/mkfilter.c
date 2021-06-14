// mkfilter() - make an IIR filter from frequency specifications
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdio.h>
#include <math.h>
#include <filters.h>
#include <siglib.h>
#include <lfulib.h>

#define TRUE	1
#define FALSE	0


void	mkfilter(
	IIRSPEC	*Filter,		// pointer to IIR filter structure
	int		Order,			// filter order
	double	SampleRate,		// sample rate (Hz)
	double	HPFreq,			// highpass frequency (Hz)
	double	LPFreq,			// lowpass frequency (Hz)
	double	*BWFreq			// effective bandwidth (Hz)
) {

	double	Freq;			// evaluation frequency
	double	h1;				// transfer function value
	double	gain;			// ratio of filter bandwidth to nominal bandwidth

	Filter->enable = FALSE;
	if(HPFreq > 0. && LPFreq > 0.) {
		if(LPFreq < HPFreq || LPFreq > SampleRate/2.)
			cleanup("bad highpass & lowpass parameters");
		Filter->order = Order;				// set filter order
		Filter->type = BANDPASS;			// bandpass filter selection
		Filter->fl = HPFreq;				// low edge frequency
		Filter->fh = LPFreq;				// high edge frequency
		Filter->rate = SampleRate;			// sample rate
		Filter->enable = TRUE;
		mkiir(Filter);
	} else if(HPFreq > 0. && LPFreq == 0.) {
		if(HPFreq > SampleRate/2.)
			cleanup("bad highpass parameter");
		Filter->order = 2 * Order;			// set filter order
		Filter->type = HIGHPASS;			// highpass filter selection
		Filter->fl = HPFreq;				// low edge frequency
		Filter->fh = SampleRate / 2.;		// high edge frequency
		Filter->rate = SampleRate;			// sample rate
		Filter->enable = TRUE;
		mkiir(Filter);
	} else if(HPFreq == 0. && LPFreq > 0.) {
		if(LPFreq >= SampleRate/2.)
			cleanup("bad lowpass parameter");
		Filter->order = 2 * Order;			// set filter order
		Filter->type = LOWPASS;				// lowpass filter selection
		Filter->fl = 0.;					// low edge frequency
		Filter->fh = LPFreq;				// high edge frequency
		Filter->rate = SampleRate;			// sample rate
		Filter->enable = TRUE;
		mkiir(Filter);
	}

	// if mkiir succeeded
	if(Filter->enable == TRUE) {

		// compute effective bandwidth of data filter in steps of 1 milliHz
		for(Freq=*BWFreq=0.; Freq<(SampleRate/4.); Freq+=0.001) {
			h1 = response(Filter, SampleRate, Freq);
			if(h1 > M_SQRT2) {
				printf("\07\nERROR: IIR filter is unstable for this sample rate\n");
				printf("\tgain=%f at %f Hz\n", h1*h1, Freq);
				printf("\tfrequency parameters are too close to limits\n");
				cleanup("data should be decimated to use these filter parameters");
			}		
			*BWFreq += 0.001 * (h1 * h1);
		}

		// test filter for low gain
		gain = *BWFreq / (LPFreq - HPFreq);
		if(gain < M_SQRT1_2) {
			printf("\07\nWARNING: IIR filter gain too low\n");
			printf("\tnominal filter bandwidth=%f\n", *BWFreq);
			cleanup("suggest increasing bandwidth of filter parameters");
		}
	} else {
		*BWFreq = 0.;
	}
}
