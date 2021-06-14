// filters.h
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#ifndef H_FILTERS
#define H_FILTERS

#define MAXORDER	20
#define ENABLE		1
#define DISABLE		0

// filter type definitions
#define LOWPASS		0
#define HIGHPASS	1
#define BANDPASS	2
#define BANDREJECT	3
#define ALLPASS		4

// FFT filter structure
typedef struct fftspec {
	double	*Gain;			// Gain[N] -- filter gain vs frequency
	int		type;			// filter type (from above)
	int		order;			// filter order
	double	rate;			// sample rate (Hz)
	int		N;				// number of samples
} FFTSPEC;

// IIR filter structure
typedef struct	iirspec {
	int		type;			// filter type
	int		enable;			// 1=enable filter, 0=disable filter
	double	fh;				// high edge filter frequency (Hz)
	double	fl;				// low edge filter frequency (Hz)
	double	rate;			// sample rate (Hz)
	double	bwn;			// noise bandwidth (Hz)
	int		order;			// filter order
	int		NC;				// number of coefficient pairs
	double	num[MAXORDER];	// numerator iir coefficients
	double	den[MAXORDER];	// denominator iir coefficients
} IIRSPEC;

// comb filter structure
typedef struct	comb {
	double	LowEdge;
	double	HighEdge;
	int		Order;
} COMB;

#endif	// H_FILTERS
