// spectrum() -- computes the spectrum from 0 to 1/4 sample rate

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <siglib.h>

#define REAL_	0
#define IMAG_	1


void	spectrum(
	double	*x,		// x[T] -- input time-series
	double	*s,		// s[F] -- spectrum
	double	srate,	// sample rate
	int		N		// number of samples
)

{
	static fftw_plan	fft_x;		// fft plan for time-series x
	static fftw_complex	*X;			// X[N] -- spectrum of time-series x
	static double		*In;		// In[N] -- windowed input to FFT
	static double		*Window;	// Window[N] -- Hanning window coefficients
	int					f;			// frequency bin index
	int					t;			// sample index
	static int			T = 0;

	// one-time initialization
	if(T == 0) {
		X = fftw_malloc(sizeof(fftw_complex) * N);
		if((In = (double *)malloc((size_t)N * sizeof(double))) == NULL)
			allocfailed("In[]");
		if((Window = (double *)malloc((size_t)N * sizeof(double))) == NULL)
			allocfailed("Window[]");
		Hanning(Window, N);
		fft_x = fftw_plan_dft_r2c_1d(N, In, X, FFTW_MEASURE);
		T = N;
	}

	// apply hanning window & execute FFT
	for(t=0; t<T; t++)
		In[t] = Window[t] * x[t];
	fftw_execute(fft_x);

	// compute spectrum
	for(f=0; f<N/2; f++)
		s[f] = sqrt((double)X[f][REAL_] * (double)X[f][REAL_] + (double)X[f][IMAG_] * (double)X[f][IMAG_]);
}
