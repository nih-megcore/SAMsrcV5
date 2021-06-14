// FFTfilter() -- filter data using gain coefficients. A forward
//	fft is applied to the data. The complex results are multiplied
//	by the gain coefficients. The complex is then transformed
//	back to real using an inverse fft. Note that the complex
//	format of the forward fft has the positive frequencies in
//	0 to len/2 and the negative frequencies in len/2+1 to len-1.
//	negative frequencies are stored in reverse (mirrored) order.
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <siglib.h>


void	FFTfilter(
	double		*data,		// data[N] -- input data
	double		*result,	// result[N] -- output data (filtered)
	FFTSPEC		*Filter,	// filter structure
	int			len			// number of samples
)
{
	static fftw_plan	p1;				// forward fft plan
	static fftw_plan	p2;				// reverse fft plan
	static fftw_complex *x;				// real data -- x[i][0]
	static fftw_complex	*X;				// complex data
	double				*p;				// result pointer
	int					i;				// i-index (positive freqs)
	int					j;				// j-index (negative freqs)
	static int			planlen = 0;	// initialization
	char				*Wisfile = NULL; // wisdom file string
	FILE				*wisdom;		// wisdom file pointer

	// sanity check -- len must match filter
	if(len != Filter->N)
		cleanup("filter length mismatch");

	// Keep the arrays and plans around from last time, since this
	// is a very common case. Reallocate them if they change.
	if (len != planlen && planlen > 0) {
		fftw_destroy_plan(p1);
		fftw_destroy_plan(p2);
		fftw_free(x);
		fftw_free(X);
		planlen = 0;
	}
	if (planlen == 0) {
		planlen = len;
		x = fftw_malloc(sizeof(fftw_complex) * len);
		X = fftw_malloc(sizeof(fftw_complex) * len);

		// Get any accumulated wisdom.
		set_wisfile();
		wisdom = fopen(Wisfile, "r");
		if(wisdom) {
			fftw_import_wisdom_from_file(wisdom);
			fclose(wisdom);
		}

		// Set up the fftw plans.
		p1 = fftw_plan_dft_1d(len, x, X, FFTW_FORWARD, FFTW_ESTIMATE);
		p2 = fftw_plan_dft_1d(len, X, x, FFTW_BACKWARD, FFTW_ESTIMATE);

		// Save the wisdom.
		wisdom = fopen(Wisfile, "w");
		if(wisdom) {
			fftw_export_wisdom_to_file(wisdom);
			fclose(wisdom);
		}
	}

	// Convert the input to complex.
	memset(x, 0, sizeof(fftw_complex) * len);
	for(i=0; i<len; i++)
		x[i][0] = data[i];

	// FFT.
	fftw_execute(p1);	// x -> X

	// multiply spectral coefficients by filter gain -- careful to include positive & negative frequencies!
	for(i=0, j=len-1; i<len/2; i++, j--) {
		X[i][0] *= Filter->Gain[i];		// positive frequency real
		X[i][1] *= Filter->Gain[i];		// positive frequency imaginary
		X[j][0] *= Filter->Gain[j];		// negative frequency real
		X[j][1] *= Filter->Gain[j];		// negative frequency imaginary
	}

	// Inverse FFT.
	fftw_execute(p2);	// X -> x

	// Fill in the rows of the result.
	for(i=0, p=result; i<len; i++)
		*p++ = x[i][0] / (double)len;
}
