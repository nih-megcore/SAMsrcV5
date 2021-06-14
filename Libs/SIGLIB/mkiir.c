// mkiir() - compute Butterworth IIR filter coefficients
//
//	Butterworth digital filter design program.
//
//	Method used is to generate an analog filter and
//	then use Bilinear transform to map to a digital
//	filter.  Filters up to 20th order may be designed.
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <math.h>
#include <complex.h>
#include <filters.h>
#include <lfulib.h>

#define P           2*(MAXORDER+1)

void	mkiir(
    IIRSPEC *fp     // pointer to iir filter
)

{
    double _Complex znpoly[P];
    double _Complex zdpoly[P];
    double _Complex numrut[P];
    double _Complex denrut[P];
    double _Complex zpole[P];
    double _Complex zzero[P];
    double _Complex root[P];
    double _Complex cgain;
    double _Complex btmp;
    double _Complex dtmp;
    double 			psi;
    double 			arg;
    double 			gain = 1.;
    double 			wo;
    double 			wh;
    double 			wc;
    double 			tmp;
    int 			order;
    int 			even;
    int 			npole;
    int 			nzero = 0;
    register int 	i;
    register int 	j;

    // set number of coefficients
    switch (fp->type) {
	    case HIGHPASS:
	    case LOWPASS:
			fp->NC = fp->order + 1;
			break;
	    case BANDPASS:
	    case BANDREJECT:
			fp->NC = 2 * fp->order + 1;
			break;
	    default:
			cleanup("unknown iir filter type");
    }

    // initialize polynomials
    for (i = 1; i < P; i++)
		znpoly[i] = zdpoly[i] = (0. + 0.0Fi);

    // compute location of Butterworth poles from filter order
    even = (fp->order / 2) * 2;
    psi = M_PI / (double)fp->order;
    arg = (even == fp->order) ? .5 * psi : psi;
    for (i = 1; i <= even; i += 2) {
		denrut[i] = -1. * (cos(arg) + I * sin(arg));
		denrut[i + 1] = conj(denrut[i]);
		arg += psi;
    }
    if (even != fp->order)
		denrut[fp->order] = (-1. + 0.0Fi);
    npole = fp->order;

    // map frequency to warped z-plane frequency
    tmp = M_PI / fp->rate;      // conversion factor for frequency to radians/sec
    wo = tan((double)(fp->fl * tmp));
    wh = tan((double)(fp->fh * tmp));

    // compute location of Butterworth zeroes from filter type
    switch (fp->type) {
	    case LOWPASS:
			nzero = 0;
			for (i = 1, gain = 1.; i <= npole; i++) {
			    denrut[i] = wh * denrut[i];
			    gain *= wh;         // gain = wh^npole
			}
			break;

	    case HIGHPASS:
			nzero = npole;
			for (i = 1; i <= nzero; i++) {
			    numrut[i] = (0. + 0.0Fi);
			    denrut[i] = wo / denrut[i];
			}
			gain = 1.;
			break;

	    case BANDPASS:
			nzero = npole;
			wc = wo * wh;
			wh -= wo;
			for (i = 1, gain = 1.; i <= npole; i++) {
			    numrut[i] = (0. + 0.0Fi);
			    gain *= wh;         // gain = wh^npole
			}
			for (i = 1, j = 1; i <= npole; i++, j += 2) {
			    btmp = -wh * denrut[i];
			    dtmp = csqrt(btmp * btmp - 4. * wc);
			    root[j] = 0.5 * (dtmp - btmp);
			    root[j + 1] = 0.5 * (-dtmp - btmp);
			}
			npole *= 2;
			for (i = 1; i <= npole; i++)
			    denrut[i] = root[i];
			break;

	    case BANDREJECT:
			gain = 1.;
			wc = wo * wh;
			wh -= wo;
			for (i = 1, j = 1; i <= npole; i++, j += 2) {
			    numrut[j] = (0. + sqrt(wc) * 1.0Fi);
			    numrut[j + 1] = conj(numrut[j]);
			    btmp = csqrt(wh * wh - 4. * wc * denrut[i] * denrut[i]);
			    dtmp = 2. * denrut[i];
			    root[j] = (wh + btmp) / dtmp;
			    root[j + 1] = (wh - btmp) / dtmp;
			}
			npole *= 2;
			nzero = npole;
			for (i = 1; i <= npole; i++)
			    denrut[i] = root[i];
			break;

	    default:
			cleanup("unknown iir filter type");
    }

    // map poles & zeroes from s to z-plane
    cgain = (gain + 0.0Fi);
    for (i = 1; i <= nzero; i++) {
		zpole[i] = (1. + denrut[i]) / (1. - denrut[i]);
		zzero[i] = (1. + numrut[i]) / (1. - numrut[i]);
		cgain *= (1. - numrut[i]) / (1. - denrut[i]);
    }

    // map any zeroes at +/- infinity to (-1, 0)
    if (npole > nzero) {
		for (i = (nzero + 1); i <= npole; i++) {
		    zzero[i] = (-1. + 0.0Fi);
		    zpole[i] = (1. + denrut[i]) / (1. - denrut[i]);
		    cgain /= 1. - denrut[i];
		}
		nzero = npole;
    }
    cgain *= conj(cgain);
    gain = sqrt(creal(cgain));

    // multiply out all sections to form numerator & denominator polynomial
    znpoly[1] = (1. + 0.0Fi);
    zdpoly[1] = (1. + 0.0Fi);
    znpoly[2] = -zzero[1];
    zdpoly[2] = -zpole[1];
    for (j = 2, order = 1; j <= npole; j++, order++) {
		for (i = (order + 2); i >= 2; i--) {
		    znpoly[i] -= zzero[j] * znpoly[i - 1];
		    zdpoly[i] -= zpole[j] * zdpoly[i - 1];
		}
    }

    // multiply numerator by gain
    for (i = 1; i <= (nzero + 1); i++)
		znpoly[i] = gain * znpoly[i];

    // output real part of complex polynomial as numerator & denominator coefficients
    for (i = 1; i <= (npole + 1); i++) {
		fp->num[i - 1] = creal(znpoly[i]);
		fp->den[i - 1] = creal(zdpoly[i]);
    }
}

