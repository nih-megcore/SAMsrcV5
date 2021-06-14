// mknotch() -- make notch filters for designated bandpass
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <filters.h>
#include <siglib.h>
#include "sam_param.h"

// notch filter frequencies
#define NOTCHW_2    1.              // half notch width -- +/- 1 Hz
#define MAXNOTCH    20              // maximum 20 notches

void mknotch(
    FFTSPEC *np,                    // notch filter pointer (output)
    double  HPFreq,                 // highpass frequency (Hz) (input)
    double  LPFreq,                 // lowpass frequency (Hz) (input)
    double  Rate,                   // sample rate (Hz) (input)
    int     N,                      // number of samples
    double  Mains                   // mains frequency
) {
    FFTSPEC         Notch;          // one of the notch filters
    double          fc;             // notch center frequency
    double          fh;             // notch highpass frequency
    double          fl;             // notch lowpass frequency
    double          BW;             // notch effective bandwidth
    int             n;              // notch index
    int             c;              // coefficient index

    // make an allpass filter first!
    Butterworth(np, N, 0., 0., Rate, ALLPASS, MAXORDER, &BW);

    // determine which notches to use according to high and lowpass
    for (n=1; n<=MAXNOTCH; n++) {
        fc = (double)n * Mains;
        if (fc >= HPFreq && fc <= LPFreq) {
            fh = fc - NOTCHW_2;
            fl = fc + NOTCHW_2;
            Butterworth(&Notch, N, fh, fl, Rate, BANDREJECT, MAXORDER, &BW);
            for (c=0; c<N; c++)
                np->Gain[c] *= Notch.Gain[c];
        }
    }
}
