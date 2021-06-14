#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "siglib.h"
#include "samutil.h"

char *Wisfile = NULL;
char *Wistemplate = "%s/.fftwis";
#define WISLEN 8

void set_wisfile(void)
{
    char *home;

    if (Wisfile) return;
    home = getenv("HOME");
    Wisfile = new_string(strlen(home) + WISLEN);
    sprintf(Wisfile, Wistemplate, home);
}

/* Create FFTW plans and work arrays. Call with do_wis TRUE to update the wisdom. */

FFTWPLAN *new_fftwplan(int len, int do_wis)
{
    FFTWPLAN *plan;
    FILE *wisdom;

    plan = new(FFTWPLAN);
    plan->len = len;
    plan->x = fftw_malloc(sizeof(fftw_complex) * len);
    plan->X = fftw_malloc(sizeof(fftw_complex) * len);

    /* Get any accumulated wisdom. */

    if (do_wis) {
        set_wisfile();
        wisdom = fopen(Wisfile, "r");
        if (wisdom) {
            fftw_import_wisdom_from_file(wisdom);
            fclose(wisdom);
        }
    }

    /* Set up the fftw plans. */

    plan->p1 = fftw_plan_dft_1d(len, plan->x, plan->X, FFTW_FORWARD, FFTW_MEASURE);
    plan->p2 = fftw_plan_dft_1d(len, plan->X, plan->x, FFTW_BACKWARD, FFTW_MEASURE);

    /* Save the wisdom. */

    if (do_wis) {
        wisdom = fopen(Wisfile, "w");
        if (wisdom) {
            fftw_export_wisdom_to_file(wisdom);
            fclose(wisdom);
        }
    }

    return plan;
}

void free_fftwplan(FFTWPLAN *plan)
{
    fftw_destroy_plan(plan->p1);
    fftw_destroy_plan(plan->p2);
    fftw_free(plan->x);
    fftw_free(plan->X);
    free(plan);
}

/* The Hilbert transform. Return the imaginary part only in the result. */

void hilbert(FFTWPLAN *plan, double *data, double *result)
{
    int i, len, l2;
    double *p;
    fftw_complex *h, *H;

    /* Convert the input to complex. */

    len = plan->len;
    h = plan->x;
    H = plan->X;
    memset(h, 0, sizeof(fftw_complex) * len);
    for (i = 0; i < len; i++) {
        h[i][0] = data[i];
    }

    /* FFT. */

    fftw_execute(plan->p1); /* h -> H */

    /* Hilbert transform. The upper half-circle gets multiplied by
    two, and the lower half-circle gets set to zero.  The real axis
    is left alone. */

    l2 = (len + 1) / 2;
    for (i = 1; i < l2; i++) {
        H[i][0] *= 2.;
        H[i][1] *= 2.;
    }
    l2 = len / 2 + 1;
    for (i = l2; i < len; i++) {
        H[i][0] = 0.;
        H[i][1] = 0.;
    }

    /* Inverse FFT. */

    fftw_execute(plan->p2); /* H -> h */

    /* Fill in the result. */

    p = result;
    for (i = 0; i < len; i++) {
        *p++ = h[i][1] / len;
    }
}
