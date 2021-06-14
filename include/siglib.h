#ifndef H_SIGLIB
#define H_SIGLIB

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <filters.h>
#include <complex.h>
#include <geoms.h>
#include <fftw3.h>
#include <time.h>
#include "samutil.h"

#define BASE_MEAN   0
#define BASE_TREND  1
#define BASE_QUAD   2
#define BASE_CUBE   3

typedef struct {
    int len;
    fftw_plan p1;
    fftw_plan p2;
    fftw_complex *x;
    fftw_complex *X;
} FFTWPLAN;

void    balanc(double **, int);
void    bdiir(double *, double *, int, IIRSPEC *);
double  betacf(double, double, double);
double  betai(double, double, double);
void    Butterworth(FFTSPEC *Filter, int N, double HighPassFreq, double LowPassFreq, double SampleRate, int FilterType, int Order, double *bw);
int     chirp(double *, int, double, double, double, double);
void    comb(double *, int, double, double);
int     corr(double **, int, int, int, int, double **);
int     cov(double **, int, int, int, int, double **);
int     cspline(double *, double *, int, int);
void    debase(double *, int, int);
void    demean(double *, int);
void    dequad(double *, int);
void    denoise(double **, int, int *, int, int *, int);
void    derivs(double x, double *y, double dydx);
void    detect(double *x, int N, IIRSPEC *Filter);
void    detrend(double *, int);
int     dft(double *, double *, int, double, double, double, double);
void    dswap(double *, double *);
void    eigen(double **A, double *B, int N);
void    eig22(double e[2][2], double l[2]);
int     eig33(double (*)[3], double *);
void    equalizer(double *data,double *result, FFTSPEC *filter, double fhi, double flo, double srate, int len);
// double  erf(double);
// double  erfc(double);
double  evlmem(double, double *, int, double);
void    evolve(double *x, double *h, int TS, int TE, double LPFreq, double tau, double rate);
int     extendts(double *, double *, int, int);
void    FFTfilter(double *data, double *result, FFTSPEC *Filter, int len);
void    FFTrand(double *In, double *Out, int len);
int     filter(double *, double *, int, IIRSPEC *);
void    fir(double *, double *, int);
int     fixrts(double *, int);
int     fmf(double *, double *, int, int);
void    ForwardDft(double _Complex *, double _Complex *, int);
double  freqinfo(double *x, double *y, int TS, int TE, double HPFreq, double LPFreq, double rate);
double  gammln(double);
double  gammp(double, double);
double  gammq(double, double);
int     gcf(double *, double, double, double *);
int     gser(double *, double, double, double *);
void    Hanning(double *, int);
FFTWPLAN *new_fftwplan(int len, int do_wis);
void    free_fftwplan(FFTWPLAN *);
void    hilbert(FFTWPLAN *, double *data, double *result);
void    iir(double *, double *, int, IIRSPEC *);
int     inv22(double A[2][2]);
int     inv33(double A[3][3]);
void    invert(gsl_matrix *A, gsl_matrix *Ainv);
double  kendall(double *x, double *y, unsigned long n);
double  kurtosis(double *, int, int);
int     laguerre(double _Complex *, int, double _Complex *, double, int);
void    mab(double *, int, double *, double *);
void    sam(double **, double *, int, int, double *);
double  legendre(int, double);
double	Levy(double x, double mu, double c);
int     lpgap(double *, double *, int, int, int, double *, int, int);
int     marker(double **, int, double, int *, int, double, double, double *);
double  median(double *, int);
void    memcof(double *, int, int, double *, double *);
void    mkfilter(IIRSPEC *, int, double, double, double, double *);
void    mkiir(IIRSPEC *);
void    mknotch(FFTSPEC *, double HPFreq, double LPFreq, double Rate, int N, double Mains);
double  modal(double *, int, int);
double  mode(double *, int, int);
void    Morlet(double *, int);
int     mr1df(double, double*, double*, int);
void    mutinfo(double *x, double *y, double *Info, int TS, int TE, double LPx, double LPy, double tau, double rate);
double  mutinfostat(double *x, double *y, double LPx, double LPy, int TS, int TE, double rate);
int     notch(double *, int, double, double);
int     outlier(double *, double *, int, double *, int);
double  pax(double *, double *, int, int);
void    pca22(double **x, int T);
double  peak(double *, int, double, double, double, double);
double  PeakSlope(double *, int);
void    phevolve(double *x, double *h, int TS, int TE, double LPFreq, double tau, double rate);
void    pinv(gsl_matrix *A, gsl_matrix *Ainv);
void    pinvN(gsl_matrix *A, gsl_matrix *Ainv, int N);
double  power(double *, int, int);
double  PpRms(double *, int);
int     predict(double *, double *, int, int, double *, int, int);
void    predictability(double *x, double *lambda, int TS, int TE, int W, int F, int P);
int     ProbF(double, int, double, int, double *, double *);
int     ProbT(double, double, int, double, double, int, double *, double *);
int     ProbZ(double, double, int, double *, double *);
void    pwr(double *x, double *p, int TS, int TE, double HPFreq, double rate);
void    QtoM(QUATERNION quaternion, double Rotate[3][3]);
void    rk4(double *y, double *dydx, int n, double x, double h, double *yout, void (*derivs)(double, double *, double *));
void    rkdumb(double *vstart, int nvar, double x1, double x2, int nstep, void (*derivs)(double, double *, double *));
int     realfft(double, double *, double *, int);
double  response(IIRSPEC *, double, double);
double  rms(double *, int);
int     robavg(double ***, double **, double **, double *, double, int, int, int, int *, int *);
int     robvar(double ***, double **, double **, double *, double, int, int, int);
void    set_wisfile(void);
double  slope(short int *, int, int, int, int);
void    SolveEigen(double **A, double **B, double *Vec, double *Val, int N,  int MinMax);
void    spectrum(double *x, double *s, double srate, int N);
//void    st(int n, int lo, int hi, double *data, double *result);
//int     st_freq(double f, int len, double srate);
void    SVDcmp(SVD *);
void    SVDallocate(SVD *);
void    SVDfree(SVD *);
void    SVDLcmp(SVDL *);
void    SVDLallocate(SVDL *);
void    SVDLfree(SVDL *);
double  vfcorr(float *, float *, int);
double  vfdist(float *, float *, int);
int     vlag(double **, double **, double *, int, int, int, int, int, double, double *);
double  vpax(double **, double **, int, int, int);
double  vxcorr(double **, double **, int, int, int);
double  vxreg(double **, double **, int, int, int);
void    whiten(double **x, double **y, int N, int *K, int T);
int     widrow(double *, int, double, double, double);
int     widrowBi(double *, int, double, double, double);
double  xCorr(double **, double **, int, int);
void    xmtx(double **, int, int, int *, int, int *, int, double **, double **, double *, gsl_matrix *);
void    YtoR(double theta, double rotation[3][3]);
int     zpf(double *, double *, int, IIRSPEC *);
int     zroots(double _Complex *, int, double _Complex *, int);

// for computing Spearman rank correlation -- absent in gsl-1.15 but present in gsl-2.5
double  spearman(double *data1, double *data2, int n);
void compute_rank(gsl_vector *v);

#if !defined(log2)
double  log2(double x); /* sometimes this is not in math.h */
#endif

#endif  // H_SIGLIB

