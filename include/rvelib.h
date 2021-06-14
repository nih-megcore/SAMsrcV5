#ifndef H_RVELIB
#define H_RVELIB

#include <geoms.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <lfulib.h>
#include <rvelib.h>
#include <TimeSegs.h>
#include "samutil.h"

#define FWD	0
#define REV	1

void	Aug2Sym(double *x, double *y, short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, int T, int *TE, int E, int *N, int advance, double LPFreq, double rate);
void	Data2CondRVE(double *x, double *h, int T, int E, double LPFreq, double tau, double advance, double rate);
void	Data2RVE(double *x, double *h, int T, int E, double LPFreq, double tau, double rate);
void	Data2Sym(double *x, double *y, short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, int T, int *TE, int E, int advance, double LPFreq, double rate);
void	Data2SymbolArray(double *x, short *Symbol, int T, int *TE, int E, double LPFreq, double rate);
void	Data2SymbolsArray(double *x, short *SymPhs, short *SymAdv, int T, int *TE, int E, int advance, double LPFreq, double rate);
void	logistic(double *x, double *y, double ratio, int T);
void	Lorenz(double *Out, double dt, int T);
void	MackeyGlassX(double *x, int T);
void	MackeyGlassXY(double *x, double *y, double coupling, int T);
double	MeshInfo(double **x, double LPFreq, double **C, HeaderInfo *Header, ChannelInfo *NewChannel, HULL *Mesh);
double	MeshTE(double **x, double LPFreq, double **C, HeaderInfo *Header, ChannelInfo *NewChannel, HULL *Mesh);
double	MGdx(double x_now, double x_tau, double time);
double	MGdy(double y_now, double y_tau, double x_now, double x_tau, double coupling, double time);
void	RandomSurrogate(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, int E, int T);
void	RandomSymbol(short *, int, int);
void	RandomSymbol4(short *Sym4, short *Sym1, int);
void	RandomSymbol5(short *Sym5, short *Sym1, int);
void	randsample(double *, double *, int);
double	rk4x(double x_now, double x_tau, double time, double step);
double	rk4y(double x_now, double x_tau, double y_now, double y_tau, double coupling, double time, double step);
double	rvestatic(double *x, int TS, int TE, int E, double LPFreq, double rate);
void	Sym2AvgSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, int E, double *Hx2y, double *Hy2x, TRIAL *Trial, int M, int W);
void	Sym2erSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, TRIAL *Trial, int M, float *Hfwd, float *Hrev, int E);
void	Sym2H(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, double rate, double *H, int Direction, int E, int T);
void	Sym2InstSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, int tau, float *H, int Direction, int E, int T);
void	Sym2InstSTE2(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, int tau, float *H, int Direction, int E, int T);
double	Sym2MeanFwdSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, double tau, double rate, int E, int T);
double	Sym2MeanRevSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, double tau, double rate, int E, int T);
double	Sym2MI(short *SymPhsX, short *SymPhsY, int TE, int advance, int E);
void	Sym2RVE(short *Sym, double *h, int T, int E);
void	Sym2TrialRVE(short	*SymX, TRIAL *Trial, int N, double *Hx, double *Hx2, int T, int E);
void	Sym2STE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, float *Hx2y, float *Hy2x, int E, int T);
void	Sym2STEi(short *SymPhsX, short *SymPhsY, short *SymAdvY, double tau, double rate, float *Hxy, int T, int E);
void	Sym2STEint(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, double tau, double rate, float *H, int Direction, int E, int T);
void	Sym2tdSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, double tau, double rate, float *H, int Direction, int E, int T);
void	Sym2tdSTE2(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, double tau, double rate, float *Hfwd, float *Hrev, int E, int T);
void	Sym2tdSTE3(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, float alpha, float *Hfwd, float *Hrev, int T);
void	Sym2Fwd_tdSTE(short *SymPhsX, short *SymPhsY, short *SymAdvY, float alpha, float *H, int T);
void	Sym2TrialSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, TRIAL *Trial, int M, float *Hx2y, float *Hy2x, int E);
void	Sym2TrialInstSTE(short *SymPhsX, short *SymAdvX, short *SymPhsY, short *SymAdvY, TRIAL *Trial, int M, float *Ax2y, float *Ay2x, int E, int T);
short	Vec2Advance(double *x, int E);
short	Vec2Augmented(double *x, int Tau);
short	Vec2Symbol(double *x, int E);
void	W4toMIt(short *SymbolX, short *SymbolY, double *info, double tau, double rate, int TL, int advance);
void	W5toMIt(short *SymbolX, short *SymbolY, double *info, double tau, double rate, int TL, int advance);
void	W6toMIt(short *SymbolX, short *SymbolY, double *info, double tau, double rate, int TL, int advance);
void	W4toRVE(double *x, double *h, int TS, int TE, double LPFreq, double tau, double rate);
void	W5toRVE(double *x, double *h, int TS, int TE, double LPFreq, double tau, double rate);
void	W6toRVE(double *x, double *h, int TS, int TE, double LPFreq, double tau, double rate);
void	W4toRVEi(short *Sym4, double tau, double rate, float *Hx, int T);
void	W5toRVEi(short *Sym5, double tau, double rate, float *Hx, int T);
void	W6toRVEi(short *Sym6, double tau, double rate, float *Hx, int T);
void	W4toSTE(short *SymX, short *SymY, double *tfwd, double *trev, int T, int delta);
double	W4toSTE_f(short *SymX, short *SymY, int T, int delta);
void	W41toSTE(short *Sym4X, short *Sym1X, short *Sym4Y, short *Sym1Y, float *tfwd, float *trev, int T);
void	W51toSTE(short *Sym5X, short *Sym1X, short *Sym5Y, short *Sym1Y, float *tfwd, float *trev, int T);
void	W41toSTEi(short *Sym4X, short *Sym4Y, short *Sym1Y, double tau, double rate, float *Hxy, int T);
void	W51toSTEi(short *Sym5X, short *Sym5Y, short *Sym1Y, double tau, double rate, float *Hxy, int T);
void	W41toSTEint(short *Sym4X, short *Sym1X, short *Sym4Y, short *Sym1Y, double tau, double rate, float *Hfwd, float *Hrev, int T);
void	W51toSTEint(short *Sym5X, short *Sym1X, short *Sym4Y, short *Sym1Y, double tau, double rate, float *Hfwd, float *Hrev, int T);
float	W41toSTEone(short *Sym4X, short *Sym4Y, short *Sym1Y, int T);
float	W51toSTEone(short *Sym5X, short *Sym5Y, short *Sym1Y, int T);
void	W41toSTEt(short *Sym4X, short *Sym4Y, short *Sym1Y, double tau, double rate, float *Fwd, int T);
void	W51toSTEt(short *Sym5X, short *Sym5Y, short *Sym1Y, double tau, double rate, float *Fwd, int T);

#if !defined(log2)
double	log2(double x);	/* sometimes this is not in math.h */
#endif

#endif	// H_RVELIB

