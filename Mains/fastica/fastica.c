// fastica -- independent components analysis, using algorithms developed at
//	Helsinki University of Technology, by Hyvarinen
//
//	command line flags:
//		-r <run name>	-- MEG run name
//		-n <rank>	-- number of independent components
//		-f <hp, lp>	-- high & lowpass frequencies
//		-g <function>	-- nonlinear function used for ICA
//		-p <precision>	-- number of digits for convergence
//		-b <order>	-- baseline removal order
//		-d		-- delete backup meg4 file (optional)
//		-v		-- verbose
//
// Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <byteswap.h>
#include <string.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>
#include <model.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <SAMfiles.h>
#include <DataFiles.h>
#include <samUtils.h>
#include <version.h>

#define MINOR_REV	0
#define IIRORDER	4
#define MAXITS		100			// maximum number of iterations
#define MAXRETRIES	50			// maximum number of retries before terminating search
#define MAXGUESS	10			// maximum number of successive guesses with declining dot product tolerance
#define DIGITS		4			// default: 4 digits of precision
#define BASELINE	0			// default: remove mean
#define SIG_ORDER	4			// 4th-order IIR filter


int
main(
	int		argc,
	char	**argv

)

{
	FFTSPEC			DataFFT;			// FFT data filter structure
	HeaderInfo		Header;				// header information data
	ChannelInfo		*Channel;			// Channel[M] -- channel information data
	EpochInfo		*Epoch;				// Epoch[E] -- epoch information data
	int32_t			*xint;				// xint[T] -- single channel, one-trial time-series data
	gsl_matrix		*B;					// B - MxM mixing & unmixing matrix
	gsl_matrix		*Binv;				// B^-1
	gsl_matrix		*A;					// A - MxM resultant of BB'
	gsl_matrix		*Tmp;				// temporary matrix
	double			**data;				// data[M][TT] -- filtered data time-series (input)
	double			**x;				// x[M][TT] -- sphered data time-series
	double			**y;				// y[M][TT] -- ICA output time-series
	double			**C;				// C[M][M] -- covariance matrix
	double			**S;				// S[M][M] -- sphering matrix
	double			*u;					// u[M] -- updated weight vector
	double			*v;					// v[M] -- intermediate weight vector
	double			*w;					// w[M] -- weight vector
	double			*hypTan;			// hypTan[TT] - hyperbolic tangent time-series
	double			*In;				// In[T] -- single-channel time-series
	double			*Out;				// Out[T] -- single-channel time-series
	double			*L;					// L[M] -- eigenvalues
	double			scale;				// SQUID conversion bits-->Tesla
	double			sum;				// summation of tanh() transformed data
	double			r;
	double			r2;
	double			rold;
	double			ntt;				// total samples
	double			Gauss;
	double			z;
	double			z2;
	double			ex;
	double			dot;
	double			max;
	double			BWFreq;				// filter effective bandwidth
	double			tolerance = .001;	// convergence tolerance
	static double	a1 = 1.;
	static double	a2 = 1.;
	static double	HighPass = 0.;		// highpass frequency
	static double	LowPass = 0.;		// lowpass frequency
	register double	tmp;				// the famous temporary variable of hacker mythology
	int				*discard;			// discard[M] -- array of IC channel numbers to be discarded
	int				D;					// number of discarded ICs
	int				E;					// number of good epochs
	int				M;					// number of good channels
	int				T;					// samples per epoch per channel
	int				TT;					// total samples
	int				d;					// discarded IC index
	int				e;					// epoch index
	int				f;					// 1st primary channel index
	int				m;					// sensor index
	int				n;					// 'n' stands for 'nuther sensor index
	int				t;					// time index
	int				tt;					// dataset time index (multiple epochs)
	int				i;					// a very popular index
	int				j;					// an oft-times useful index
	int				c;					// a character
	int				its;				// iterations
	int				guess;				// number of bad guesses
	int				retry;				// number of retries
	unsigned char	**Bad;				// Bad[E][T] -- flags for bad samples
	static int		digits = DIGITS;	// default number of digits for convergence
	static int		function = 1;		// default nonlinear function is kurtosis (1)
	static int		K = 0;				// rank
	static int		baseline = BASELINE;	// default trend removal
	static int		rflg = FALSE;		// dataset name flag
	static int		vflg = FALSE;		// verbose mode flag
	int32_t			Ltmp;				// temporary long integer variable
	extern char		*optarg;
	extern int		optind;
	extern int		opterr;
	char			DSPath[256];		// data set name full path
	char			DSName[256];		// MEG dataset name
	char			fpath[256];			// general path name
	char			OldName[256];		// old .meg4 file name
	char			NewName[256];		// new .meg4 file name
	char			StringField[256];	// reserved for searches & response
	char			FileID[9];			// file ID string
	static char		MEG4Hdr[] =			// MEG 4.1 header string
		"MEG41CP";
	FILE			*fin;				// input file pointer
	FILE			*fout;				// output file pointer
	FILE			*fp;				// file pointer for log
	off_t			offset;				// number of bytes offset

	// parse command line parameters
	opterr = 1;		// enable error reporting
	while((c = getopt(argc, argv, "r:f:g:p:b:v")) != EOF) {
		switch(c) {

			case 'r':	// data set directory name (w/o '.ds')
				sprintf(DSName, "%s", optarg);
				rflg = TRUE;
				break;

			case 'f':	// output high & low pass filter frequencies
				if(sscanf(optarg, "%lf%lf", &HighPass, &LowPass) != 2)
					Cleanup("'-f' flag requires high and lowpass frequencies, in quotes");
				if(HighPass >= LowPass)
					Cleanup("highpass frequency must be less than lowpass frequency");
				break;

			case 'g':	// nonlinear function
				function = atoi(optarg);
				if(function < 1 || function > 3)
					Cleanup("nonlinear function must be in range 1 to 3");
				break;

			case 'p':	// precision
				digits = atoi(optarg);
				if(digits < 2 || digits > 10)
					Cleanup("number of significant digits for convergence must be in range 2 to 10");
				break;

			case 'b':	// baseline removal order
				baseline = atoi(optarg);
				if(baseline < 0 || baseline > 4)
					Cleanup("baseline removal order must be in range from 0 to 4");
				break;

			case 'v':	// verbose operation
				vflg = TRUE;
				break;

			default:
				Cleanup("unknown option flag");
				break;
		}
	}
	if(rflg == FALSE) {
		fprintf(stderr, "usage: fastica\t-r <dataset name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "\t\t-f <highpass, lowpass (Hz)>\n");
		fprintf(stderr, "\t\t-g <nonlinear function value>\n");
		fprintf(stderr, "\t\t\t1 -- kurtosis\n");
		fprintf(stderr, "\t\t\t2 -- hyperbolic tangent\n");
		fprintf(stderr, "\t\t\t3 -- gaussian\n");
		fprintf(stderr, "\t\t-p <number of digits for convergence tolerance>\n");
		fprintf(stderr, "\t\t-b <baseline removal order>\n");
		fprintf(stderr, "\t\t\t0 -- remove mean\n");
		fprintf(stderr, "\t\t\t1 -- remove trend\n");
		fprintf(stderr, "\t\t\t2 -- remove quadratic\n");
		fprintf(stderr, "\t\t\t3 -- remove cubic\n");
		fprintf(stderr, "\t\t\t4 -- remove quartic\n");
		fprintf(stderr, "\t\t-v -- verbose mode\n");
		exit(-1);
	}

	// announce & get data with sensor structures
	if(vflg == TRUE) {
		printf("F A S T   I N D E P E N D E N T   C O M P O N E N T S   A N A L Y S I S\t%s%d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		printf("\tusing ");
		switch(function) {
			case 1: printf("kurtosis"); break;
			case 2: printf("hyperbolic tangent"); break;
			case 3: printf("gaussian"); break;
			default: Cleanup("no such nonlinear function"); break;
		}
		printf(" nonlinearity\n");
		tolerance = pow((double)10., (double)(-digits));
		printf("\ttolerance for convergence will be %e (%d digits)\n", (double)tolerance, digits);
		printf("opening data file");
		fflush(stdout);
	}

	// get data with sensor structures
	GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);
	sprintf(DSPath, "%s/%s.ds", Header.DsPath, Header.SetName);

	// extract constants
	E = Header.NumEpochs;
	M = Header.NumPri;
	T = Header.MaxSamples;
	TT = T * E;
	ntt = (double)TT;
	tolerance = pow((double)10., (double)(-digits));
	f = Header.PsIndex[0];		// 1st MEG primary sensor channel index
	if(HighPass == 0. && LowPass == 0.) {
		HighPass = Channel[f].AnalogHpFreq;
		LowPass = Channel[f].AnalogLpFreq;
	}

	// open ica logfile
	sprintf(fpath, "%s-fastica.log", DSName);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open logfile");
	fprintf(fp, "Fastica Log File for dataset: %s\n", DSName);
	fprintf(fp, "Number of Channels: %d\n", M);
	fprintf(fp, "Number of Epochs: %d\n", E);
	fprintf(fp, "Number of Samples/Epoch: %d\n", T);
	fprintf(fp, "Sample Rate: %7.2f Hz\n", Header.SampleRate);
	if(LowPass > 0. || HighPass > 0.) {
		fprintf(fp, "Bandpass: %4.1f to %5.1f Hz\n", HighPass, LowPass);	
	} else {
		fprintf(fp, "(No Bandpass Filters)\n");
	}
	fprintf(fp, "Using ");
	switch(function) {
		case 1: fprintf(fp, "Kurtosis"); break;
		case 2: fprintf(fp, "Hyperbolic Tangent"); break;
		case 3: fprintf(fp, "Gaussian"); break;
		default: Cleanup("no such nonlinear function"); break;
	}
	fprintf(fp, " Nonlinearity\n");
	fprintf(fp, "Tolerance for Convergence: %e (%d digits)\n", (double)tolerance, digits);
	fflush(fp);

	// allocate memory for arrays
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("allocating memory");
		fflush(stdout);
	}
	B = gsl_matrix_alloc(M, M);
	Binv = gsl_matrix_alloc(M, M);
	A = gsl_matrix_alloc(M, M);
	Tmp = gsl_matrix_alloc(M, M);
	if((xint = (int32_t *)malloc(T * sizeof(int32_t))) == NULL)
		allocfailed("xint[]");
	if((data = (double **)malloc(M * sizeof(double *))) == NULL)
		allocfailed("data[]");
	for(i=0; i<M; i++)
		if((data[i] = (double *)malloc(TT * sizeof(double))) == NULL)
			allocfailed("data[][]");
	if((x = (double **)malloc(M * sizeof(double *))) == NULL)
		allocfailed("x[]");
	for(i=0; i<M; i++)
		if((x[i] = (double *)malloc(TT * sizeof(double))) == NULL)
			allocfailed("x[][]");
	if((y = (double **)malloc(M * sizeof(double *))) == NULL)
		allocfailed("y[]");
	for(i=0; i<M; i++)
		if((y[i] = (double *)malloc(TT * sizeof(double))) == NULL)
			allocfailed("y[][]");
	if((C = (double **)malloc(M * sizeof(double *))) == NULL)
		allocfailed("C[]");
	for(i=0; i<M; i++)
		if((C[i] = (double *)malloc(M * sizeof(double))) == NULL)
			allocfailed("C[][]");
	if((S = (double **)malloc(M * sizeof(double *))) == NULL)
		allocfailed("S[]");
	for(i=0; i<M; i++)
		if((S[i] = (double *)malloc(M * sizeof(double))) == NULL)
			allocfailed("S[][]");
	if((u = (double *)malloc(M * sizeof(double))) == NULL)
		allocfailed("u[]");
	if((v = (double *)malloc(M * sizeof(double))) == NULL)
		allocfailed("v[]");
	if((w = (double *)malloc(M * sizeof(double))) == NULL)
		allocfailed("w[]");
	if((L = (double *)malloc(M * sizeof(double))) == NULL)
		allocfailed("L[]");
	if((hypTan = (double *)malloc(TT * sizeof(double))) == NULL)
		allocfailed("hypTan[]");
	if((In = (double *)malloc(T * sizeof(double))) == NULL)
		allocfailed("In[]");
	if((Out = (double *)malloc(T * sizeof(double))) == NULL)
		allocfailed("Out[]");
	if((discard = (int *)malloc(M * sizeof(int))) == NULL)
		allocfailed("discard[]");

	// make filter
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("making %6.3f to %6.3f Hz data filter", HighPass, LowPass);
		fflush(stdout);
	}
	Butterworth(&DataFFT, T, HighPass, LowPass, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

	// for each epoch
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("reading data");
		fflush(stdout);
	}

	// read & filter MEG data
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("reading & filtering MEG data");
		fflush(stdout);
	}
	for(e=0; e<E; e++) {
		for(i=0; i<M; i++) {

			// get one epoch of data
			n = Header.PsIndex[i];
			scale = Channel[n].Scale;
			GetDsData(&Header, e, n, scale, In);

			// apply bandpass filter & remove baseline
			FFTfilter(In, Out, &DataFFT, T);
			debase(Out, T, baseline);

			// move data to y[][]
			for(t=0, tt=e*T; t<T; t++, tt++)
				data[i][tt] = Out[t];
		}
		if(vflg == TRUE) {			// if req'd - epoch heartbeat
			printf(".");
			fflush(stdout);
		}
	}

	// sphere the data & remove baseline
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("sphering data");
		fflush(stdout);
	}
	whiten(data, x, M, &K, TT);			// x is the whitened (sphered) data
	for(i=0; i<M; i++)
		debase(x[i], TT, baseline);

	// initialize mixing matrix
	gsl_matrix_set_identity(B);

	// recursively find each component
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("finding independent components:\n");
		fflush(stdout);
	}
	for(m=0; m<M; m++) {	// assume finding all M source components
		if(vflg == TRUE) {
			printf("component %d:\n", m+1);
			fflush(stdout);
		}
		fprintf(fp, "\nIndependent Component %d:\n", m+1);
		fflush(fp);
		retry = 0;
restart:

		// compute random unit weight vector
		for(i=0, tmp=0.; i<M; i++) {
			w[i] = RandomGauss(0., 1.);
			tmp += w[i] * w[i];
		}
		r = sqrt(tmp);
		for(i=0; i<M; i++)
			w[i] /= r;

		// compute A = I - B B', where B is the current mixing matrix
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., B, B, 0., Tmp);
		gsl_matrix_set_identity(A);
		gsl_matrix_sub(A, Tmp);

		// remove projections of previous independent components
		for(i=0, r2=0.; i<M; i++) {
			for(j=0, tmp=0.; j<M; j++)
				tmp += w[j] * gsl_matrix_get(A, j, i);
			u[i] = tmp;
			r2 += tmp * tmp;
		}

		// make w a unit vector
		r = sqrt(r2);
		for(i=0; i<M; i++)
			w[i] = u[i] / r;

		// log seed weights
		fprintf(fp, "\tTry %d, Seed Weight Vector:\n", retry);
		for(i=0; i<M; i++) {
			fprintf(fp, "W[%d]=%e\n", i, w[i]);
			fflush(fp);
		}

		// iterate weight vector
		rold = 0.;
		its = guess = 0;
		do {

			// select nonlinearity
			switch(function) {

				case 1:		// compute kurtosis nonlinearity
					for(i=0; i<M; i++) {
						for(t=0, u[i]=0.; t<TT; t++) {
							for(j=0, tmp=0.; j<M; j++)
								tmp += w[j] * x[j][t];
							u[i] += tmp * tmp * tmp * x[i][t];
						}
						u[i] = u[i] / ntt - 3. * w[i];
					}
					break;

				case 2:		// compute hyperbolic tangent nonlinearity
					for(t=0, sum=0.; t<TT; t++) {
						for(i=0, tmp=0.; i<M; i++)
							tmp += x[i][t] * w[i];
						hypTan[t] = tanh(a1 * tmp);
						sum += 1. - hypTan[t] * hypTan[t];
					}

					// compute w = (x * hypTan - a1 * sum(1 - hypTan^2)' * w) / ntt
					for(i=0; i<M; i++) {
						for(t=0, tmp=0.; t<TT; t++)
							tmp += x[i][t] * hypTan[t];
						dot = tmp;
						u[i] = (dot - a1 * sum * w[i]) / ntt;
					}
					break;

				case 3:		// compute Gaussian nonlinearity
					for(i=0; i<M; i++)
						u[i] = 0.;
					for(t=0, sum=0.; t<TT; t++) {
						for(i=0, tmp=0.; i<M; i++)
							tmp += x[i][t] * w[i];
						z = tmp;
						z2 = tmp * tmp;
						ex = exp(-a2 * z2 / 2.);
						Gauss = z * ex;				// Gauss[t]
						sum += (1. - a2 * z2) * ex;	// dGauss[t]
						for(i=0; i<M; i++)
							u[i] += x[i][t] * Gauss;
					}
					for(i=0; i<M; i++) {
						u[i] -= sum * w[i];
						u[i] /= ntt;
					}
					break;

				default:
					Cleanup("no such nonlinear function");
					break;
			}

			// remove projections of previous independent components
			for(i=0, r2=0.; i<M; i++) {
				for(j=0, tmp=0.; j<M; j++)
					tmp += u[j] * gsl_matrix_get(A, j, i);
				v[i] = tmp;
				r2 += tmp * tmp;
			}

			// normalize u, compute dot product with w, & update w
			r = sqrt(r2);
			for(i=0, tmp=0.; i<M; i++) {
				u[i] = v[i] / r;
				tmp += u[i] * w[i];
				w[i] = u[i];
			}

			r = sqrt(fabs(tmp));
			its++;
			if(vflg == TRUE) {
				printf("\tpass %d: norm=%13.10f\n", its, r);
				fflush(stdout);
			}

			// check for bad guess
			if(r <= rold)
				guess++;
			else
				guess = 0;
			if(guess > MAXGUESS) {
				if(retry > MAXRETRIES) {
					K = m - 1;
					printf("\t\t(terminating search at %d components)\n", K);
					fflush(stdout);
					fprintf(fp, "search terminated at %d components\n", K);
					fflush(fp);
					break;
				}
				if(vflg == TRUE) {
					printf("\t\t(restart component %d after %d retries)\n", m+1, retry);
					fflush(stdout);
				}
				retry++;
				goto restart;
			} else
				rold = r;

		} while((1. - r) > tolerance && its < MAXITS);

		// update mixing matrix B
		for(i=0; i<M; i++)
			gsl_matrix_set(B, i, m, w[i]);

		// log final weights
		fprintf(fp, "\tComponent %d ICA Weight Vector:\n", m+1);
		for(i=0; i<M; i++) {
			fprintf(fp, "W[%d]=%e\n", i, w[i]);
			fflush(fp);
		}
	}	// end M components loop

	// project sphered data through rotation matrix
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("projecting sphered data through unmixing matrix");
		fflush(stdout);
	}
	for(t=0; t<TT; t++)
		for(i=0; i<M; i++) {
			for(j=0, tmp=0.; j<M; j++)
				tmp += gsl_matrix_get(B, j, i) * x[j][t];
			y[i][t] = tmp;
		}

	// find maximum/minimum for each IC & scale to fit INT_MAX/2
	for(i=0; i<M; i++) {
		scale = 0.5 * (double)INT_MAX * Channel[m].Scale;	// 1/2 full scale value at current gain
		for(t=0, max=0.; t<TT; t++) {
			tmp = fabs(y[i][t]);
			if(tmp > max)
				max = tmp;
		}
		scale /= max;				// output scale factor
		for(t=0; t<TT; t++)
			y[i][t] *= scale;
	}

	if(vflg == TRUE) {
		printf(" - done\n");
		printf("opening data files");
		fflush(stdout);
	}
	sprintf(fpath, "%s/%s.meg4", DSPath, DSName);
	if((fin = fopen(fpath, "r")) == NULL)
		Cleanup("can't open '.meg4' file for read");
	if(fread((void *)FileID, 8, 1, fin) != 1)
		Cleanup("can't read header from '.meg4' file");
	if(strncmp(FileID, "MEG4", 4))
		Cleanup("unknown header ID in '.meg4' file");
	sprintf(fpath, "%s/%s.meg4.ICA", DSPath, DSName);
	if((fout = fopen(fpath, "w")) == NULL)
		Cleanup("can't open '.meg4.ICA' file for write");
	if(fwrite((void *)MEG4Hdr, 8, 1, fout) != 1)
		Cleanup("can't write header to '.meg4.ICA' file");

	// write each epoch of data
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing ICs");
		fflush(stdout);
	}
	for(e=0; e<E; e++) {
		for(m=i=0; m<Header.NumChannels; m++) {
			if(m == Header.PsIndex[i]) {
				for(t=0, tt=e*T; t<T; t++, tt++) {
					Ltmp = (int32_t)(y[i][tt] / Channel[m].Scale);
					xint[t] = bswap_32(Ltmp);
				}
				if(fwrite((void *)xint, sizeof(int32_t), T, fout) != T)
					Cleanup("can't write data to '.meg4' file");
				i++;
			} else {		// channel isn't in primary channels list -- just copy it...
				offset = (off_t)((T * sizeof(int32_t)) * (m + e * Header.NumChannels) + 8);
				if(fseek(fin, offset, SEEK_SET) != 0)
					Cleanup("'fseek()' to channel failed");
				if(fread((void *)xint, sizeof(int32_t), T, fin) != T)
					Cleanup("can't read data from 'meg4' file");
				if(fwrite((void *)xint, sizeof(int32_t), T, fout) != T)
					Cleanup("can't write data to '.meg4' file");
			}
		}

		// hearbeat
		if(vflg == TRUE) {
			printf(".");
			fflush(stdout);
		}
	}
	fclose(fin);
	fclose(fout);
	fclose(fp);

	// rename original '.meg4' file
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("moving & renaming '.meg4' file");
		fflush(stdout);
	}
	sprintf(OldName, "%s/%s.meg4", DSPath, DSName);
	sprintf(NewName, "%s/%s.meg4.bak", DSPath, DSName);
	if(rename(OldName, NewName) == -1)
		Cleanup("backup of original '.meg4' file failed");
	sprintf(OldName, "%s/%s.meg4.ICA", DSPath, DSName);
	sprintf(NewName, "%s/%s.meg4", DSPath, DSName);
	if(rename(OldName, NewName) == -1)
		Cleanup("rename of '.meg4.ICA' file failed");

	// invert unmixing matrix
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("inverting unmixing matrix");
		fflush(stdout);
	}
	pinv(B, Binv);

	// give instructions to remove ICs
	printf(" - done\n");
	printf("(1) open %s.ds using 'DataEditor'\n", DSName);
	printf("(2) inspect transformed channels and mark channels/components\n");
	printf("\tto be removed as bad, using 'Set Good/Bad'\n");
	printf("(3) when done, select 'File' and then 'Save Dataset (Edit Info)'\n");
	printf("(4) hit 'y' when ready to proceed\n");
	fflush(stdout);

	// wait for approval to remove ICs
	printf("\nProceed? (y/n)\n");
	scanf("%s", StringField);

	// open 'BadChannels' & read names of ICs to discard
	for(m=0; m<M; m++)			// initialize list
		discard[m] = FALSE;
	sprintf(fpath, "%s/BadChannels", DSPath);
	d = 0;
	if((fp = fopen(fpath, "r")) != NULL) {
		while(fgets(StringField, 128, fp) != NULL) {
			j = strlen(StringField) - 1;
			for(m=0; m<M; m++) {			// compare name with channel name
				if(!strncmp(StringField, Channel[m].ChannelName, j)) {
					discard[m] = TRUE;		// set flag to discard this IC
					d++;					// count number of discards
					break;
				}
			}
		}
	} else {
		Cleanup("'BadChannels' file not found");
	}
	fclose(fp);
	D = d;								// save number of discarded ICs

	// make sure there are ICs!
	if(D == 0)
		Cleanup("There were no ICA components identified in 'BadChannels'");

	// echo discard list
	printf("\nremoving IC projections from %d channels:\n", D);
	for(m=0; m<M; m++)
		if(discard[m] == TRUE)
			printf("[%d] %s\n", m+1, Channel[m].ChannelName);
	fflush(stdout);

	// for each channel, sum the artifact-free ICs through the mixing matrix
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("reconstructing data without artifact ICs");
		fflush(stdout);
	}
	for(t=0; t<TT; t++)
		for(i=0; i<M; i++) {
			for(j=0, tmp=0.; j<M; j++)
				if(discard[j] == FALSE)		// sum good ICs, only
					tmp += gsl_matrix_get(Binv, i, j) * data[i][t];
			x[i][t] = tmp;
		}

	// write data back to .meg4 file
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing cleaned '.meg4' file");
		fflush(stdout);
	}
	sprintf(fpath, "%s/%s.meg4.bak", DSPath, DSName);
	if((fin = fopen(fpath, "r")) == NULL)
		Cleanup("can't open '.meg4.bak' file for read");
	if(fread((void *)FileID, 8, 1, fin) != 1)
		Cleanup("can't read header from '.meg4' file");
	if(strncmp(FileID, "MEG4", 4))
		Cleanup("unknown header ID in '.meg4' file");
	sprintf(fpath, "%s/%s.meg4", DSPath, DSName);
	if((fout = fopen(fpath, "w")) == NULL)
		Cleanup("can't open '.meg4' file for write");
	if(fwrite((void *)MEG4Hdr, 8, 1, fout) != 1)
		Cleanup("can't write header to '.meg4' file");
	for(e=0; e<E; e++) {
		for(m=i=0; m<Header.NumChannels; m++) {
			if(m == Header.PsIndex[i]) {
				for(t=0, tt=e*T; t<T; t++, tt++) {
					Ltmp = (int32_t)(x[i][tt] / Channel[m].Scale);
					xint[t] = bswap_32(Ltmp);
				}
				if(fwrite((void *)xint, sizeof(int32_t), T, fout) != T)
					Cleanup("can't write data to '.meg4' file");
				i++;
			} else {		// channel isn't in primary channels list -- just copy it...
				offset = (off_t)((T * sizeof(int32_t)) * (m + e * Header.NumChannels) + 8);
				if(fseek(fin, offset, SEEK_SET) != 0)
					Cleanup("'fseek()' to channel failed");
				if(fread((void *)xint, sizeof(int32_t), T, fin) != T)
					Cleanup("can't read data from 'meg4' file");
				if(fwrite((void *)xint, sizeof(int32_t), T, fout) != T)
					Cleanup("can't write data to '.meg4' file");
			}
		}

		// hearbeat
		if(vflg == TRUE) {
			printf(".");
			fflush(stdout);
		}
	}
	fclose(fout);

	if(vflg == TRUE) {
		printf(" - done\n");
		printf("fastica: done\n");
		fflush(stdout);
	}
	exit(0);
}
