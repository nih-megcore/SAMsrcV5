// SAMersc -- synthetic aperture magnetometry imaging for source event related
//		signals. SAMersc & computes the noise-normalized source waveforms that
//		can be transformed to Hilbert envelope, rank vector entropy, or conditonal
//		rank vector entropy. The resulting time series is averaged over all
//		trials and the baseline subtracted. If two markers are designated, the
//		difference is computed. The result is an event-related averaged time-series
//		for each voxel. The average is mapped to a 3d+time image (.nii) file
//
//		flags:
//		-r <dataset name>			-- MEG run name
//		-d <data file name>			-- MEG data file name (4D, only)
//		-m <parameter file name>	-- name of trigger marker to average
//		-z <snr threshold>			-- threshold for SNR of virtual sensor
//		-v							-- verbose mode
//
//		Author: Stephen E. Robinson
//				MEG Core Facility
//				NIMH
//


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <string.h>
#include <fcntl.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <rvelib.h>
#include <Params.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <samUtils.h>
#include <filters.h>
#include <nifti1.h>
#include <version.h>

#define MINOR_REV		5
#define TRUE			1
#define FALSE			0
#define SIG_ORDER		4


typedef struct {							// structure for listing signal averaging trials
	int				Epoch;					// epoch number
	int				TS;						// average starting sample
	int				Marker;					// either Marker1 (or, optionally, Marker2) -- value is 1 or 2
} AVG_TRIAL;


int	main(
	int		argc,
	char	**argv
)

{
	IIRSPEC			DataIIR;				// IIR data filter structure
	IIRSPEC			SmoothIIR;				// IIR smooting filter structure
	FFTSPEC			DataFFT;				// FFT data filter structure
	FFTSPEC			SmoothFFT;				// FFT smoothing filter structure
	FFTSPEC			NotchFFT;				// FFT notch filter structure
	PARMINFO		Params;					// analysis parameters
	AVG_TRIAL		*Trial;					// Trial[N]
	SAM_MARKS		*Marker;				// Marker[N] -- array of time-ordered markers
	HeaderInfo		Header;					// MEG data header
	nifti_1_header	NiiHdr;					// NIFTI header
	unsigned char	*ExtHdr;				// NiiHdr[N] -- extended NIFTI header
	ChannelInfo		*Channel;				// MEG channel info
	EpochInfo		*Epoch;					// MEG epoch info
	gsl_matrix		*Wgt;					// Wgt - VxM] SAM weights
	gsl_vector		*W;						// W - Mx1 SAM weight vector
	double			***x;					// x[E][T][M] -- multisensor data
	double			*avg0;					// avg0[T] -- final voxel average
	double			*avg1;					// avg1[T] -- voxel average for marker1
	double			*avg2;					// avg2[T] -- voxel average for marker2
	double          *In;                    // In[TT] -- input time-series
	double			*Out;					// Out[TT] -- output time-series
	double			baseline;				// baseline for SAM average
	double			rms_noise;				// rms noise of primary sensors projected through SAM weights
	double          BWFreq;                 // effective bandwidth
	double			BWratio;				// ratio of image bandwidth to covariance bandwidth
	double          sigma2;                 // noise power in bandwidth
	double          ta;                     // number of samples in average (float)
	double          tb;                     // number of samples in baseline (float)
	double			tp;						// number of sanokes in polarity (float)
	double			time;					// average or baseline window duration (seconds)
	double			sign;					// polarity of response window
	double          percent;                // percent done
	double 			tmp;                    // what else?
	float			**Image;				// Image[D][V] -- image voxels for each latency
	int             M;                      // number of MEG or reference channels used
	int				N;						// number of characters in extended NIFTI header
	int				E;						// number of epochs
	int             T;                      // number of samples
	int				TT;						// E*T
	int             TA;                     // number of samples in average
	int             TB;                     // number of samples in baseline
	int				TP;						// number of samples in polarity
	int             AS;                     // average starting sample
	int             AE;                     // average end sample
	int             BS;                     // baseline starting sample
	int             BE;                     // baseline ending sample
	int				PS;						// polarity starting time
	int				PE;						// polarity ending time
	int             V;                      // number of voxels
	int             c;                      // input character & virtual channel
	int             pct;                    // integer percentage
	int             m;                      // channel index
	int             n;                      // trigger marker index
	int				d;						// image latency index
	int             e;                      // epoch index
	int				f;						// 1st primary channel index
	int             t;                      // sample index
	int             tt;                     // alternate sample index
	int				ts;						// starting sample index for response
	int             k;                      // k-index
	int				na1;					// number of trials for marker1
	int				na2;					// number of trials for marker2
	int             sample;					// sample index
	int             NumMarkers;				// number of markers
	int             NumTriggers;			// number of triggers matching value specified in command line
	int				TotTrials;				// total number of trials for marker 1 & optionally marker 2
	int             NumTrials1;				// number of trials in average 1
	int				NumTrials2;				// number of trials in average 2
	int				NumImg;					// number of images
	int             LongIndex;              // position of arguments
	static int		aflg = FALSE;			// take absolute value of voxels
	static int		bflg = FALSE;			// baseline removal flag
	static int      eflg = FALSE;           // command-line error flag
	static int      mflg = FALSE;           // parameter file name flag
	static int		pflg = FALSE;			// polarity correction flag
	static int      rflg = FALSE;           // dataset name flag
	static int      vflg = FALSE;           // verbose mode flag
	extern char     *optarg;
	extern int      opterr;
	char            fpath[256];             // general path name
	char            DSName[256];            // MEG dataset name
	char            DSpath[256];            // MEG dataset path
	char            SAMpath[256];           // SAM subdirectory path
	char			ImgPath[256];			// image path
	char			Prefix[64];				// output file prefix
	char			Suffix[10];				// output file suffix
	char			ParmName[256];			// parameter file name
	char			WgtDir[256];			// weight file directory pathname
	char			FilterName[8];			// filter name -- IIR | FFT
	char			extension[4];			// extension
#if BTI
	static int      dflg = FALSE;           // pdf file name flag
	char            PDFName[256];           // name of data file
	char            ShortOpts[] = "r:d:m:av";
	static struct option    LongOpts[] = {  // command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"pdf_name", required_argument, NULL, 'd'},
		{"analysis_name", required_argument, NULL, 'm'},
		{"absolute_value", no_argument, NULL, 'a'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};
#else
	double          scale;                  // bits-->Tesla conversion
	unsigned char   **Bad;                  // Bad[E][T] -- flags for bad samples
	char            ShortOpts[] = "r:m:av";
	static struct option    LongOpts[] = {  // command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"analysis_name", required_argument, NULL, 'm'},
		{"absolute_value", no_argument, NULL, 'a'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};
#endif
	FILE            *fp;					// output file pointer
	FILE			*np;					// noise file pointer

	// parse command line parameters
	opterr = 1;             // enable error reporting
	while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
		switch(c) {
			case 'r':       // run name
				sprintf(DSName, "%s", optarg);
				rflg = TRUE;
				break;
#if BTI
			case 'd':       // data file name
				sprintf(PDFName, "%s", optarg);
				dflg = TRUE;
				break;
#endif
			case 'm':       // analysis parameter file name
				sprintf(ParmName, "%s", optarg);
				mflg = TRUE;
				break;
			case 'a':		// absolute value of voxels
				aflg = TRUE;
				break;
			case 'v':       // verbose mode
				vflg = TRUE;
				break;
			default:
				eflg = TRUE;
				break;
		}
	}
#if BTI
	if(eflg == TRUE || rflg == FALSE || dflg == FALSE || mflg == FALSE) {
		fprintf(stderr, "SAMersc\t-r <run name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "\t-d <data file name>\n");
#else
		if(eflg == TRUE || rflg == FALSE || mflg == FALSE) {
			fprintf(stderr, "SAMersc\t-r <run name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
#endif
		fprintf(stderr, "\t-m <parameter file name>\n");
		fprintf(stderr, "\t-a -- absolute voxel values\n");
		fprintf(stderr, "\t-v -- verbose mode\n");
		exit(-1);
	}

	// get dataset information structures
	if(vflg == TRUE) {
		printf("S A M   E V E N T - R E L A T E D   S O U R C E   I M A G I N G\t%s\n", PRG_REV);
		printf("opening data file");
		fflush(stdout);
	}

#if BTI
	// get data with sensor structures
	if(GetMagnesInfo(DSName, PDFName, &Header, &Channel, &Epoch) == -1)
		Cleanup("can't read Magnes information");
	sprintf(DSpath, "%s/%s", Header.DsPath, Header.SetName);
#else
	// get data with sensor structures
	GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);
	sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
#endif
	sprintf(SAMpath, "%s/SAM", DSpath);

	// extract constants
	M = Header.NumPri;
	E = Header.NumEpochs;
	T = Header.MaxSamples;
	TT = E * T;
	f = Header.PsIndex[0];		// 1st MEG primary sensor channel index

	// parse analysis parameter specification file in current working directory
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("parsing '%s' parameter file", ParmName);
		fflush(stdout);
	}
	sprintf(fpath, "%s.param", ParmName);
	GetParams(fpath, &Params);
	if(vflg == TRUE) {
		printf(" - done\n");
		switch(Params.ImageMetric) {
			case MOMENT: printf("imaging averaged signal"); break;
			case POWER: printf("imaging averaged power"); break;
			case RV_ENTROPY: printf("imaging averaged rank vector entropy"); break;
			case RVC_ENTROPY: printf("imaging averaged conditional rank vector entropy"); break;
			default: Cleanup("illegal image metric"); break;
		}
		fflush(stdout);
	}
	if(Params.SignSegStart != -999 && Params.SignSegEnd != -999. && Params.ImageMetric == MOMENT) {	// determine polarity for moment, only!
		if(vflg == TRUE) {
			printf(" with polarity correction\n");
			fflush(stdout);
		}
		pflg = TRUE;
	} else {
		if(vflg == TRUE) {
			printf("\n");
			fflush(stdout);
		}
	}

	// sanity checks
	if(Params.BaseStart != -999. && Params.BaseEnd != -999.)
		bflg = TRUE;
	if(Params.CovHP == -999. || Params.CovLP == -999.)
		Cleanup("CovBand parameters req'd");
	if(Params.ImageHP == -999. || Params.ImageLP == -999.)
		Cleanup("ImageBand parameters req'd");
	if(Params.ImageHP < Channel[f].AnalogHpFreq)
		Cleanup("ImageBand highpass cutoff frequency out of recorded range");
	if(Params.ImageLP > Channel[f].AnalogLpFreq)
		Cleanup("ImageBand lowpass cutoff frequency out of recorded range");
	BWratio = (Params.ImageLP - Params.ImageHP) / (Params.CovLP - Params.CovHP);
	if(BWratio > 1.)
		Cleanup("covariance bandwidth must be >= image bandwidth");
	if((Params.ImageMetric == POWER || Params.ImageMetric == MOMENT) && (Params.SmoothHP == -999. || Params.SmoothLP == -999.))
		Cleanup("SmoothBand parameters req'd");
	switch(Params.NumMark) {
		case 1:
			if(!strncmp("NULL", Params.Marker[0].MarkName, 4))
				Cleanup("Marker1 name req'd");
			break;
		case 2:
			if(!strncmp("NULL", Params.Marker[0].MarkName, 4) || !strncmp("NULL", Params.Marker[1].MarkName, 4))
				Cleanup("both Marker1 & Marker2 names req'd");
			if(Params.Marker[0].MarkStart != Params.Marker[1].MarkStart || Params.Marker[0].MarkEnd != Params.Marker[1].MarkEnd)
				Cleanup("Marker1 & Marker2 must have identical start & end times");
			break;
		default:
			break;
	}
	if(Params.ImageMetric == RV_ENTROPY)
		if(Params.Dims < 3 || Params.Dims > 6)
			Cleanup("Rank vector entropy embedding dimension must be 3 to 6");
	if(Params.ImageMetric == RVC_ENTROPY)
		if(Params.Lags < 1. || Params.Lags > 100.)
			Cleanup("Conditional entropy advance must be from 1-100 ms");
	if(Params.FilterType == IIR)
		strcpy(FilterName, "IIR");
	else
		strcpy(FilterName, "FFT");

	// set prefix to NumPrefix characters (default is 8-character dataset hashcode)
	memset(Prefix, 0, 64);
#if BTI
	memcpy(Prefix, "BIOMAG4D", 8);				// 4D datasets just have serial numbers
#else
	memcpy(Prefix, DSName, Params.NumPrefix);	// default: copy the eight-letter NIMH hash code in dataset name
#endif
	memset(Suffix, 0, 10);
	switch(Params.ImageMetric) {
		case MOMENT: memcpy(Suffix, "MOM", 3); break;
		case POWER: memcpy(Suffix, "PWR", 3); break;
		case RV_ENTROPY: memcpy(Suffix, "RVE", 3); break;
		case RVC_ENTROPY: memcpy(Suffix, "RVC", 3); break;
		default: break;
	}

	// read markers
#if BTI                                                         // parse trigger channel for markers
	if(vflg == TRUE) {
		printf("parsing trigger markers");
		fflush(stdout);
	}
	if(GetMagnesMarkers(DSpath, PDFName, &Header, Channel, &Marker, &NumMarkers) == -1)
		Cleanup("cannot get trigger markers");
#else                                                           // try to read 'MarkerFile.mrk'
	if(vflg == TRUE) {
		printf("reading 'MarkerFile.mrk'");
		fflush(stdout);
	}
	sprintf(fpath, "%s/MarkerFile.mrk", DSpath);
	GetMarkers(fpath, &Marker, &NumMarkers);
#endif

	// summarize the number of epochs matching the marker name
	for(n=NumTriggers=0; n<NumMarkers; n++) {
		if(!strcmp(Marker[n].Name, Params.Marker[0].MarkName))
			NumTriggers++;
		if(Params.NumMark == 2)
			if(!strcmp(Marker[n].Name, Params.Marker[1].MarkName))
				NumTriggers++;
	}
	if(NumTriggers == 0)
		Cleanup("specified trigger marker(s) not found -- use valid trigger marker name");

	// set average & baseline window sizes
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("%d markers matched the designated trigger marker(s)\n", NumTriggers);
		printf("parsing time windows");
		fflush(stdout);
	}
	time = Params.Marker[0].MarkEnd - Params.Marker[0].MarkStart;		// duration of average window (seconds)
	TA = (int)rint(time * Header.SampleRate);							// number of samples in average window
	ta = (double)TA;
	NumImg = (Params.TimeStep > 0.)? time / Params.TimeStep: time * Header.SampleRate;
	AS = (int)rint(Params.Marker[0].MarkStart * Header.SampleRate);		// average window starting sample index, relative to marker
	AE = AS + TA;														// average window ending sample index
	if(bflg == TRUE) {
		time = Params.BaseEnd - Params.BaseStart;						// duration of baseline window (seconds)
		TB = (int)rint(time * Header.SampleRate);						// number of samples in baseline window
		tb = (double)TB;
		BS = (int)rint(Params.BaseStart * Header.SampleRate);			// baseline starting sample relative to trigger marker (can be negative!)
		BE = BS + TB;													// baseline ending sample relative to trigger marker (can be negative!)
		if(BS < AS || BS > AE || BE < AS || BE > AE)
			Cleanup("baseline starting & ending samples out of range");
	}
	if(pflg == TRUE) {													// optionally, determine time window for polarity
		time = Params.SignSegEnd - Params.SignSegStart;					// duration of sign segment window (seconds)
		TP = (int)rint(time * Header.SampleRate);						// number of samples in polarity window
		tp = (double)TP;
		PS = (int)rint(Params.SignSegStart * Header.SampleRate);
		PE = PS + TP;
		if(PS < AS || PS > AE || PE < AS || PE > AE)
			Cleanup("polarity segment starting & ending samples out of range");
	}

	// read trigger markers & parse
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("parsing trials");
		fflush(stdout);
	}
	if((Trial = (AVG_TRIAL *)malloc((size_t)(NumTriggers * sizeof(AVG_TRIAL)))) == NULL)
		Cleanup("'Trial[N]' allocation failed");
	for(n=k=NumTrials1=NumTrials2=0; n<NumMarkers; n++) {
		if(!strcmp(Marker[n].Name, Params.Marker[0].MarkName)) {
			sample = (int)(Marker[n].Latency * Header.SampleRate);	// sample for this trigger marker
			Trial[k].Epoch = Marker[n].Epoch;
			Trial[k].TS = sample + AS;								// offset trial start time relative to sample marker
			if(Trial[k].TS >= 0 || (Trial[k].TS + TA) < T) {		// increment if this trial is within the epoch samples
				Trial[k].Marker = 1;								// 1st marker used for ERS mode
				NumTrials1++;
				k++;
			}
		}
		if(Params.NumMark == 2) {
			if(!strcmp(Marker[n].Name, Params.Marker[1].MarkName)) {	// for DERS, read control marker name
				sample = (int)(Marker[n].Latency * Header.SampleRate);	// sample for this trigger marker
				Trial[k].Epoch = Marker[n].Epoch;
				Trial[k].TS = sample + AS;							// offset trial start time relative to sample marker
				if(Trial[k].TS >= 0 || (Trial[k].TS + TA) < T) {	// increment if this trial is within the epoch samples
					Trial[k].Marker = 2;							// 2nd marker used for DERS mode
					NumTrials2++;
					k++;
				}
			}
		}
	}
	TotTrials = k;
	if(TotTrials < 4)
		Cleanup("insufficient trials - check trigger & time windows");
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("found %d trials for Marker1", NumTrials1);
		if(Params.NumMark == 2)
			printf(" & %d trials for Marker2", NumTrials2);
		printf(" (total trials %d)\n", k);
		fflush(stdout);
	}

	// revise start & end times relative to start of average so that we can deal with averaged data...
	BS -= AS;
	BE -= AS;
	PS -= AS;
	PE -= AS;
	AE -= AS;
	AS = 0;							// average window now starts at sample zero!
	if(vflg == TRUE) {
		printf("average: sample %d to %d (%d samples, %5.3f sec)\n", AS, AE, TA, ta/Header.SampleRate);
		if(bflg == TRUE)
			printf("baseline: sample %d to %d (%d samples, %5.3f sec)\n", BS, BE, TB, tb/Header.SampleRate);
		if(pflg == TRUE)
			printf("sign: sample %d to %d (%d samples, %5.3f sec)\n", PS, PE, TP, tp/Header.SampleRate);
		fflush(stdout);
	}

	// weight name must be extend using the params...
	if(vflg == TRUE) {
		printf("reading ");
		fflush(stdout);
	}
	sprintf(WgtDir, "%s/%s,%-d-%-dHz", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
	switch(Params.CovType) {
		case GLOBAL_:	// 'Global' weights
			if(vflg == TRUE) {
				printf("'Global'");
				fflush(stdout);
			}
			sprintf(fpath, "%s/Global_Noise", WgtDir);
			if((np = fopen(fpath, "r")) == NULL)
				Cleanup("can't open noise file");
			if(fscanf(np, "%le", &sigma2) != 1)
				Cleanup("can't read noise file");
			fclose(np);
			if(Params.ImageFormat == TLRC)
				sprintf(fpath, "%s/Global_at.nii", WgtDir);
			else
				sprintf(fpath, "%s/Global.nii", WgtDir);
			Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
			V = Wgt->size1;
			if(Wgt->size2 != M)
				Cleanup("number of channels in dataset does not match number in weight file");
			break;
		case SUM_:		// 'Sum' weights
			if(vflg == TRUE) {
				printf("'Sum'");
				fflush(stdout);
			}
			sprintf(fpath, "%s/Sum_Noise", WgtDir);
			if((np = fopen(fpath, "r")) == NULL)
				Cleanup("can't open noise file");
			if(fscanf(np, "%le", &sigma2) != 1)
				Cleanup("can't read noise file");
			fclose(np);
			if(Params.ImageFormat == TLRC)
				sprintf(fpath, "%s/Sum_at.nii", WgtDir);
			else
				sprintf(fpath, "%s/Sum.nii", WgtDir);
			Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
			V = Wgt->size1;
			if(Wgt->size2 != M)
				Cleanup("number of channels in dataset does not match number in weight file");
			break;
		default:
			Cleanup("CovType must be either GLOBAL or SUM");
			break;
	}
	if(vflg == TRUE) {
		printf(" SAM weights");
		fflush(stdout);
	}

	// allocate memory for arrays
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("allocating memory");
		fflush(stdout);
	}
	W = gsl_vector_alloc(M);
	if((x = (double ***)malloc((size_t)E * sizeof(double **))) == NULL)
		allocfailed("x[]");
	for(e=0; e<E; e++) {
		if((x[e] = (double **)malloc((size_t)T * sizeof(double *))) == NULL)
			allocfailed("x[][]");
		for(t=0; t<T; t++)
			if((x[e][t] = (double *)malloc((size_t)M * sizeof(double))) == NULL)
				allocfailed("x[][][]");
	}
	if((In = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
			allocfailed("In[]");
	if((Out = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
			allocfailed("Out[]");
	if((avg0 = (double *)malloc((size_t)TA * sizeof(double))) == NULL)
		allocfailed("avg0[]");
	if((avg1 = (double *)malloc((size_t)TA * sizeof(double))) == NULL)
		allocfailed("avg1[]");
	if((avg2 = (double *)malloc((size_t)TA * sizeof(double))) == NULL)
		allocfailed("avg2[]");
	if((Image = (float **)malloc((size_t)NumImg * sizeof(float *))) == NULL)
		allocfailed("Image[]");
	for(d=0; d<NumImg; d++)
		if((Image[d] = (float *)malloc((size_t)V * sizeof(float))) == NULL)
			allocfailed("Image[][]");

	// determine filter type, & determine band-limited noise
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("making %6.3f to %6.3f Hz %s filter", Params.ImageHP, Params.ImageLP, FilterName);
		fflush(stdout);
	}
	if(Params.FilterType == IIR)
		mkfilter(&DataIIR, SIG_ORDER, Header.SampleRate, Params.ImageHP, Params.ImageLP, &BWFreq);
	else
		Butterworth(&DataFFT, T, Params.ImageHP, Params.ImageLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

	// if req'd, compute lowpass smoothing filter
	if(Params.ImageMetric == POWER || Params.ImageMetric == MOMENT) {
		if(vflg == TRUE) {
			printf(" - done\n");
			printf("making %6.3f Hz lowpass %s filter", Params.SmoothLP, FilterName);
			fflush(stdout);
		}
		if(Params.FilterType == IIR)
			mkfilter(&SmoothIIR, SIG_ORDER, Header.SampleRate, 0., Params.SmoothLP, &BWFreq);
		else
			Butterworth(&SmoothFFT, TT, 0., Params.SmoothLP, Header.SampleRate, LOWPASS, MAXORDER, &BWFreq);
		}

	// compute notch filter
	if(Params.Notch == TRUE) {
		if(Params.FilterType == FFT) {
			if(vflg == TRUE) {
				printf(" - done\n");
				printf("making FFT notch filter");
				fflush(stdout);
			}
			mknotch(&NotchFFT, Params.CovHP, Params.CovLP, Header.SampleRate, T);
		} else {
			fprintf(stderr, "\nreminder: when using IIR filters, data should be notch filtered for power mains prior to analysis\n");
		}
	}


	// read & filter MEG data
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("reading & filtering MEG data");
		fflush(stdout);
	}
	for(e=0; e<E; e++) {
		for(m=0; m<M; m++) {
#if BTI
			if(GetMagnesChan(DSpath, PDFName, &Header, Channel, e, m, In) == -1)
				Cleanup("'GetMagnesChan()' failed");
#else
			n = Header.PsIndex[m];
			scale = Channel[n].Scale;
			GetDsData(&Header, e, n, scale, In);
#endif
			demean(In, T);
			if(Params.FilterType == IIR) {
				bdiir(In, Out, T, &DataIIR);
				for(t=0; t<T; t++)
					x[e][t][m] = Out[t];
			} else {
				FFTfilter(In, Out, &DataFFT, T);
				if(Params.Notch == TRUE) {
					FFTfilter(Out, In, &NotchFFT, T);
					for(t=0; t<T; t++)
						x[e][t][m] = In[t];
				} else {
					for(t=0; t<T; t++)
						x[e][t][m] = Out[t];
				}
			}
		}
		if(vflg == TRUE) {			// if req'd - heartbeat
			printf(".");
			fflush(stdout);
		}
	}

	// for each voxel...
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("computing event-related image voxels:\n");
		fflush(stdout);
	}
	for(c=0, pct=(-1); c<V; c++) {

		// get voxel weight vector & determine if within head boundary
		gsl_matrix_get_row(W, Wgt, c);
		if(gsl_vector_isnull(W) == 0) {		// compute image only for non-zero weights

			// compute band-limited noise of this virtual sensor
			for(m=0, tmp=0.; m<M; m++)		// compute noise of unaveraged primary sensors
				tmp += gsl_vector_get(W, m) * sigma2 * gsl_vector_get(W, m);
			tmp *= BWratio;
			rms_noise = sqrt(tmp);			// rms noise estimate

			// compute voxel time-series for entire dataset
			for(e=tt=0; e<E; e++)
				for(t=0; t<T; t++, tt++) {
					for(m=0, tmp=0.; m<M; m++)
						tmp += x[e][t][m] * gsl_vector_get(W, m);
					In[tt] = tmp / rms_noise;
				}

			// compute metric for entire voxel time-series
			switch(Params.ImageMetric) {
				case MOMENT:			// event-related moment (normalized S/N)
					if(Params.FilterType == IIR)
						bdiir(In, Out, TT, &SmoothIIR);				// lowpass to smooth moment
					else
						FFTfilter(In, Out, &SmoothFFT, TT);			// lowpass to smooth moment
					break;
				case POWER:				// Hilbert envelope of power (normalized S/N)
					hilbert(TT, In, Out);
					for(t=0; t<TT; t++)								// sum the square of the in-phase & quadrature signal components
						In[t] = 0.5 * (In[t] * In[t] + Out[t] * Out[t]);
					if(Params.FilterType == IIR)					// lowpass to smooth Hilbert envelope
						bdiir(In, Out, TT, &SmoothIIR);
					else
						FFTfilter(In, Out, &SmoothFFT, TT);
					break;
				case RV_ENTROPY:		// rank vector entropy
					Data2RVE(In, Out, TT, Params.Dims, Params.ImageLP, Params.Tau, Header.SampleRate);
					break;
				case RVC_ENTROPY:		// rank vector conditional entropy
					Data2CondRVE(In, Out, TT, Params.Dims, Params.ImageLP, Params.Tau, Params.Lags, Header.SampleRate);
					break;
				default:
					Cleanup("ImageMetric not implemented");
					break;
			}

			// zero averages
			for(t=0; t<TA; t++) {
				avg1[t] = 0.;
				avg2[t] = 0.;
			}

			// for each trial...
			for(k=0, na1=na2=0; k<TotTrials; k++) {

				// set up segment initialization
				e = Trial[k].Epoch;
				ts = e * T + Trial[k].TS;

				// update average(s)
				switch(Trial[k].Marker) {
					case 1:		// marker 1
						for(t=0, tt=ts; t<TA; t++, tt++)
							avg1[t] += Out[tt];
						na1++;
						break;
					case 2:		// marker 2
						for(t=0, tt=ts; t<TA; t++, tt++)
							avg2[t] += Out[tt];
						na2++;
						break;
					default:
						break;
				}
			}

			// compute mean values
			if(na1 > 0)
				for(t=0; t<TA; t++)
					avg1[t] /= (double)na1;
			if(na2 > 0)
				for(t=0; t<TA; t++)
					avg2[t] /= (double)na2;

			// compute difference
			for(t=0; t<TA; t++)
				avg0[t] = avg1[t] - avg2[t];			// simple difference (if only 1 marker, it still works)

			// if req'd multiply RVE by 1000
			if(Params.ImageMetric == RV_ENTROPY || Params.ImageMetric == RVC_ENTROPY) {
				for(t=0; t<TA; t++)
					avg0[t] *= 1000.;
			}

			// if req'd, compute & subtract baseline
			if(bflg == TRUE) {
				for(t=BS, tmp=0.; t<BE; t++)
					tmp += avg0[t];
				baseline = tmp / tb;
				for(t=0; t<TA; t++)
					avg0[t] -= baseline;
			}

			// if req'd, estimate polarity
			if(strncmp("NULL", Params.Atlas, 4) || pflg == TRUE) {
				for(t=PS, tmp=0.; t<PE; t++)
					tmp += avg0[t];
				sign = (tmp > 0.)? 1.: -1.;
			} else {
				sign = 1.;
			}

			// for each latency...
			for(d=0, time=0.; d<NumImg; d++, time+=Params.TimeStep) {
				t = (int)rint(time * Header.SampleRate);
				if(aflg == TRUE)
					Image[d][c] = (float)fabs(avg0[t]);
				else
					Image[d][c] = (float)(avg0[t] * sign);
			}
		} else {
			for(d=0; d<NumImg; d++)
				Image[d][c] = 0.;
		}

		// percent done
		if(vflg == TRUE) {
			percent = 100. * (double)c / (double)V;
			if(rint(percent) > (double)pct) {
				pct = (int)rint(percent);
				printf("\r\t%d percent done", pct);
				fflush(stdout);
			}
		}
	}

	// open output SAM image file, & write header
	if(vflg == TRUE) {
		printf("\nSAM(ers) image done\n");
		printf("writing SAM(ers) images");
		fflush(stdout);
	}
	if(!strncmp(Params.DirName, "NULL", 4))
		sprintf(ImgPath, "%s/%s,%s,", SAMpath, Prefix, ParmName);
	else
		sprintf(ImgPath, "%s/%s,%s,", Params.DirName, Prefix, ParmName);
	switch(Params.NumMark) {
		case 1: sprintf(fpath, "%s%s,%s,ERSC.nii", ImgPath, Params.Marker[0].MarkName, Suffix); break;
		case 2: sprintf(fpath, "%s%s-%s,%s,dERSC.nii", ImgPath, Params.Marker[0].MarkName, Params.Marker[1].MarkName, Suffix); break;
		default: break;
	}
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open .nii file for write");

	// modify NIFTI header
	memset(extension, 0, 4);									// for now, ignore the fact that we may have an extended header
	NiiHdr.dim[0] = 4;											// 3d+time
	NiiHdr.dim[4] = NumImg;										// number of 3d images
	NiiHdr.toffset = (float)Params.Marker[0].MarkStart;			// time offset shift (markers 1 & 2 must be identical)
	NiiHdr.slice_duration = (float)Params.TimeStep;				// time step size (s)
	NiiHdr.pixdim[4] = (float)Params.TimeStep;
	NiiHdr.slice_code = NIFTI_SLICE_SEQ_INC;
	NiiHdr.xyzt_units = (NIFTI_UNITS_MM | NIFTI_UNITS_SEC);		// 3d+time units
	memset(NiiHdr.intent_name, 0, 16);
	switch(Params.NumMark) {
		case 1: memcpy(NiiHdr.intent_name, "ERS", 3); break;
		case 2: memcpy(NiiHdr.intent_name, "Differential ERS", 16); break;
		default: break;
	}

	// write header & extension
	if(fwrite((void *)&NiiHdr, sizeof(nifti_1_header), 1, fp) != 1)
		Cleanup("can't write nifti header to .nii file");
	if(fwrite((void *)extension, 4, 1, fp) != 1)
		Cleanup("can't write 'extension' to .nii file");
	if((int)NiiHdr.vox_offset != 352 && ExtHdr != NULL) {
		N = (int)NiiHdr.vox_offset - 352;
		if(fwrite((void *)ExtHdr, N, 1, fp) != 1)
			Cleanup("can't write extended header to .nii file");
	}
	for(d=0; d<NumImg; d++)		// write time as BRIKs
		if(fwrite((void *)Image[d], sizeof(float), V, fp) != V)
			Cleanup("can't write data image file");
	fclose(fp);

	// we're done now!
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("'SAMersc' done\n");
		fflush(stdout);
	}
	exit(0);
}
