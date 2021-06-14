// SAMhfo -- synthetic aperture magnetometry for imaging epiletic activity.
//	SAMhfo reads SAM weights for volume image, & computes the source time-series
//	for each voxel.
//
//	A bandpass filter for HFOs is applied to the data. For each voxel, the virtual
//	sensor time-series is computed & a Hilbert envelope applied to the HFO source.
//	The excess kurtosis of the Hilbert envelope is used to localize "events" such
//	as short bursts of HFOs.
//
//	flags:
//	-r <run name>		-- MEG run name
//	-d <data name>		-- MEG data file name (4D, only)
//	-m <name>			-- parameter file name
//	-v					-- verbose mode
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <samUtils.h>
#include <filters.h>
#include <nifti1.h>
#include <version.h>

#define MINOR_REV	5
#define TRUE		1
#define FALSE		0
#define SIG_ORDER	4


int	main(
	int		argc,
	char	**argv
)

{
	IIRSPEC			HFO_IIR;		// IIR HFO data filter structure
	IIRSPEC			SmoothIIR;		// IIR smoothing filter structure
	FFTSPEC			HFO_FFT;		// HFOFilter -- FFT filter structure for HFOs
	FFTSPEC			SmoothFFT;		// Smooth -- FFT filter structure for smoothing Hilbert envelope
	PARMINFO		Params;			// analysis parameters
	HeaderInfo		Header;			// MEG data header
	ChannelInfo		*Channel;		// MEG channel info
	EpochInfo		*Epoch;			// MEG epoch info
	nifti_1_header	NiiHdr;			// NIFTI header
	unsigned char	*ExtHdr;		// NiiHdr[N] -- extended NIFTI header
	gsl_matrix		*Wgt;			// Wgt - VxM SAM coefficients
	gsl_vector		*W;				// W - Mx1 SAM coefficient vector
	double			**x;			// x[TT][M] -- multi-sensor signal
	double			*In;			// In[TT] -- input time-series
	double			*Out;			// Out[TT] -- output time-series
	double			p2;				// p^2
	double			p4;				// p^4
	double			var;			// 2nd moment (variance)
	double			sx2;			// sum of 2nd power
	double			sx4;			// sum of 4th power
	double			BWFreq;			// filter bandwidth
	double			tw;				// total samples in Window
	double			scale;			// bits-->T conversion
	double			percent;		// percent done
	register double	tmp;			// what else?
	float			*Image;			// Image[V] -- SAM(g2) image
	int				E;				// number of good epochs
	int				M;				// number of MEG or reference channels used
	int				N;				// number of characters in extended NIFTI header
	int				T;				// number of samples
	int				TT;				// total sample size
	int				D;				// offset from start & end of envelope (to avoid Hilbert artifacts)
	int				V;				// number of voxels
	int				c;				// input character & virtual channel
	int				pct;			// integer percentage
	int				m;				// channel index
	int				n;				// alternate channel index
	int				e;				// epoch index
	int				f;				// 1st primary channel index
	int				t;				// sample index
	int				tt;				// alternate time index
	int				stride;			// sample stride for decimating Hilbert envelope
	static int		eflg = FALSE;	// command-line error flag
	static int		rflg = FALSE;	// dataset name flag
	static int		vflg = FALSE;	// verbose mode flag
	static int		wflg = FALSE;	// SAM weight flag
	extern char		*optarg;
	extern int		opterr;
	int				LongIndex;		// position of arguments
	char			fpath[256];		// general path name
	char			DSName[256];	// MEG run name
	char			DSpath[256];	// MEG path name
	char			SAMpath[256];	// SAM subdirectory path
	char			ImgPath[256];	// output path for writing image
	char			ParmName[256];	// parameter file name
	char			Prefix[64];		// image file prefix
	char			FilterName[8];	// filter name -- IIR | FFT
	char			extension[4];	// extension
#if BTI
	static int		dflg = FALSE;	// pdf file error flag
	char			PDFName[256];	// pdf file name
	char			ShortOpts[] = "r:d:m:v";
	static struct option	LongOpts[] = {	// command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"pdf_name", required_argument, NULL, 'd'},
		{"parameter_name", required_argument, NULL, 'm'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};
#else
	char			ShortOpts[] = "r:m:v";
	static struct option	LongOpts[] = {	// command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"parameter_name", required_argument, NULL, 'm'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};
	unsigned char	**Bad;			// Bad[E][T] -- flags for bad samples
	unsigned char	*Flags;			// Flags[TT] -- transformed flags
#endif
	FILE			*fp;			// output file pointer (to SAM image)

	// parse command line parameters
	opterr = 1;		// enable error reporting
	while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
		switch(c) {
			case 'r':	// run name
				sprintf(DSName, "%s", optarg);
				rflg = TRUE;
				break;
#if BTI
			case 'd':	// data file name
				sprintf(PDFName, "%s", optarg);
				dflg = TRUE;
				break;
#endif
			case 'm':	// parameter file name
				sprintf(ParmName, "%s", optarg);
				wflg = TRUE;
				break;
			case 'v':	// verbose mode
				vflg = TRUE;
				break;
			default:
				eflg = TRUE;
				break;
		}
	}
#if BTI
	if(eflg == TRUE || rflg == FALSE || dflg == FALSE || wflg == FALSE) {
		fprintf(stderr, "SAMhfo\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "mandatory arguments:\n");
		fprintf(stderr, "\t-r <run name>\n");
		fprintf(stderr, "\t-d <data file name>\n");
#else
	if(eflg == TRUE || rflg == FALSE || wflg == FALSE) {
		fprintf(stderr, "SAMhfo\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "mandatory arguments:\n");
		fprintf(stderr, "\t-r <run name>\n");
#endif
		fprintf(stderr, "\t-m <parameter file name>\n");
		fprintf(stderr, "\t-v -- verbose mode\n");
		exit(-1);
	}

	// get dataset information structures
	if(vflg == TRUE) {
		printf("S A M   H F O   I M A G I N G\t%s\n", PRG_REV);
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

	// establish epoch, channel, & sample count
	E = Header.NumEpochs;
	T = Header.MaxSamples;
	TT = E * T;
	M = Header.NumPri;
	f = Header.PsIndex[0];		// 1st MEG primary sensor channel index

	// check to see if this is a continuous dataset
	if(E != 1)
		fprintf(stderr, " (reminder: multi-epoch data must be contiguous)");

	// parse analysis parameter specification file in current working directory
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("parsing '%s' parameter file", ParmName);
		fflush(stdout);
	}
	sprintf(fpath, "%s.param", ParmName);
	GetParams(fpath, &Params);

	// check if this is the right program for you!
	if(Params.CovHP == -999. || Params.CovLP == -999.)
		Cleanup("CovBand specifications required");
	if(Params.ImageHP == -999. || Params.ImageLP == -999.)
		Cleanup("ImageBand specifications required");
	if(Params.ImageHP < Channel[0].AnalogHpFreq)
		Cleanup("ImageBand highpass cutoff frequency out of recorded range");
	if(Params.ImageLP > Channel[0].AnalogLpFreq)
		Cleanup("ImageBand lowpass cutoff frequency out of recorded range");
	if(Params.SmoothHP == -999. || Params.SmoothLP == -999.)
		Cleanup("SmoothBand specifications required");
	if(Params.SmoothHP < Channel[f].AnalogHpFreq)
		Cleanup("SmoothBand highpass cutoff frequency out of recorded range");
	if(Params.SmoothLP > Channel[f].AnalogLpFreq)
		Cleanup("SmoothBand lowpass cutoff frequency out of recorded range");
	if(Params.FilterType == IIR)
		strcpy(FilterName, "IIR");
	else
		strcpy(FilterName, "FFT");
	if(Params.ImageMetric != KURTOSIS)
		Cleanup("ImageMode must be 'Kurtosis'");
	if(Params.CovType != GLOBAL_)
		Cleanup("CovType must be GLOBAL");
	if(strncmp(Params.DirName, "NULL", 4))
		sprintf(ImgPath, "%s", Params.DirName);
	else
		sprintf(ImgPath, "%s", SAMpath);

	// set prefix to NumPrefix characters (default is 8-character dataset hashcode)
	memset(Prefix, 0, 64);
#if BTI
	memcpy(Prefix, "BIOMAG4D", 8);				// 4D datasets just have serial numbers
#else
	memcpy(Prefix, DSName, Params.NumPrefix);

	// copy flags
	if((Flags = (unsigned char *)malloc((size_t)TT * sizeof(unsigned char))) == NULL)
		allocfailed("Flags[]");
	for(e=0; e<E; e++)
		for(t=0; t<T; t++) {
			tt = e * T + t;
			Flags[tt] = Bad[e][t];
		}
#endif

	// get weight name from parameters & read SAM weights
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("reading Global SAM weights");
		fflush(stdout);
	}
	sprintf(fpath, "%s/%s,%-d-%-dHz/Global.nii", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
	Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
	V = Wgt->size1;
	if(Wgt->size2 != M)
		Cleanup("number of channels in dataset does not match number in weight file");

	// set up imaging constants
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("determining imaging constants:\n");
		fflush(stdout);
	}
	stride = (int)floor(Header.SampleRate / (2. * Params.SmoothLP));
	D = 3 * (int)floor(Header.SampleRate / Params.ImageHP);
	if(vflg == TRUE) {
		printf("\tstride=%d samples\n", stride);
		printf("\toffset=%d samples\n", D);
		printf("\tepochs=%d, samples/epoch=%d (total %d)\n", E, T, TT);
		fflush(stdout);
	}

	// allocate memory for arrays
	if(vflg == TRUE) {
		printf("allocating memory");
		fflush(stdout);
	}
	W = gsl_vector_alloc(M);
	if((x = (double **)malloc((size_t)TT * sizeof(double *))) == NULL)
		allocfailed("x[]");
	for(t=0; t<TT; t++)
		if((x[t] = (double *)malloc((size_t)M * sizeof(double))) == NULL)
				allocfailed("x[][]");
	if((In = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
		allocfailed("In[]");
	if((Out = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
		allocfailed("Out[]");
	if((Image = (float *)malloc((size_t)V * sizeof(float))) == NULL)
		allocfailed("Image[]");

	// generate HFO filter
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("making %6.3f to %6.3f Hz %s filter for HFO bandpass", Params.ImageHP, Params.ImageLP, FilterName);
		fflush(stdout);
	}
	if(Params.FilterType == IIR)
		mkfilter(&HFO_IIR, SIG_ORDER, Header.SampleRate, Params.ImageHP, Params.ImageLP, &BWFreq);
	else
		Butterworth(&HFO_FFT, TT, Params.ImageHP, Params.ImageLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

	// generate smoothing filter
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("making %6.3f Hz lowpass %s filter for smoothing Hilbert envelope", Params.SmoothLP, FilterName);
		fflush(stdout);
	}
	if(Params.FilterType == IIR)
		mkfilter(&SmoothIIR, SIG_ORDER, Header.SampleRate, 0., Params.SmoothLP, &BWFreq);
	else
		Butterworth(&SmoothFFT, TT, 0., Params.SmoothLP, Header.SampleRate, LOWPASS, MAXORDER, &BWFreq);

	// read & filter MEG data
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("reading & filtering MEG data");
		fflush(stdout);
	}
	for(m=0; m<M; m++) {							// read & process data in channel order
		n = Header.PsIndex[m];
		scale = Channel[n].Scale;
		for(e=0; e<E; e++) {						// multi-epoch data is presumed to be contiguous!
#if BTI
			if(GetMagnesChan(DSpath, PDFName, &Header, Channel, e, m, In) == -1)
				Cleanup("'GetMagnesChan()' failed");
#else
			GetDsData(&Header, e, n, scale, In);	// get one epoch of data for this channel -- length T samples
#endif
			// make continuous data out of multi-epochs
			for(t=0; t<T; t++) {
				tt = e * T + t;
				Out[tt] = In[t];
			}
		}
		demean(Out, TT);							// remove the mean
		if(Params.FilterType == IIR)
			bdiir(Out, In, TT, &HFO_IIR);			// now, filter the 'Out' data
		else
			FFTfilter(Out, In, &HFO_FFT, TT);		// now, filter the 'Out' data
		for(t=0; t<TT; t++)							// move to 'x[][]' array
			x[t][m] = In[t];

		// heartbeat
		if(vflg == TRUE) {
			printf(".");
			fflush(stdout);
		}
	}

	// now, step through all coordinates
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("stepping through voxels to localize HFOs:\n");
		fflush(stdout);
	}
	for(c=0, pct=(-1); c<V; c++) {
		gsl_matrix_get_row(W, Wgt, c);
		if(gsl_vector_isnull(W) == 0) {

			// (1) compute virtual sensor time series
			for(t=0; t<TT; t++) {
				for(m=0, tmp=0.; m<M; m++)
					tmp += gsl_vector_get(W, m) * x[t][m];
				In[t] = tmp;
			}

			// (2) compute smoothed Hilbert envelope of HFOs
			hilbert(TT, In, Out);
			for(t=0; t<TT; t++)
				Out[t] = (In[t] * In[t] + Out[t] * Out[t]);
			if(Params.FilterType == IIR)
				bdiir(Out, In, TT, &SmoothIIR);		// lowpass to smooth Hilbert envelope
			else
				FFTfilter(Out, In, &SmoothFFT, TT);	// lowpass to smooth Hilbert envelope

			// (3) compute g2 of smoothed Hilbert envelope (with decimation relative to filter lowpass)
#if BTI
			for(t=D, sx2=sx4=tw=0.; t<(TT-D); t+=stride) {
				p2 = In[t] * In[t];
				p4 = p2 * p2;
				sx2 += p2;
				sx4 += p4;
				tw += 1.;
			}
#else
			for(t=D, sx2=sx4=tw=0.; t<(TT-D); t+=stride)
				if(Bad[t] == GOOD) {
					p2 = In[t] * In[t];
					p4 = p2 * p2;
					sx2 += p2;
					sx4 += p4;
					tw += 1.;
				}
#endif
			if(tw > 12.) {
				var = sx2 / tw;
				Image[c] = sx4 / (tw * var * var) - 3.;
			} else {
				Image[c] = 0.;
			}
		} else {
			Image[c] = 0.;
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

	// write out SAM(g2) image
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing SAM(Epi) image");
		fflush(stdout);
	}
	sprintf(fpath, "%s/%s,%s,HFO.nii", ImgPath, Prefix, ParmName);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open .nii file for write");

	// modify NIFTI header
	memset(extension, 0, 4);				// for now, ignore the fact that we may have an extended header
	NiiHdr.dim[0] = 3;						// 3d volume
	NiiHdr.dim[4] = 0;						// number of 3d images
	NiiHdr.toffset = 0.;					// time offset shift
	NiiHdr.slice_duration = 1.;
	NiiHdr.slice_code = NIFTI_SLICE_UNKNOWN;	
	NiiHdr.xyzt_units = NIFTI_UNITS_MM;		// units
	memset(NiiHdr.intent_name, 0, 16);
	memcpy(NiiHdr.intent_name, "Epilepsy HFO", 12);

	if(fwrite((void *)&NiiHdr, sizeof(nifti_1_header), 1, fp) != 1)
		Cleanup("can't write nifti header to .nii file");
	if(fwrite((void *)extension, 4, 1, fp) != 1)
		Cleanup("can't write 'extension' to .nii file");
	if((int)NiiHdr.vox_offset != 352 && ExtHdr != NULL) {
		N = (int)NiiHdr.vox_offset - 352;
		if(fwrite((void *)ExtHdr, N, 1, fp) != 1)
			Cleanup("can't write extended header to .nii file");
	}
	if(fwrite((void *)Image, sizeof(float), V, fp) != V)
		Cleanup("can't write voxel values to .nii file");
	fclose(fp);

	// that's all, folks!
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("'SAMhfo' done\n");
		fflush(stdout);
	}
	exit(0);
}
