// SAMvs -- SAMvs reads SAM weights for a specified list of virtual sensors,
//	and computes the source time-series for each. The output is an ascii
//	file.
//
//	flags:
//	-r <run name>				-- MEG run name
//	-d <data file name>			-- MEG data file name (4D, only)
//	-m <parameter file name>	-- parameter file name
//	-v							-- verbose mode
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

#define MINOR_REV	0
#define SIG_ORDER	4

int
main(
	int			argc,
	char		**argv
)

{
	IIRSPEC		DataIIR;		// IIR filter structure
	FFTSPEC		DataFFT;		// FFT filter structure
	PARMINFO	Params;			// analysis parameters
	HeaderInfo	Header;			// MEG data header
	ChannelInfo	*Channel;		// MEG channel info
	EpochInfo	*Epoch;			// MEG epoch info
	nifti_1_header	NiiHdr;		// NIFTI header
	unsigned char	*ExtHdr;	// NiiHdr[N] -- extended NIFTI header
	gsl_matrix	*Wgt;			// Wgt - VxM SAM coefficients
	gsl_vector	*W;				// W - Mx1 SAM coefficient vector
	double		**x;			// x[TT][M] -- MEG signal
	double		*In;			// In[T] -- input time-series
	double		*Out;			// Out[T] -- output time-series
	double		BWFreq;			// filter bandwidth
	double		tmp;			// what else?
	float		*vs;			// vs[TT] -- virtual sensor time series
	float		xp;
	float		yp;
	float		zp;
	int			E;				// number of good epochs
	int			M;				// number of MEG or reference channels used
	int			T;				// number of samples
	int			TT;				// total samples E*T
	int			V;				// number of virtual sensor channels
	int			c;				// input character & virtual channel
	int			stride;			// output stride
	int			m;				// channel index
	int			e;				// epoch index
	int			f;				// 1st primary sensor index
	int			t;				// sample index
	int			tt;				// alternate sample index
	int			i;				// matrix index
	int			j;				// another index
	int			k;				// last index
	int			LongIndex;		// position of arguments
	static int	eflg = FALSE;	// command-line error flag
	static int	rflg = FALSE;	// run name flag
	static int	vflg = FALSE;	// verbose mode flag
	static int	mflg = FALSE;	// SAM weight flag
	extern char	*optarg;
	extern int	opterr;
	char		SAMpath[256];	// SAM path
	char		fpath[256];		// general path name
	char		DSName[256];	// MEG dataset name
	char		DSpath[256];	// dataset path
	char		ParmPath[256];	// parameter pathname
	char		*ParmName;		// name of parameter file
	char		ImgPath[256];	// output path for writing image
	char		Prefix[64];		// image file prefix
	char		extension[4];	// extension
#if BTI
	static int		dflg = FALSE;	// dataset flag
	char			ShortOpts[] = "r:d:m:v";
	static struct option	LongOpts[] = {	// command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"pdf_name", required_argument, NULL, 'd'},
		{"analysis_parameters", required_argument, NULL, 'm'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};
	char			PDFName[256];	// pdf file name
#else
	double			scale;		// bits-->Tesla conversion
	int				n;			// alternate sensor index
	unsigned char	**Bad;		// Bad[E][T] -- flags for bad samples
	char			ShortOpts[] = "r:m:v";
	static struct option	LongOpts[] = {	// command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"analysis_parameters", required_argument, NULL, 'm'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};
#endif	
	FILE		*fp;			// file pointer for reading time segment

	// parse command line parameters
	opterr = 1;			// enable error reporting
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
			case 'm':	// parameter pathname
				sprintf(ParmPath, "%s", optarg);
				mflg = TRUE;
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
	if(eflg == TRUE || rflg == FALSE || dflg == FALSE || mflg == FALSE) {
		fprintf(stderr, "SAMvs\t-r <run name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "\t-d <data file name>\n");
#else
		if(eflg == TRUE || rflg == FALSE || mflg == FALSE) {
			fprintf(stderr, "SAMvs\t-r <run name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
#endif
		fprintf(stderr, "\t-m <parameter file name>\n");
		fprintf(stderr, "\t-v -- verbose mode\n");
		exit(-1);
	}

	// get dataset information structures
	if(vflg == TRUE) {
		printf("S A M   V I R T U A L   S E N S O R   T O   T E X T\t%s\n", PRG_REV);
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

	// determine channels, epochs, & sample count
	M = Header.NumPri;
	E = Header.NumEpochs;
	T = Header.MaxSamples;
	TT = E * T;
	f = Header.PsIndex[0];		// 1st MEG primary sensor channel index

	// parse analysis parameter specification file in current working directory
	if(rindex(ParmPath, '/') == NULL)
		ParmName = ParmPath;
	else
		ParmName = rindex(ParmPath, '/') + 1;
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("parsing '%s' parameter file", ParmName);
		fflush(stdout);
	}
	sprintf(fpath, "%s.param", ParmPath);
	GetParams(fpath, &Params);

	// check if this is the right program for you!
	if(Params.CovHP == -999. || Params.CovLP == -999.)
		Cleanup("CovBand specifications required");
	if(Params.CovHP < Channel[f].AnalogHpFreq)
		Cleanup("CovBand highpass cutoff frequency out of recorded range");
	if(Params.CovLP > Channel[f].AnalogLpFreq)
		Cleanup("CovBand lowpass cutoff frequency out of recorded range");
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
#endif

	// set up imaging constants
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("determining stride:\n");
		fflush(stdout);
	}
	stride = (int)floor(Header.SampleRate / (2. * Params.CovLP));
	if(vflg == TRUE) {
		printf("\tstride=%d samples\n", stride);
		printf("\tepochs=%d, samples/epoch=%d (total %d)\n", E, T, TT);
		fflush(stdout);
	}

	// get weight name from parameters & read SAM weights
	if(vflg == TRUE) {
		printf("reading Global SAM weights");
		fflush(stdout);
	}
	if(Params.ImageFormat == TLRC)
		sprintf(fpath, "%s/%s,%-d-%-dHz/Global_at.nii", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
	else
		sprintf(fpath, "%s/%s,%-d-%-dHz/Global.nii", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
	Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
	V = Wgt->size1;
	if(Wgt->size2 != M)
		Cleanup("number of channels in dataset does not match number in weight file");

	// allocate memory for arrays
	if(vflg == TRUE) {
		printf(" - done\n");
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
	if((vs = (float *)malloc((size_t)TT * sizeof(float))) == NULL)
		allocfailed("vs[]");

	// compute filter
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("making %6.3f to %6.3f Hz filter", Params.CovHP, Params.CovLP);
		fflush(stdout);
	}
	if(Params.FilterType == IIR)
		mkfilter(&DataIIR, SIG_ORDER, Header.SampleRate, Params.CovHP, Params.CovLP, &BWFreq);
	else
	Butterworth(&DataFFT, TT, Params.CovHP, Params.CovLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

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
			bdiir(Out, In, TT, &DataIIR);			// now, filter the 'Out' data
		else
			FFTfilter(Out, In, &DataFFT, TT);		// now, filter the 'Out' data
		for(t=0; t<TT; t++)							// move to 'x[][]' array
			x[t][m] = In[t];

		// heartbeat
		if(vflg == TRUE) {
			printf(".");
			fflush(stdout);
		}
	}

	// open file for write
	sprintf(fpath, "%s/%s,%s,VS_Data", ImgPath, Prefix, ParmName);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open data file for write");

	// create virtual sensor data
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("generating virtual sensors");
		fflush(stdout);
	}
	for(c=0; c<V; c++) {
		gsl_matrix_get_row(W, Wgt, c);
		if(gsl_vector_isnull(W) == 0) {

			// compute virtual sensor time series
			for(t=i=0; t<TT; t+=stride, i++) {
				for(m=0, tmp=0.; m<M; m++)
					tmp += gsl_vector_get(W, m) * x[t][m];
				vs[i] = (float)tmp;
			}

			// write a row
			if(fwrite((void *)vs, sizeof(float), i, fp) != i)
				Cleanup("can't write to data values");
		}
	}
	fclose(fp);

	// open file for write
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing time track");
		fflush(stdout);
	}
	sprintf(fpath, "%s/%s,%s,VS_Time", ImgPath, Prefix, ParmName);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open time file for write");

	// write the time track
	for(t=i=0; t<TT; t+=stride, i++)
		vs[i] = (float)t / (float)Header.SampleRate;
	if(fwrite((void *)vs, sizeof(float), i, fp) != i)
		Cleanup("can't write time values");
	fclose(fp);

	// write the coordinates out
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing coordinates");
		fflush(stdout);
	}
	sprintf(fpath, "%s/%s,%s,VS_Coords", ImgPath, Prefix, ParmName);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open time file for write");
	i = NiiHdr.dim[1];
	j = NiiHdr.dim[2];
	k = NiiHdr.dim[3];
	for(k=c=0; k<NiiHdr.dim[3]; k++)
		for(j=0; j<NiiHdr.dim[2]; j++)
			for(i=0; i<NiiHdr.dim[1]; i++, c++) {
				gsl_matrix_get_row(W, Wgt, c);
				xp = NiiHdr.srow_x[0] * i + NiiHdr.srow_x[1] * j + NiiHdr.srow_x[2] * k + NiiHdr.srow_x[3];
				yp = NiiHdr.srow_y[0] * i + NiiHdr.srow_y[1] * j + NiiHdr.srow_y[2] * k + NiiHdr.srow_y[3];
				zp = NiiHdr.srow_z[0] * i + NiiHdr.srow_z[1] * j + NiiHdr.srow_z[2] * k + NiiHdr.srow_z[3];
				for(m=0, tmp=0.; m<M; m++)
					tmp += gsl_vector_get(W, m) * gsl_vector_get(W, m);
				if(tmp > 0.)
					fprintf(fp, "%6.3f\t%6.3f\t%6.3f\n", 0.1*yp, 0.1*xp, 0.1*zp);
			}
	fclose(fp);

	// that's all there is to it!
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("SAMvs done\n");
		fflush(stdout);
	}
}
