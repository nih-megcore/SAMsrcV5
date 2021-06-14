// MUTimg -- synthetic aperture magnetometry imaging for source event related
//	time-locked signals. MUTimg computes an array of symbols for each voxel.
//	It then computes all combinations of mutual information between voxels &
//	writes out a 3D image file.
//
//	flags:
//		-r <dataset name>		-- MEG dataset name
//		-m <trigger name>		-- name of window parameter file
//		-v						-- verbose mode
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

// include...
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <fcntl.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <time.h>
#include <sys/time.h>
#include <Params.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <samUtils.h>
#include <rvelib.h>
#include <filters.h>
#include <samlib.h>
#include <nifti1.h>
#include <version.h>

// internal constants
#define MINOR_REV	0
#define THRESH		3.0					// surrogate data threshold
#define ORIG		0
#define TLRC		1
#define TRUE		1
#define FALSE		0
#define SIG_ORDER	4


int	main(
	int		argc,
	char	**argv)
{
	IIRSPEC			DataIIR;			// IIR filter structure
	FFTSPEC			DataFFT;			// FFT data filter structure
	PARMINFO		Params;				// analysis parameters
	HeaderInfo		Header;				// MEG data header
	ChannelInfo		*Channel;			// MEG channel info
	EpochInfo		*Epoch;				// MEG epoch info
	nifti_1_header	NiiHdr;				// NIFTI header
	unsigned char	*ExtHdr;			// NiiHdr[N] -- extended NIFTI header
	gsl_matrix		*Wgt;				// Wgt - VxM ROI SAM weights
	gsl_vector		*W;					// W - Mx1 SAM weight vector
	double			**x;				// x[TT][M] -- filtered data
	double			*sam;				// sam[TT] -- virtual sensor time-series
	double			*In;				// In[TT] -- input time-series
	double			*Out;				// Out[TT] -- output time-series
	double			BWFreq;				// effective bandwidth
	double			percent;			// percent done
	double			scale;				// bits-->Tesla conversion
	double			max;				// maximum distance & maximum seed image value
	double 			tmp;				// what else?
	float			**Info;				// Info[V][V] -- mutual information matrix
	float			Surrogate;			// surrogate value
	float			voxel;				// voxel value
	short			**SymData;			// SymData[V][TT] -- symbols for data at each voxel & each time step
	short			**SymSurr;			// SymSurr[V][TT] -- symbols for surrogate data
	int				*VoxROI;			// VoxROI[V] -- flags for ROI
	int				*VoxIndex;			// VoxIndex[N] -- index of true ROI voxel indices
	int				E;					// number of good epochs
	int				M;					// number of MEG or reference channels used
	int				N;					// number of bytes in extended header
	int				S;					// number of symbols
	int				T;					// number of samples in epoch
	int				TT;					// total number of samples
	int				TL;					// length of symbol array
	int				V;					// number of voxels in image
	int				NumVoxels;			// number of true ROI voxels
	int				NumMask;			// number of voxels within mask
	int				lags;				// number of sample lags
	int				stride;				// embedding space lags
	int				c;					// input character & voxel index
	int				cc;					// alternate voxel index
	int				cx;					// voxel index
	int				cy;					// 2nd voxel index
	int				pct;				// integer percentage
	int				i;					// general index
	int				j;					// another index
	int				e;					// epoch index
	int				f;					// 1st primary channel index
	int				m;					// channel index
	int				n;					// alternate channel index
	int				s;					// seed index
	int				t;					// sample index
	int				tt;					// alternate sample index
	int				v;					// voxel index
	int 			nvx;				// number of voxels in mask
	int				LongIndex;			// position of arguments
	static int		bflg = FALSE;		// read binary flag
	static int		eflg = FALSE;		// command-line error flag
	static int		kflg = FALSE;		// mask flag
	static int		mflg = FALSE;		// marker parameter file name flag
	static int		rflg = FALSE;		// run name flag
	static int		vflg = FALSE;		// verbose mode flag
	extern char		*optarg;
	extern int		opterr;
	char			fpath[256];			// general path name
	char			DSName[256];		// MEG dataset name
	char			DSpath[256];		// MEG dataset path
	char			SAMpath[256];		// SAM subdirectory path
	char			ImgPath[256];		// save image path (optional)
	char			ParmName[256];		// name of marker parameter file
	char			Prefix[64];			// incomplete image name
	char			extension[4];		// extension
	unsigned char	**Bad;				// Bad[E][T] -- flags for bad samples
	unsigned char	*Mask;				// mask data
	FILE			*fp;				// data output file pointer
	char			ShortOpts[] = "r:m:bv";
	static struct option	LongOpts[] = {	// command line options
		{"dataset", required_argument, NULL, 'r'},
		{"analysis_name", required_argument, NULL, 'm'},
		{"read_binary", no_argument, NULL, 'b'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};

	// parse command line parameters
	opterr = 1;				// enable error reporting
	while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
		switch (c) {
			case 'r':		// datset name
				sprintf(DSName, "%s", optarg);
				rflg = TRUE;
				break;
			case 'm':		// time-segment parameter file name
				sprintf(ParmName, "%s", optarg);
				mflg = TRUE;
				break;
			case 'b':		// read binary instead of recomputing
				bflg = TRUE;
				break;
			case 'v':		// verbose mode
				vflg = TRUE;
				break;
			default:
				eflg = TRUE;
				break;
		}
	}
	if(eflg == TRUE || rflg == FALSE || mflg == FALSE) {
		fprintf(stderr, "MUTimg [arguments]\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "\t-r <dataset name>\n");
		fprintf(stderr, "\t-m <parameter file name>\n");
		fprintf(stderr, "\t-b -- read binary file\n");
		fprintf(stderr, "\t-v -- verbose mode\n");
		exit(-1);
	}

	// get dataset information structures
	if(vflg == TRUE) {
		printf("S A M   M U T U A L   I N F O R M A T I O N   I M A G I N G\t%s\n", PRG_REV);
		printf("opening data file");
		fflush(stdout);
	}

	// get data with sensor structures
	GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);
	sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
	sprintf(SAMpath, "%s/SAM", DSpath);

	// set constants
	M = Header.NumPri;
	E = Header.NumEpochs;
	T = Header.MaxSamples;
	TT = E * T;
	f = Header.PsIndex[0];		// 1st MEG primary sensor channel index
	if(E != 1)
		fprintf(stderr, "reminder: multi-epoch data must be contiguous\n");

	// parse analysis parameter specification file in current working directory
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("parsing '%s' parameter file:\n", ParmName);
		fflush(stdout);
	}
	sprintf(fpath, "%s.param", ParmName);
	GetParams(fpath, &Params);
	if(Params.ImageMetric != MUT_INFO)
		Cleanup("ImageMetric 'MutualInfo' req'd");
	if(Params.NumMark != 0)
		Cleanup("'NumMarkers' specification must be 0");
	if(Params.CovHP == -999. || Params.CovLP == -999.)
		Cleanup("CovBand specifications required");
	if(Params.ImageHP == -999. || Params.ImageLP == -999.)
		Cleanup("ImageBand specifications required");
	if(Params.ImageHP < Channel[f].AnalogHpFreq)
		Cleanup("ImageBand highpass cutoff frequency out of recorded range");
	if(Params.ImageLP > Channel[f].AnalogLpFreq)
		Cleanup("ImageBand lowpass cutoff frequency out of recorded range");
	if(Params.CovType == -1)
		Cleanup("'CovType' is mandatory");
	if(strncmp(Params.DirName, "NULL", 4))
		sprintf(ImgPath, "%s", Params.DirName);
	else
		sprintf(ImgPath, "%s", SAMpath);
	if(Params.Dims < 4 || Params.Dims > 5)
		Cleanup("embedding space dimension must be 4 or 5");
	switch(Params.Dims) {
		case 3: S = 6;
		case 4: S = 24;
		case 5: S = 120;
		case 6: S = 720;
		default: cleanup("embedding space dimension out of range"); break;
	}

	// set prefix to NumPrefix characters (default is 8-character dataset hashcode)
	memset(Prefix, 0, 64);
	memcpy(Prefix, DSName, Params.NumPrefix);

	// weight name is merely the name of the trigger -- it must be extended using the params...
	if(vflg == TRUE) {
		switch(Params.ImageFormat) {
			case ORIG: printf("reading Global weights"); break;
			case TLRC: printf("reading Global weights (Talairach)"); break;
			default: break;
		}
		fflush(stdout);
	}
	switch(Params.ImageFormat) {
		case ORIG: sprintf(fpath, "%s/%s,%-d-%-dHz/Global.nii", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP); break;
		case TLRC: sprintf(fpath, "%s/%s,%-d-%-dHz/Global_at.nii", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP); break;
		default: break;
	}
	Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
	V = Wgt->size1;
	if(Wgt->size2 != M)
		Cleanup("number of channels in dataset does not match number in weight file");
	if(Params.Lags < 0. || Params.Lags > 20.)
		Cleanup("sample lags must be between 0 & 20 ms");
	lags = (int)rint(0.001 * Params.Lags * Header.SampleRate);
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("computing symbolic mutual information between pairs of voxels\n");
		printf("\tfor advance of %d samples (%4.1f ms)\n", lags, Params.Lags);
		switch(Params.Dims) {
			case 4: printf("\t4-dimensional embedding space"); break;
			case 5: printf("\t5-dimensional embedding space"); break;
			default: break;
		}
		fflush(stdout);
	}

	// if req'd, read mask
	if(Params.ImageFormat == TLRC && strncmp(Params.Mask, "NULL", 4)) {
		if(vflg == TRUE) {
			printf(" - done\n");
			printf("reading mask");
			fflush(stdout);
		}
		GetNIFTIMask(Params.Mask, &Mask, &nvx);
		if(nvx != V)
			Cleanup("total number of voxels in mask file does not agree with number in weight file");
		for(v=0, NumMask=0; v<nvx; v++)
			if(Mask[v] == 1)
				NumMask++;
		kflg = TRUE;
	} else {
		if((Mask = (unsigned char *)malloc((size_t)V * sizeof(unsigned char))) == NULL)
			allocfailed("Mask[]");
		for(v=0; v<V; v++)
			Mask[v] = 1;
		NumMask = V;
		kflg = FALSE;
	}


	// important! Mask the beamformer coefficients
	W = gsl_vector_alloc(M);
	gsl_vector_set_zero(W);
	for(c=0; c<V; c++)
		if(Mask[c] == 0)
			gsl_matrix_set_row(Wgt, c, W);

	// if req'd, summarize epoch dimensions
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("dataset '%s' duration %6.3f seconds\n", DSName, (double)TT/Header.SampleRate);
		fflush(stdout);
	}

	// count non-zero SAM coefficients
	if(vflg == TRUE) {
		printf("counting non-zero voxels");
		fflush(stdout);
	}
	if((VoxROI = (int *)malloc((size_t)V * sizeof(int))) == NULL)	// allocate space for _all_ voxels
		allocfailed("VoxROI[]");
	for(c=NumVoxels=0; c<V; c++) {
		gsl_matrix_get_row(W, Wgt, c);
		if(gsl_vector_isnull(W) == 0) {
			VoxROI[c] = TRUE;
			NumVoxels++;				// this is the count of ROI voxels where coefficients are non-zero
		} else {
			VoxROI[c] = FALSE;
		}
	}
	if((VoxIndex = (int *)malloc((size_t)NumVoxels * sizeof(int))) == NULL)
		allocfailed("VoxIndex[]");
	for(c=cc=0; c<V; c++)				// set an index of non-zero voxels
		if(VoxROI[c] == TRUE) {
			VoxIndex[cc] = c;
			cc++;
		}

	// allocate memory
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("number of voxels inside hull boundary: %d out of %d", NumVoxels, V);
		if(kflg == TRUE)
			printf(" (%d inside mask)", NumMask);
		printf("\nallocating memory");
		fflush(stdout);
	}
	if((SymData = (short **)malloc((size_t)NumVoxels * sizeof(short *))) == NULL)
		allocfailed("SymData[]");
	for(c=0; c<NumVoxels; c++)
		if((SymData[c] = (short *)malloc((size_t)TT * sizeof(short))) == NULL)
			allocfailed("SymData[][]");
	if((SymSurr = (short **)malloc((size_t)NumVoxels * sizeof(short *))) == NULL)
		allocfailed("SymSurr[]");
	for(c=0; c<NumVoxels; c++)
		if((SymSurr[c] = (short *)malloc((size_t)TT * sizeof(short))) == NULL)
			allocfailed("SymSurr[][]");
	if((x = (double **)malloc((size_t)TT * sizeof(double *))) == NULL)
		allocfailed("x[]");
	for(t=0; t<TT; t++)
		if((x[t] = (double *)malloc((size_t)M * sizeof(double))) == NULL)
			allocfailed("x[][]");
	if((Info = (float **)malloc((size_t)NumVoxels * sizeof(float *))) == NULL)
		allocfailed("Info[]");
	for(c=0; c<NumVoxels; c++)
		if((Info[c] = (float *)malloc((size_t)NumVoxels * sizeof(float))) == NULL)
			allocfailed("Info[][]");
	if((sam = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
		allocfailed("sam[]");
	if((In = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
		allocfailed("In[]");
	if((Out = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
		allocfailed("Out[]");

	// make filter
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("making %6.3f to %6.3f Hz filter", Params.ImageHP, Params.ImageLP);
		fflush(stdout);
	}
	if(Params.FilterType == IIR)
		mkfilter(&DataIIR, SIG_ORDER, Header.SampleRate, Params.ImageHP, Params.ImageLP, &BWFreq);
	else
		Butterworth(&DataFFT, TT, Params.ImageHP, Params.ImageLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

	// read & process data
	if(bflg == FALSE) {

		// read & filter MEG data
		if(vflg == TRUE) {
			printf(" - done\n");
			printf("reading & filtering MEG data");
			fflush(stdout);
		}
		for(m=0; m<M; m++) {							// read & process data in channel order
			n = Header.PsIndex[m];
			scale = Channel[n].Scale;
			for(e=0; e<E; e++) {
				GetDsData(&Header, e, n, scale, In);	// get one epoch of data for this channel -- length T samples
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
		}

		// loop on all coordinates in the ROI
		if(vflg == TRUE) {
			printf(" - done\n");
			printf("computing symbols for ROI voxels:\n");
			fflush(stdout);
		}
		for(cx=0, pct=(-1); cx<NumVoxels; cx++) {

			// get actual voxel index for cx
			c = VoxIndex[cx];
			gsl_matrix_get_row(W, Wgt, c);

			// compute SAM time-series for voxel
			for(t=0; t<TT; t++) {
				for(m=0, tmp=0.; m<M; m++)
					tmp += x[t][m] * gsl_vector_get(W, m);
				sam[t] = tmp;
			}

			// compute array of symbols
			Data2SymbolArray(sam, SymData[cx], TT, &TL, Params.Dims, Params.ImageLP, Header.SampleRate);
			RandomSymbol(SymSurr[cx], TT, S);

			// percent done
			if(vflg == TRUE) {
				percent = 100. * (double)cx / (double)NumVoxels;
				if(rint(percent) > (double)pct) {
					pct = (int)rint(percent);
					printf("\r\t%d percent done", pct);
					fflush(stdout);
				}
			}
		}

		stride = Header.SampleRate / (2. * Params.ImageLP);		// the lags depend upon the sample rate & upper frequency limit
		i = Params.Dims * stride;
		TL = TT - (i - 1) * lags - 1;

		// free up the memory allocated for symbols
		if(vflg == TRUE) {
			printf(" - done\n");
			printf("freeing unused memory allocations");
			fflush(stdout);
		}
		for(t=0; t<TT; t++)
			free((void *)x[t]);
		free((void *)x);
		free((void *)sam);
		free((void *)In);
		free((void *)Out);

		// open file for seed maxima
		sprintf(fpath, "%s/%s,%-d-%-dHz,maxima.txt", SAMpath, ParmName, (int)Params.ImageHP, (int)Params.ImageLP);
			if((fp = fopen(fpath, "w")) == NULL)
				Cleanup("can't open seed maxim list file for write");

		// loop on all coordinates in the ROI & compute symbolic mutual information
		if(vflg == TRUE) {
			printf(" - done\n");
			printf("computing mutual information matrix:\n");
			fflush(stdout);
		}
		for(cx=0, pct=(-1); cx<NumVoxels; cx++) {		// 'cx' is the seed voxel -- 'cy' is advanced relative to 'cx'
			for(cy=0, max=0.; cy<NumVoxels; cy++) {
				Info[cx][cy] = Sym2MI(SymData[cx], SymData[cy], TL, lags, Params.Dims);
				Surrogate = Sym2MI(SymData[cx], SymSurr[cy], TL, lags, Params.Dims);
				if(THRESH*Surrogate > Info[cx][cy])
					Info[cx][cy] = 0.;
				if(Info[cx][cy] > max)
					max = Info[cx][cy];
			}

			// write seed maximum
			fprintf(fp, "seed %d: max=%f\n", cx, max);

			// percent done
			if(vflg == TRUE) {
				percent = 100. * (double)cx / (double)NumVoxels;
				if(rint(percent) > (double)pct) {
					pct = (int)rint(percent);
					printf("\r\t%d percent done", pct);
					fflush(stdout);
				}
			}
		}
		fclose(fp);

		// write mutual information matrix
		printf("\nwriting mutual information matrix");
		sprintf(fpath, "%s/%s,%-d-%-dHz,Lag%-d,Connections.bin", SAMpath, ParmName, (int)Params.ImageHP, (int)Params.ImageLP, lags);
		if((fp = fopen(fpath, "w")) == NULL)
			Cleanup("can't open mutual information matrix file for write");
		for(i=0; i<NumVoxels; i++)
			for(j=0; j<NumVoxels; j++)
				if(fwrite((void *)&Info[i][j], sizeof(float), 1, fp) != 1)
					Cleanup("can't write data image file");
		fclose(fp);

		// free up the memory allocated for symbols
		if(vflg == TRUE) {
			printf(" - done\n");
			printf("freeing unused memory allocations");
			fflush(stdout);
		}
		for(c=0; c<NumVoxels; c++) {
			free((void *)SymData[c]);
			free((void *)SymSurr[c]);
		}
		free((void *)SymData);
		free((void *)SymSurr);

	} else {		// read data from binary file

		// read mutual information matrix
		printf("\nreading mutual information matrix");
		sprintf(fpath, "%s/%s,%-d-%-dHz,Lag%-d,Connections.bin", SAMpath, ParmName, (int)Params.ImageHP, (int)Params.ImageLP, lags);
		if((fp = fopen(fpath, "r")) == NULL)
			Cleanup("can't open mutual information matrix file for read");
		for(i=0; i<NumVoxels; i++)
			for(j=0; j<NumVoxels; j++)
				if(fread((void *)&Info[i][j], sizeof(float), 1, fp) != 1)
					Cleanup("can't read data image file");
		fclose(fp);
	}

	// write out all images
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("opening image file");
		fflush(stdout);
	}

	// generate output file name & open image file
	sprintf(fpath, "%s/%s,%s,MUT.nii", ImgPath, Prefix, ParmName);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open .nii file for write");

	// modify NIFTI header
	memset(extension, 0, 4);								// for now, ignore the fact that we may have an extended header
	NiiHdr.dim[4] = (kflg == TRUE)? NumMask: NumVoxels;		// number of 3d images
	NiiHdr.toffset = 0.;									// time offset shift
	NiiHdr.slice_duration = 1.;
	NiiHdr.slice_start = 0;
	NiiHdr.slice_end = (kflg == TRUE)? NumMask - 1: NumVoxels - 1;
	NiiHdr.slice_code = NIFTI_SLICE_SEQ_INC;
	memset(NiiHdr.intent_name, 0, 16);
	memcpy(NiiHdr.intent_name, "Mutual Info", 11);

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

	// write out seed images
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing seed images");
		fflush(stdout);
	}

	// write out seed images
	for(v=s=0; v<V; v++)					// loop on all voxels -- where V is the number within a 3d image at 8mm
		if(Mask[v] == 1) {					// is this voxel within the mask?
			if(VoxROI[v] == TRUE) {			// yes! now test to see if this has valid coefficients
				for(c=cc=0; c<V; c++) {		// yes! must write all V voxels to 3d image
					if(VoxROI[c] == TRUE) {
						voxel = Info[s][cc];
						cc++;
					} else {
						voxel = 0.;
					}
					if(fwrite((void *)&voxel, sizeof(float), 1, fp) != 1)
						Cleanup("can't write data image file");
				}
				s++;						// increment the seed index into the Info matrix
			} else {						// write an empty 3d seed image if within the mask but without coefficients
				voxel = 0.;
				for(c=0; c<V; c++)
					if(fwrite((void *)&voxel, sizeof(float), 1, fp) != 1)
						Cleanup("can't write data image file");
			}
		}
	fclose(fp);

	// is that all there is?
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("\n'MUTimg' done\n");
		fflush(stdout);
	}
	exit(0);
}
