// SAMNwts -- linear synthetic aperture magnetometry coefficient solution.
//	Linear solution for the Cartesian moment vector. Optional 2nd covariance
//	matrix for determining the beamformer orientation vector.
//
//	flags:
//		-r <run name>	-- MEG dataset name
//		-m <parameters>	-- parameter file name
//		-Z				-- noise normalization
//		-v				-- Verbose mode
//
//	Author: Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <strings.h>
#include <fcntl.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <model.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <samUtils.h>
#include <byteswap.h>
#include <version.h>

#define MINOR_REV	6

// global variable
int		ORDER;							// spherical harmonic expansion order


int	main(
	int		argc,
	char	**argv)
{
	PARMINFO		Params;				// analysis parameters
	HEAD_MODEL		Model;				// head model structure for *.hdm file
	HULL			Hull;				// cortical hull model -- vertices & normals
	VOXELINFO		*Target;			// Target[N] -- structure for target parameters
	HeaderInfo		Header;				// MEG data header
	ChannelInfo		*Channel;			// MEG channel info
	EpochInfo		*Epoch;				// MEG epoch info
    struct dirent   *dp;                // directory entry structure
    DIR             *dirp;              // directory pointer
	COV_HDR			WgtCovHdr;			// weight covariance header
	COV_HDR			OriCovHdr;			// orientation covariance header
	N_LEADFIELD     *Lfld;              // Lfld[N] -- N Mx3 gsl lead-field matrices
	gsl_matrix		*Co;				// Co - MxM orientation covariance matrix
	gsl_matrix		*Cf;				// Cf - MxM solution covariance matrix
	gsl_matrix		*Wgt;				// Wgt - MxN n-target weights
	gsl_matrix		*Wt;				// Wt - NxM n-target weights (for NIFTI)
	gsl_matrix		*B;                 // B - MxN n-target forward solution matrix
	gsl_vector		*Wone;				// Wone - Mx1 one weight
	double			HPFreq = 0.;		// highpass frequency
	double			LPFreq = 0.;		// lowpass frequency
	double			Noise;				// effective rank-limited noise
	double			Mu_f;				// weight solution regularization
	double			Mu_o;				// orientation regularization
	double			norm;				// noise normalization
	double			tmp;				// what else?
	int				*ChanIndex;			// ChanIndex[M] -- index of channels used for covariance
	int				LongIndex;			// position of arguments
	int				M;					// number of MEG or reference channels used
	int				N;					// number of targets
	int				R;					// number of physical reference SQUID channels
	int				S;					// number of physical SQUID channels
	int				opt;				// option
	int				v;					// vector index
	int				i;					// matrix index
	int				j;					// matrix index
	int				length;				// string length
	int				found;				// string found flag
	static int		eflg = FALSE;		// command-line error flag
	static int		iflg = FALSE;		// alternate MRI directory path flag
	static int		pflg = FALSE;		// parameter file flag
	static int		rflg = FALSE;		// dataset name flag
	static int		vflg = FALSE;		// verbose mode flag (default: silent)
	static int		zflg = FALSE;		// noise normalization flag
	extern char		*optarg;
	extern int		opterr;
	extern int		CoeffFlag;			// for initializing Nolte solution
	extern int		ORDER;				// spherical harmonic expansion order
	char			fpath[256];			// general path name
	char			DSName[256];		// MEG dataset name
	char			DSpath[256];		// MEG dataset path
	char			SAMpath[256];		// SAM subdirectory path
	char			MRIpath[256];		// alternate MRI directory path
	char            TargetPath[256];    // target path
	char            WgtPath[256];       // weight file path
	char            FileName[64];       // target file name
	char			Prefix[64];			// dataset prefix for finding MRI
	char			CovDirName[256];	// covariance subdirectory name
	char			ParmPath[256];		// parameter file path name
	char			ParmName[256];		// parameter file name

	// system-dependent variables
#if BTI
	static int		dflg;				// pdf file name flag
	char			PDFName[256];		// name of data file
	char			ShortOpts[] = "r:d:m:i:Zv";
	static struct option LongOpts[] = {	// command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"pdf_name", required_argument, NULL, 'd'},
		{"analysis_name", required_argument, NULL, 'm'},
		{"mri_path", required_argument, NULL, 'i'},
		{"normalize", no_argument, NULL, 'Z'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};	
#else
	char			ShortOpts[] = "r:m:i:Zv";
	static struct option LongOpts[] = {	// command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"analysis_name", required_argument, NULL, 'm'},
		{"mri_path", required_argument, NULL, 'i'},
		{"normalize", no_argument, NULL, 'Z'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};	
#endif


	// parse command line parameters
	opterr = 1;					// enable error reporting
    while((opt = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
		switch(opt) {
			case 'r':			// run name
				sprintf(DSName, "%s", optarg);
				rflg = TRUE;
				break;
#if BTI
			case 'd':			// data file name
				sprintf(PDFName, "%s", optarg);
				dflg = TRUE;
				break;
#endif
			case 'm':			// analysis name
				sprintf(ParmPath, "%s", optarg);
				pflg = TRUE;
				break;
			case 'i':			// MRI directory path
				sprintf(MRIpath, "%s", optarg);
				iflg = TRUE;
				break;
			case 'Z':			// noise normalization
				zflg = TRUE;
				break;
			case 'v':			// verbose mode
				vflg = TRUE;
				break;
			default:
				eflg = TRUE;
				break;
		}
	}
#if BTI
    if(eflg == TRUE || rflg == FALSE || dflg == FALSE || pflg == FALSE) {
		fprintf(stderr, "SAMwts\t-r <dataset name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "\t-d <data file name>\n");
#else
    if(eflg == TRUE || rflg == FALSE || pflg == FALSE) {
		fprintf(stderr, "SAMwts\t-r <dataset name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
#endif
		fprintf(stderr, "\t-m <analysis parameter file name\n");
		fprintf(stderr, "\t\t(optional parameters)\n");
		fprintf(stderr, "\t-i <MRI directory path\n");
		fprintf(stderr, "\t-Z -- noise normalization\n");
		fprintf(stderr, "\t-v -- verbose mode\n");
		exit(-1);
	}

    // announce
	if(vflg == TRUE) {
		printf("N - S O U R C E   S A M   C O E F F I C I E N T   G E N E R A T O R\t%s\n", PRG_REV);
		printf("opening data file");
		fflush(stdout);
	}

    // get data with sensor structures
#if BTI
	if(GetMagnesInfo(DSName, PDFName, &Header, &Channel, &Epoch) == -1)
		Cleanup("can't read Magnes information");
	sprintf(DSpath, "%s/%s", Header.DsPath, Header.SetName);
#else
    GetDsInfo(DSName, &Header, &Channel, &Epoch, NULL, TRUE);
	sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
#endif
	sprintf(SAMpath, "%s/SAM", DSpath);

	// set ORDER
	ORDER = 20;				// 6th-order spherical harmonic expansion

	// set channel & epoch count constants
	M = Header.NumPri;
	R = Header.NumRef;
	S = Header.NumPri + Header.NumRef;

	// parse analysis parameter specification file in current working directory
	ParseFileName(ParmPath, ParmName);
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("parsing '%s' parameter file", ParmName);
		fflush(stdout);
	}
	sprintf(fpath, "%s.param", ParmPath);
	GetParams(fpath, &Params);
	ORDER = Params.Order;		// default should be 20th-order spherical harmonic expansion

	// added optional MRI directory path from command line
	if(iflg == TRUE) {
		if(access(MRIpath, F_OK) != 0)
			Cleanup("can't access optional MRI directory");
		strcpy(Params.MRIdirectory, MRIpath);
	}

    // generate prefix
    memset(Prefix, 0, 64);
    memcpy(Prefix, DSName, Params.NumPrefix);

	// define forward solution model
	if(Params.Model == -1)
		Cleanup("'Model' specification missing from parameter file");
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("modeling forward solution using ");
		fflush(stdout);
    }

	switch(Params.Model) {

		// single sphere origin -- note that we must add a small random offset so that the eigensystem solution for moment vector succeeds!
		case SSPHERE:
			if(vflg == TRUE) {
				printf("single local sphere origin: { %6.3f, %6.3f, %6.3f } (cm)", 100. * Params.Sp[X_], 100. * Params.Sp[Y_], 100. * Params.Sp[Z_]);
				fflush(stdout);
			}
			model = SPHERE;
			for(i=0; i<M; i++) {
				j = Header.PsIndex[i];
				for(v=X_; v<=Z_; v++)
					Channel[j].Geom.MEGSensor.LocalSphere[v] = Params.Sp[v] + RandomGauss(0., 1.0e-6);
			}
			for(i=0; i<R; i++) {
				j = Header.RsIndex[i];
				for(v=X_; v<=Z_; v++)
					Channel[j].Geom.MEGSensor.LocalSphere[v] = Params.Sp[v] + RandomGauss(0., 1.0e-6);
			}
			Hull.nv = 1;
			if((Hull.vertex = (FIELD *)malloc((size_t)2 * sizeof(FIELD))) == NULL)
				allocfailed("Hull.vertex[1]");
			for(v=X_; v<=Z_; v++)
				Hull.Vo[v] = Params.Sp[v];
			break;

		// multiple local spheres
		case MSPHERE:

			if(vflg == TRUE) {
				printf("multiple local spheres");
				fflush(stdout);
			}

			// read .hdm file into 'Model'
        	sprintf(fpath, "%s/default.hdm", DSpath);
        	GetHDM(fpath, &Model);
			
			// test if the model matches the number of sensors
			if(Model.NumSensors != S)
				Cleanup("number of sensors in '.hdm' file does not match that of the dataset");
			
			// set multiple sphere origins
			for(i=0; i<S; i++) {
				length = (int)strlen(Model.MultiSphere[i].ChanName) - 1;
				for(j=0, found=FALSE; j<Header.NumChannels; j++)
					if(!strncmp(Model.MultiSphere[i].ChanName, Channel[j].ChannelName, length)) {	
						for(v=X_; v<=Z_; v++)
							Channel[j].Geom.MEGSensor.LocalSphere[v] = Model.MultiSphere[i].Origin[v];
						found = TRUE;
					}
				if(found == FALSE)
					Cleanup("could not find matching channel name in '.HDM' file");
			}
			break;
			
		// Nolte realistic head model
		case NOLTE:

			if(vflg == TRUE) {
				printf("'hull.shape'");
				fflush(stdout);
			}
            sprintf(fpath, "%s/%s/hull.shape", Params.MRIdirectory, Prefix);
            GetHull(fpath, &Hull);      // read hull model from MRI directory
			CoeffFlag = TRUE;
			break;
			
		default:
			break;
	}

	// first try to use CovBand frequencies -- use ImageBand only if CovBand has not been specified!
	if(Params.CovHP >= 0. && Params.CovLP >= 0.) {
		HPFreq = Params.CovHP;
		LPFreq = Params.CovLP;
	} else if(Params.ImageHP >= 0. && Params.ImageLP >= 0.) {
		HPFreq = Params.ImageHP;
		LPFreq = Params.ImageLP;
	} else {
		Cleanup("Neither DataBand nor ImageBand were specified");
	}

	// generate CovDirName
	sprintf(CovDirName, "%s,%-d-%-dHz", ParmName, (int)HPFreq, (int)LPFreq);

	// test for target file & verify path
    if(strncmp("NULL", Params.Target, 4)) {     // target is in MRI directory
        sprintf(fpath, "%s/%s/%s", Params.MRIdirectory, Prefix, Params.Target);
        if(access(fpath, F_OK) != 0)
            Cleanup("can't access target file");
	    sprintf(TargetPath, "%s/%s/%s", Params.MRIdirectory, Prefix, Params.Target);
	} else {                                    // target is in Max subdirectory
        // read directory entries for target
        if(vflg == TRUE) {
            printf(" - done\n");
            printf("searching for target file:");
            fflush(stdout);
        }
        sprintf(fpath, "%s/Max", SAMpath);          // set path to Max subdirectory in dataset
        dirp = opendir(fpath);
        found = FALSE;
        while((dp = readdir(dirp)) != NULL) {
            strcpy(FileName, dp->d_name);
            i = strlen(FileName);
            if(!strcmp(&FileName[i-4], ".max")) {   // Yes, this is a .max target file
                found = TRUE;
                break;
            }
        }
        if(found == FALSE)
            Cleanup("no target file found");
        if(vflg == TRUE) {
            printf("\tfound '%s'\n", FileName);
            fflush(stdout);
        }
        closedir(dirp);
        sprintf(TargetPath, "%s/Max/%s", SAMpath, FileName);
        FileName[i-4] = '\0';                       // strip off .max
        sprintf(WgtPath, "%s/%s,%s.wts", SAMpath, ParmName, FileName);
    }

	if(vflg == TRUE) {
		printf("reading target voxel list");
		fflush(stdout);
	}
	GetTargets(TargetPath, &Target, &N);
	if(vflg == TRUE) {
		printf("found %d targets", N);
		fflush(stdout);
	}
		// copy N moment vectors to voxels
	// allocate memory for arrays
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("allocating memory");
		fflush(stdout);
	}
	if((Lfld = (N_LEADFIELD *)malloc(N * sizeof(N_LEADFIELD))) == NULL)
	    allocfailed("Lfld[]");
	for(i=0; i<N; i++)
    	Lfld[i].L = gsl_matrix_alloc(M, 3);
	Wgt = gsl_matrix_alloc(M, N);
	Cf = gsl_matrix_alloc(M, M);
	Co = gsl_matrix_alloc(M, M);
	B = gsl_matrix_alloc(M, N);
	Wone = gsl_vector_alloc(M);

	// read required covariance files
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("reading covariance file(s)");
		fflush(stdout);
	}
	sprintf(fpath, "%s/%s/Global.cov", SAMpath, CovDirName);
	GetCov(fpath, &WgtCovHdr, &ChanIndex, Cf);
	sprintf(fpath, "%s/%s/Orient.cov", SAMpath, CovDirName);
	GetCov(fpath, &OriCovHdr, &ChanIndex, Co);

	// announce regularization used
	if(vflg == TRUE && Params.Operator != -1) {
	    printf(" - done\n");
		printf("regularizing covariance matrices by ");
		switch(Params.Operator) {
			case ADD_MU:
				printf("adding %f fT/root-Hz to diagonal\n", Params.Mu);
				break;
			case MPY_MU:
				printf("multiplying diagonal by %f\n", Params.Mu);
				break;
			default:
				break;
		}
		fflush(stdout);
	}

	// if req'd, regularize
	switch(Params.Operator) {
		case ADD_MU:
			tmp = 1.0e-30 * Params.Mu * Params.Mu;		// convert fT/root-Hz to T^2/Hz
			Mu_f = WgtCovHdr.BWFreq * tmp;
			Mu_o = OriCovHdr.BWFreq * tmp;
			gsl_matrix_add_diag(Cf, Mu_f);
			gsl_matrix_add_diag(Co, Mu_o);
			Noise += Mu_f;
			break;
		case MPY_MU:
			gsl_matrix_mul_diag(Cf, Params.Mu);
			gsl_matrix_mul_diag(Co, Params.Mu);
			Noise *= Params.Mu;
			break;
		default:
			break;	// do nothing
	}

	// apply simultaneous weight solution
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("solving simultaneous beamformer weights");
		fflush(stdout);
	}
    NLeadField(Target, Lfld, &Header, Channel, &Hull, N, Params.Model);
	NSolveFwd(Target, B, Lfld, Co, N);
	NSolveWts(Wgt, B, Cf);

	// if req'd, normalize weights
	if(zflg == TRUE)
		for(i=0; i<N; i++) {	// get one weight vector at a time
			gsl_matrix_get_col(Wone, Wgt, i);
			for(j=0, tmp=0.; j<M; j++)
				tmp += gsl_vector_get(Wone, j) * Noise * gsl_vector_get(Wone, j);
			norm = 1. / sqrt(tmp);
			gsl_vector_scale(Wone, norm);
			gsl_matrix_set_col(Wgt, i, Wone);
		}

	// create weight file
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing SAM legacy weight file");
		fflush(stdout);
	}
	Wt = gsl_matrix_alloc(N, M);
	gsl_matrix_transpose_memcpy(Wt, Wgt);
	PutWts(WgtPath, Wt, Header, WgtCovHdr);

	// now write NIFTI weight file
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("writing NIFTI weight file");
		fflush(stdout);
	}
	Params.SAMStart[X_] = 0.05;
	Params.SAMEnd[X_] = 0.05;
	Params.SAMStart[Y_] = 0.;
	Params.SAMEnd[Y_] = 0.;
	Params.SAMStart[Z_] = 0.001;
	Params.SAMEnd[Z_] = (double)N * 0.001;
	Params.SAMStep = 0.001;
	sprintf(fpath, "%s/%s,%-d-%-dHz/Global.nii", SAMpath, ParmName, (int)HPFreq, (int)LPFreq);
	PutNIFTIWts(fpath, &Params, Wt ,Noise);

	// finit
	gsl_matrix_free(Wt);
	gsl_matrix_free(Cf);
	gsl_matrix_free(Co);
	gsl_vector_free(Wone);
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("'SAMNwts' done\n");
		fflush(stdout);
	}
	exit(0);
}
