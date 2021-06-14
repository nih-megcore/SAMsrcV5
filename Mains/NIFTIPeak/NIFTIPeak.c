// NIFTIPeak -- find the peaks in all NIFTI volumetric image files
//	in the SAM subdirectory of the named dataset. This version
//  finsa local maxima using a 3x3x3 kernel. If no image directory
//  is specified, then the NIFTI files will be located in
//  <datasetname>/SAM/Image. In this case, the .max files are
//  written to <datasetname/SAM/Max. Note that, SAMwts & SAMNwts
//  searches the Max subdirectory _if_ no ROI and no targets are
//  specified. NIFTI target files are identified by the dataset
//  name & its parameter file name.
//
// This version write a .tag file into the subject's MRI directory.
//  Opening this in AFNI enables selection of local maxima V1 through
//  Vn.
//
//	flags:
//	-r <run name>	-- MEG run name
//	-t <threshold>	-- threshold value (usually z-deviate)
//	-n <max>		-- maximum number of maxima
//	-v				-- verbose mode
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <model.h>
#include <samlib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <nifti1.h>
#include <version.h>

#define MINOR_REV	2

typedef struct {
    float   Vc[3];
    float   val;
    int     rank;
} MAXIMA;


int	main(
	int		argc,
	char	**argv
)

{
    HeaderInfo      Header;         // MEG data header
    ChannelInfo     *Channel;       // MEG channel info
    EpochInfo       *Epoch;         // MEG epoch info
    PARMINFO        Params;         // analysis parameters
    MAXIMA          *Maximum;       // Maximum[N] -- maxima coordinates & value
	nifti_1_header	NiiHdr;			// .nii header
	float			**Image;		// Image[D][V] -- image voxels
	float           kernel[27];     // 3x3x3 kernel
	float			Vs[3];			// starting coordinates
	float			Vstep[3];		// step size for each axis
	float           cnt;            // count of non-zero voxels
	float			rms;			// rms non-zero voxel values
	float			sum;			// sum of non-zero voxel values^2
	float			max;			// maximum voxel value
	float			voxel;			// voxel value;
	float			ImgThresh;		// threshold for current image
    static float    threshold = 0.; // image threshold
    int             V;              // number of maxima voxels
    int             c;              // input character & virtual channel
    int             cc;             // secondary voxel index
    int             eflg;           // command-line error flag
	int				flgcnt;			// command line flag count
	int             d;              // kernel d-index
	int             n;              // maxima index
	int             index;          // kernel maximum index
	int				i;              // kernel i-index
	int				j;              // kernel j-index
	int				k;              // kernel k-index
	int             x;              // image x-index
	int             xx;             // alternate x-index
	int             y;              // image y-index
	int             yy;             // alternate y-index
	int             z;              // image z-index
	int             zz;             // alternate z-index
	int				D;				// number of 3d images
	int				I;				// number of i-axis voxels
	int				J;				// number of j-axis voxels
	int				K;				// number of k-axis voxels
	int				N;				// number of maxima
	int				LongIndex;		// position of arguments
    static int      M = 50;			// default number of maxima
    static int      vflg = FALSE;	// verbose mode flag
    extern char     *optarg;
    extern int      opterr;
    char            fpath[256];     // general path name
    char            DSName[256];    // MEG run name
    char            DSpath[256];    // dataset path
    char            SAMDir[256];    // SAM subdirectory pathname
    char            NiiDir[256];    // directory for finding NIFTI files
    char            ParmName[256];  // parameter file name
    char            NiiName[256];   // SAM volume image name
    char            MaxDir[256];    // path for writing .max files
    char            Prefix[64];     // subject prefix
    char            ShortOpts[] = "r:m:t:s:n:pv";
    static struct option	LongOpts[] = {	// command line options
        {"dataset_name", required_argument, NULL, 'r'},
        {"parameter_name", required_argument, NULL, 'm'},
        {"image_threshold", required_argument, NULL, 't'},
        {"number_of_maxima", required_argument, NULL, 'n'},
        {"verbose", no_argument, NULL, 'v'},
        {NULL, no_argument, NULL, 0}
    };
    unsigned char   **Bad;          // Bad[E][T] -- flags for bad samples
    FILE            *fp;            // maximum list file

    // parse command line parameters
    opterr = 1;		// enable error reporting
    eflg = flgcnt = 0;
    while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
        switch(c) {
            case 'r':   // dataset name
                sprintf(DSName, "%s", optarg);
                flgcnt++;
                break;
            case 'm':   // parameter file name
                sprintf(ParmName, "%s", optarg);
                flgcnt++;
                break;
            case 't':   // image threshold
                threshold = atof(optarg);
                break;
            case 'n':   // limit maxima
                M = atoi(optarg);
                break;
            case 'v':   // verbose mode
                vflg = TRUE;
                break;
            default:
                eflg++;
                break;
        }
    }
    if(eflg > 0 || flgcnt != 2) {
        fprintf(stderr, "NIFTIPeak\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
        fprintf(stderr, "\t-r <dataset name>\n");
        fprintf(stderr, "\t-m <parameter file name>\n");
        fprintf(stderr, "\t-t <image threshold value (negative value over-rides rms threshold)>\n");
        fprintf(stderr, "\t-n <limit number of maxima (default 50)>\n");
        fprintf(stderr, "\t-v -- verbose mode\n");
        exit(-1);
    }

    // announce
    if(vflg == TRUE) {
        printf("N I F T I   V O L U M E   F I L E   I M A G E   P E A K   F I N D E R\t%s\n", PRG_REV);
        fflush(stdout);
    }

    // get data with sensor structures
    GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);
    sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
    sprintf(SAMDir, "%s/SAM", DSpath);

	// parse analysis parameter specification file in current working directory
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("parsing '%s' parameter file", ParmName);
		fflush(stdout);
	}
	sprintf(fpath, "%s.param", ParmName);
	GetParams(fpath, &Params);

    // set prefix
    memset(Prefix, 0, 64);
	memcpy(Prefix, Header.SetName, Params.NumPrefix);

	// NIFTI files are found either in the dataset Image subdirectory or in the image directory
	if(strncmp(Params.DirName, "NULL", 4)) {
		sprintf(NiiDir, "%s", Params.DirName);
	} else {
		sprintf(NiiDir, "%s/Image", SAMDir);
	}

	// maxima files are written to either the dataset Max subdirectory or to the MRI subdirectory
	if(strncmp(Params.MRIdirectory, "NULL", 4) && strncmp(Params.DirName, "NULL", 4)) {
		sprintf(MaxDir, "%s/%s", Params.MRIdirectory, Prefix);
	} else {
        sprintf(MaxDir, "%s/Max", SAMDir);
        if(mkdir(MaxDir, S_IRWXU | S_IRWXG | S_IRWXO) == -1)
            if(errno != EEXIST)
                Cleanup("can't create 'Max' subdirectory");
	}

    // allocate arrays for maxima
    if((Maximum = (MAXIMA *)malloc(1000 * sizeof(MAXIMA))) == NULL)
        allocfailed("Maximum[]");

    // loop through directory entries
    sprintf(NiiName, "%s,%s.nii", Header.SetName, ParmName);
    sprintf(fpath, "%s/%s", NiiDir, NiiName);
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("reading NIFTI image %s", fpath);
        fflush(stdout);
    }
    GetNIFTI(fpath, &D, &Image, &NiiHdr);
    if(D != 1)
        Cleanup("3d+time images not supported");

    // set up arrays for voxel coordinates & values (convert from mm)
    Vs[X_] = 0.001 * NiiHdr.srow_y[3];
    Vs[Y_] = -0.001 * NiiHdr.srow_x[3];
    Vs[Z_] = 0.001 * NiiHdr.srow_z[3];
    Vstep[X_] = 0.001 * NiiHdr.srow_y[2];   // i-axis SAM step size
    Vstep[Y_] = -0.001 * NiiHdr.srow_x[1];  // j-axis "
    Vstep[Z_] = 0.001 * NiiHdr.srow_z[0];   // k-axis "
    I = NiiHdr.dim[3];                      // slowest changing index (x)
    J = NiiHdr.dim[2];                      // next (y)
    K = NiiHdr.dim[1];                      // most rapidily changing index (z)
    V = I * J * K;

    // determine rms non-zero voxel absolute value
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("finding rms value of all non-zero voxels");
        fflush(stdout);
    }
    for(c=0, sum=cnt=0.; c<V; c++) {
        voxel = fabsf(Image[0][c]);
        if(voxel > 0.) {
            sum += voxel * voxel;
            cnt++;
        }
    }
    if(cnt > 1.) {
        rms = sqrtf(sum / cnt);
        if(threshold == 0.)
            ImgThresh = rms;
        else
            ImgThresh = threshold;
    }

    // search for maxima
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("rms = %7.3f rms\n", rms);
        printf("searching for maxima using threshold = %7.3f", ImgThresh);
        fflush(stdout);
    }

    // scan kernel through image voxels
    for(x=n=0; x<I-3; x++)
        for(y=0; y<J-3; y++)
            for(z=0; z<K-3; z++) {

                // copy to 3x3x3 kernel
                for(i=0, xx=x; i<3; i++, xx++)
                    for(j=0, yy=y; j<3; j++, yy++)
                        for(k=0, zz=z; k<3; k++, zz++) {
                            d = i * 9 + j * 3 + k;  // kernel index
                            cc = xx * J * K + yy * K + zz;
                            kernel[d] = Image[0][cc];
                        }

                // find maximum within 3x3x3 kernel
                for(d=0, max=0., index=-1; d<27; d++)
                if(kernel[d] > max) {
                    max = kernel[d];
                    index = d;
                }

                // is this a local maximum?
                if(index == 13 && kernel[13] >= ImgThresh) {    // middle kernel value is greatest! -- mark as a maximum
                    Maximum[n].Vc[X_] = Vs[X_] + (x + 1) * Vstep[X_];
                    Maximum[n].Vc[Y_] = Vs[Y_] + (y + 1) * Vstep[Y_];
                    Maximum[n].Vc[Z_] = Vs[Z_] + (z + 1) * Vstep[Z_];
                    Maximum[n].val = kernel[index];
                    Maximum[n].rank = -1;
                    n++;
                }
            } // end x, y, z
    N = n;
    if(M > N)
        M = N;

    // were any maxima found above threshold?
    if(n == 0) {
        if(vflg == TRUE) {
            printf(" - done\n");
            printf("no maxima found above threshold\n");
            fflush(stdout);
        }
    } else {
        if(vflg == TRUE) {
            printf("\tfound %d local maxima\n", N);
            fflush(stdout);
        }

        // sort maxima in descending order
        for(i=0; i<N; i++) {
            for(n=0, max=0.; n<N; n++)
                if(Maximum[n].rank == -1 && Maximum[n].val > max) {
                    max = Maximum[n].val;
                    j = n;
                }
            Maximum[j].rank = i;
        }

        // open file to list maxima
        sprintf(fpath, "%s/%s.max", MaxDir, NiiName);
        if((fp = fopen(fpath, "w")) == NULL)
            Cleanup("'MaxList' open for write failed");

        // write maxima
        if(vflg == TRUE) {
            printf("saving %d out of %d maxima", M, N);
            fflush(stdout);
        }
        for(i=0; i<M; i++)
            for(n=0; n<N; n++)
                if(Maximum[n].rank == i) {
                    fprintf(fp, "%7.3f\t%7.3f\t%7.3f\t%7.3f\n",
                        100.*Maximum[n].Vc[X_], 100.*Maximum[n].Vc[Y_], 100.*Maximum[n].Vc[Z_], Maximum[n].val);
                    fflush(fp);
                }
        fclose(fp);

        // write tag file
        if(vflg == TRUE) {
            printf(" - done\n");
            printf("writing AFNI tag file");
            fflush(stdout);
        }
        sprintf(fpath, "%s/%s/%s.tag", Params.MRIdirectory, Prefix, DSName);
        if((fp = fopen(fpath, "w")) == NULL)
            Cleanup("tagfile open for write failed");
        for(i=0; i<M; i++)
            for(n=0; n<N; n++)
                if(Maximum[n].rank == i) {
                    fprintf(fp, "'V%02d' %d %d %d\n", i+1,
                        (int)rint(1000.*Maximum[n].Vc[Y_]), (int)rint(-1000.*Maximum[n].Vc[X_]), (int)rint(1000.*Maximum[n].Vc[Z_]));
                    fflush(fp);
                }
        fclose(fp);
    }

    // free allocations
    free((void *)Maximum);

    if(vflg == TRUE)
        printf("\n'NIFTIPeak' done\n");
    exit(0);
}
