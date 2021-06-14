// SAMepi -- synthetic aperture magnetometry for imaging epiletic activity.
//  SAMepi reads SAM weights for volume image, & computes the source time-series
//  for each voxel.
//
//  Unlike prior versions, the excess kurtosis is computed for short overlapping
//  segments.
//
//  This algorithm is dedicated to the memory of Max Burbank. Max had
//  encouraged development of the original SAM(g2) algorithm as a means
//  of enabling automated studies at clinical labs everywhere in the
//  world.
//
//  Yuval Harpaz (BIU) suggested counting spikes in overlapping segments.
//  A more general method is to compute excess kurtosis (g2) for overlapping
//  segments of ~1 or 2 seconds - depending on spike frequency. The SAMepi
//  voxel values are the rms value of g2 for all positive non-zero segments.
//
//  flags:
//      -r <run name>   -- MEG run name
//      -d <data name>  -- MEG data file name (4D, only)
//      -m <name>       -- parameter file name
//      -t <threshold>  -- g2 threshold (optional)
//      -v              -- verbose mode
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <filters.h>
#include <nifti1.h>
#include <version.h>
#include "samutil.h"
#include "sam_parse.h"

// #define MINOR_REV   0       // initial version using GSL
// #define MINOR_REV   1       // added flag to use negative g2 with positive
// #define MINOR_REV   2       // changed -n to -t flag to set g2 threshold
// #define MINOR_REV   3       // automatic g2 threshold -- decrease threshold by 0.707 if all zero voxels
// #define MINOR_REV   4       // added stride, based on 2x lowpass
#define MINOR_REV   5       // new command line parser
#define G2_THRESH   10.0    // g2 threshold for window integration

int main(
    int     argc,
    char    **argv
) {
    FFTSPEC         DataFFT;                // FFT data filter structure
    FFTSPEC         NotchFFT;               // optional notch filter
    PARMINFO        Params;                 // analysis parameters
    HeaderInfo      Header;                 // MEG data header
    ChannelInfo     *Channel;               // MEG channel info
    EpochInfo       *Epoch;                 // MEG epoch info
    nifti_1_header  NiiHdr;                 // NIFTI header
    unsigned char   *ExtHdr;                // NiiHdr[N] -- extended NIFTI header
    gsl_matrix      *x;                     // x - MxTT multi-sensor signal
    gsl_matrix      *Wgt;                   // Wgt - VxM SAM coefficients
    gsl_vector      *Wv;                    // Wv - Mx1 SAM coefficient vector
    gsl_vector      *vs;                    // vs - TTx1 virtual sensor time series
    double          *In;                    // In[TT] -- input time-series
    double          *Out;                   // Out[TT] -- output time-series
    double          *Data;                  // Data[T] -- one epoch dataset time-sers
    double          *Seg;                   // Seg[D] -- g2 value for each segment
    double          *Win;                   // Win[TW] - window of data for segment
    double          g2;                     // excess kurtosis
    double          BWFreq;                 // filter bandwidth
    double          scale;                  // bits-->T conversion
    double          percent;                // percent done
    double          g2_thresh = G2_THRESH;  // initialize with default g2 threshold
    float           *Image;                 // Image[V] -- SAM(g2) image
    int             D;                      // number of segments
    int             E;                      // number of good epochs
    int             M;                      // number of MEG or reference channels used
    int             N;                      // number of characters in extended NIFTI header
    int             T;                      // number of samples
    int             TW;                     // number of samples in window
    int             TT;                     // total sample size
    int             V;                      // number of voxels
    int             c;                      // input character & virtual channel
    int             dd;                     // number of segments used for g2
    int             pct;                    // integer percentage
    int             m;                      // channel index
    int             n;                      // alternate channel index
    int             d;                      // segment index
    int             e;                      // epoch index
    int             f;                      // 1st primary channel index
    int             t;                      // sample index
    int             ts;                     // starting sample for segment
    int             tt;                     // alternate time index
    int             to;                     // time offset in epoch
    int             stride;                 // stride for gsl_stats_kurtosis
    int             flag;                   // flag for good/bad window
    int             eflg = FALSE;           // command-line error flag
    int             mflg = FALSE;           // parameter file specification flag
    int             rflg = FALSE;           // dataset name flag
    int             vflg = FALSE;           // verbose mode flag
    int             wflg = FALSE;           // SAM weight flag
    int             zflg = FALSE;           // zero voxel flag
    char            fpath[256];             // general path name
    char            *DSName = NULL;         // dataset name
    char            DSpath[256];            // MEG path name
    char            SAMDir[256];            // SAM subdirectory path
    char            ImgDir[256];            // directory for output of SAMepi NIFTI image
    char            ImgPath[256];           // full path for output NIFTI image
    char            ImgName[256];
    char            OutName[256];           // output file name
    char            WtsName[256];
    char            Prefix[64];             // image file prefix
    char            FilterName[8];          // filter name -- IIR | FIR
    char            extension[4];           // extension
    FILE            *fp;                    // output file pointer (to SAM image)

#if BTI
    static int      dflg = FALSE;           // pdf file error flag
    char            PDFName[256];           // pdf file name
#else
    unsigned char   **Bad;                  // Bad[E][T] -- flags for bad samples
#endif

#if 0
            case 't':       // set g2 threshold (optional)
                g2_thresh = atof(optarg);
#endif

    // Register the parameters used by this program.

    reg_usage(argv[0], MINOR_REV, __DATE__);
    reg_std_parm();

    reg_parm("CovType");
    set_parm_arg_help("CovType", "GLOBAL");
    set_parm_help("CovType", "GLOBAL is always used");

    reg_parm("CovBand");
    reg_parm("ImageBand");
    //reg_parm("FilterType");
    reg_parm("Notch");
    reg_parm("Hz");

    reg_parm("TimeInt");
    reg_parm("ImageDirectory");

    reg_parm("MRIDirectory");
    reg_parm("MRIPattern");
    reg_parm("PrefixLength");

    reg_parm("ImageMetric");
    set_parm_arg_help("ImageMetric", "Kurtosis");
    set_parm_help("ImageMetric", "Kurtosis is always used");

    //reg_parm("g2_thresh");

    reg_parm("WtsName");
    reg_parm("OutName");

    // Parse the arguments, environment variables, and parameter files.

    do_parse_args(argc, argv);

    // Extract specified parameters.

    eflg = get_params(&Params);

    p = get_parm("DataSet");
    if (p->set) {
        rflg = TRUE;
        DSName = copy_string(Params.DataSetName);   // DSName is modified
    }
#if BTI
    p = get_parm("PDFName");
    if (p->set) {
        dflg = TRUE;
        PDFName = (char *)p->ptr;
    }
#endif
    p = get_parm("verbose");
    vflg = p->set;
    p = get_parm("param");
    mflg = p->set;

#if BTI
    if (!rflg || !dflg || !mflg || eflg) {
#else
    if (!rflg || !mflg || eflg) {
#endif
        msg("dataset (-r) and parameter file (-m) are required\n");
        do_help();  // doesn't return
    }

    OutName = Params.ParmName;
    p = get_parm("OutName");
    if (p->set) {
        OutName = (char *)p->ptr;
    }
    WtsName = OutName;
    p = get_parm("WtsName");
    if (p->set) {
        WtsName = (char *)p->ptr;
    }

    Params.CovType = GLOBAL;        // global is the only choice
    Params.ImageMetric = KURTOSIS;  // kurtosis is the only choice

    // announce
    if (vflg) {
        printf("S A M   S P I K E   I M A G I N G\t%s\n", PRG_REV);
        printf("opening data file");
        fflush(stdout);
    }

#if BTI
    // get data with sensor structures
    if (GetMagnesInfo(DSName, PDFName, &Header, &Channel, &Epoch) == -1)
        Cleanup("can't read Magnes information");
    sprintf(DSpath, "%s/%s", Header.DsPath, Header.SetName);
#else
    // get data with sensor structures
    GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);
    sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
#endif
    GetFilePath(DSDIR, SAMDir, sizeof(SAMDir), &Params, "SAM", 0);

    // establish epoch, channel, & sample count
    E = Header.NumEpochs;
    T = Header.MaxSamples;
    TT = E * T;
    M = Header.NumPri;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index

    // check to see if this is a continuous dataset
    if (E != 1 && vflg)
        fprintf(stderr, " (reminder: multi-epoch data must be contiguous)");

    // check if this is the right program for you!
    if (Params.CovHP == -999. || Params.CovLP == -999.)
        Cleanup("CovBand parameters required");
    if (Params.ImageHP == -999. || Params.ImageLP == -999.)
        Cleanup("ImageBand parameters required");
    if (Params.ImageHP < Channel[f].AnalogHpFreq)
        Cleanup("ImageBand highpass cutoff frequency out of recorded range");
    if (Params.ImageLP > Channel[f].AnalogLpFreq)
        Cleanup("ImageBand lowpass cutoff frequency out of recorded range");
    if (Params.TimeInt == -999.)
        Cleanup("TimeInt specification missing from parameter file");
    if (Params.FilterType == IIR)
        strcpy(FilterName, "IIR");
    else
        strcpy(FilterName, "FFT");

    // set prefix to NumPrefix characters (default is 8-character dataset hashcode)
    memset(Prefix, 0, sizeof(Prefix));
    memcpy(Prefix, DSName, Params.NumPrefix);

    // set ImgDir & ImgPath
    if (!strncmp(Params.DirName, "NULL", 4)) {      // if no Image directory name, use MRI directory
        sprintf(ImgPath, "%s/%s,%s", Params.MRIdirectory, Prefix, OutName);
    } else {                                        // else, use the image directory
        sprintf(ImgPath, "%s/%s,%s", Params.DirName, Prefix, OutName);
    }

    // get weight name from parameters & read SAM weights
    if (vflg) {
        printf(" - done\n");
        printf("reading Global SAM weights");
        fflush(stdout);
    }
    sprintf(fpath, "%s/%s,%-d-%-dHz/Global.nii", SAMDir, WtsName, (int)Params.CovHP, (int)Params.CovLP);
    Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
    V = Wgt->size1;
    if (Wgt->size2 != M)
        Cleanup("number of channels in dataset does not match number in weight file");

    // set up imaging constants
    if (vflg) {
        printf(" - done\n");
        printf("determining imaging constants:\n");
        fflush(stdout);
    }
    stride = (int)floor(Header.SampleRate / (2. * Params.ImageLP));
    TW = (int)floor(Params.TimeInt * Header.SampleRate);        // length of each g2 integration segment
    D = 2 * TT / TW - 1;                                        // number of overlapping integration segments
    if (vflg) {
        printf("\tepochs = %d, samples/epoch = %d (total %d)\n", E, T, TT);
        printf("\tintegration window = %d samples\n", TW);
        printf("\tstride = %d samples\n", stride);
        printf("\tnumber of overlapping windows = %d\n", D);
        printf("\tg2 threshold = %4.1f\n", g2_thresh);
        fflush(stdout);
    }

    // allocate memory for arrays
    if (vflg) {
        printf("allocating memory");
        fflush(stdout);
    }
    x = gsl_matrix_alloc(M, TT);
    Wv = gsl_vector_alloc(M);
    vs = gsl_vector_alloc(TT);
    Seg = new_arrayE(double, D, "Seg[]");
    Win = new_arrayE(double, TW, "Win[]");
    In = new_arrayE(double, TT, "In[]");
    Out = new_arrayE(double, TT, "Out[]");
    Data = new_arrayE(double, T, "Data[]");
    Image = new_arrayE(float, V, "Image[]");

    // generate Butterworth FFT filters
    if (vflg) {
        printf(" - done\n");
        printf("making %6.3f to %6.3f Hz FFT filter for image bandpass", Params.ImageHP, Params.ImageLP);
        fflush(stdout);
    }
    Butterworth(&DataFFT, TT, Params.ImageHP, Params.ImageLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);
    if (Params.Notch) {
        if (vflg) {
            printf(" - done\n");
            printf("making FFT notch filter");
            fflush(stdout);
        }
        mknotch(&NotchFFT, Params.ImageHP, Params.ImageLP, Header.SampleRate, TT, Params.Hz);
    }

    // read & filter MEG data
    if (vflg) {
        printf(" - done\n");
        printf("reading & filtering MEG data");
        fflush(stdout);
    }
    for(m=0; m<M; m++) {                                // read & process data in channel order
        n = Header.PsIndex[m];
        scale = Channel[n].Scale;
        for(e=tt=0; e<E; e++) {                         // assemble multi-epoch contiquous data for this channel
            GetDsData(&Header, e, n, scale, Data);      // get one epoch of data for this channel -- length T samples
            for(t=0; t<T; t++, tt++)
                In[tt] = Data[t];                       // 'In' is contiguous with no filters and offsets
        }
        demean(In, TT);
        FFTfilter(In, Out, &DataFFT, TT);               // apply bandpass filter
        if(Params.Notch == FALSE) {
            for(t=0; t<TT; t++)
                gsl_matrix_set(x, m, t, Out[t]);
        } else {
            FFTfilter(Out, In, &NotchFFT, TT);          // optionally, apply notch filters
            for(t=0; t<TT; t++)
                gsl_matrix_set(x, m, t, In[t]);
        }

        // channel heartbeat
        if(vflg == TRUE) {
            printf(".");
            fflush(stdout);
        }
    }
    if(vflg == TRUE) {
        printf(" - done\n");
        fflush(stdout);
    }

    // if all voxels are zero, decrease g2 threshold & repeat...
    zflg = TRUE;    // zflg is true until we find at least one non-zero voxels!
    do {

        // now, step through all coordinates
        if(vflg == TRUE) {
            printf("stepping through voxels to determine g2:\n");
            fflush(stdout);
        }
        for(c=0, pct=(-1); c<V; c++) {

            // copy SAM coefficients for this voxel
            gsl_matrix_get_row(Wv, Wgt, c);

            // test for non-zero coefficients
            if(!gsl_vector_isnull(Wv)) {

                // (1) compute virtual sensor time-series
                gsl_blas_dgemv(CblasTrans, 1., x, Wv, 0., vs);      // vs = TTxM * Mx1 == TTx1

                // (2) for each segment, compute excess kurtosis
                for(d=dd=0; d<D; d++) {
                    ts = d * TW / 2;
                    for(t=0, flag=FALSE; t<TW; t++) {
                        tt = ts + t;
                        Win[t] = gsl_vector_get(vs, tt);
#if BTI
                        flag = TRUE;        // all samples in segment TRUE for BTI/4D -- no bad segment list
#else
                        // convert tt to epoch & time offset to check Bad samples
                        e = tt / T;
                        to = e * T - tt;
                        if(Bad[e][to] == GOOD)
                            flag = TRUE;    // segment is TRUE if _all_ samples in Win are good
                        else
                            flag = FALSE;
#endif
                    }
                    if(flag == TRUE) {
                        g2 = gsl_stats_kurtosis(Win, stride, TW/stride);
                        if(g2 >= g2_thresh) {               // ignore g2 values below threshold (default is 0.)
                            Seg[dd] = g2;
                            dd++;
                        }
                    }
                }
                if(dd > 2) {
                    Image[c] = gsl_stats_sd(Seg, 1, dd);    // rms value of all non-zero segments
                    zflg = FALSE;                           // non-zero voxel -- g2 threshold is okay if any voxel is non-zero
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

        // check whether to decrease threshold & repeat
        if(zflg == TRUE) {
            if(g2_thresh > 0.) {
                g2_thresh *= M_SQRT1_2; // decrease g2 threshold by 0.707...
                if(vflg == TRUE) {
                    printf("\ndecreasing g2 threshold to %6.3f & repeating analysis\n", g2_thresh);
                    fflush(stdout);
                }
            } else {
                zflg = FALSE;
            }
        }

    } while(zflg == TRUE);

    // write out SAM(g2) image
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("writing SAM(Epi) image");
        fflush(stdout);
    }
    sprintf(fpath, "%s/%s,%s.nii", ImgPath, DSName, ParmName);
    if((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open .nii file for write");

    // modify NIFTI header
    memset(extension, 0, 4);                                    // for now, ignore the fact that we may have an extended header
    NiiHdr.dim[0] = 3;                                          // 3d volume
    NiiHdr.dim[4] = 0;                                          // number of 3d images
    NiiHdr.toffset = 0.;                                        // time offset shift
    NiiHdr.slice_duration = 1.;
    NiiHdr.slice_code = NIFTI_SLICE_UNKNOWN;
    NiiHdr.xyzt_units = NIFTI_UNITS_MM;                         // units
    memset(NiiHdr.intent_name, 0, 16);
    memcpy(NiiHdr.intent_name, "Epilepsy Spike", 14);

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
    if(fwrite((void *)Image, sizeof(float), V, fp) != V)
        Cleanup("can't write voxel values to .nii file");
    fclose(fp);

    // create symbolic link from subject's MRI directory to Image subdirectory (this gives NIFTIpeak access to the file)
    sprintf(ImgDir, "%s/Image/", SAMDir);
    if(mkdir(ImgDir, S_IRWXU | S_IRWXG | S_IRWXO) == -1)
        if(errno != EEXIST)
            Cleanup("can't create 'Max' subdirectory");
    sprintf(ImgName, "%s/%s,%s.nii", ImgDir, DSName, ParmName);
    symlink(fpath, ImgName);

    // that's all, folks!
    gsl_matrix_free(x);
    gsl_vector_free(Wv);
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("'SAMepi' done\n");
        fflush(stdout);
    }
    exit(0);
}
