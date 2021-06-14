// SAMpower -- synthetic aperture magnetometry imaging for single-state
//  signals. SAMpower reads SAM weights for volume image, & computes
//  the source waveform for each voxel. SAMpower computes 3D + time
//  images of output (S/N power) (pseudo-Z). The output is the change
//  in power relative to the mean power.
//
//  flags:
//      -r <dataset name>   -- MEG dataset name
//      -d <data file name> -- MEG data file name (4D, only)
//      -m <analysis name>  -- name of analysis file
//      -v                  -- verbose mode
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
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
#include <byteswap.h>
#include <samlib.h>
#include <siglib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <filters.h>
#include <TimeSegs.h>
#include <version.h>
#include "samutil.h"
#include "sam_param.h"
#include "sam_parse.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// internal constants
#define MINOR_REV   1           // new command line parser
#define SIG_ORDER   4
#define NYQUIST     0.5

int main(
    int             argc,
    char            **argv
) {
    IIRSPEC         DataIIR;            // IIR filter structure
    FFTSPEC         DataFFT;            // FFT filter structure
    IIRSPEC         SmoothIIR;          // IIR smoothing filter for Hilbert envelope
    FFTSPEC         SmoothFFT;          // FFT smoothing filter structure
    FFTSPEC         NotchFFT;           // FFT notch filter structure (sum of notches)
    FFTWPLAN        *fftwplan;          // fftw plans and work arrays
    PARMINFO        Params;             // analysis parameters
    PARM            *p;                 // pointer into active parameter list
    HeaderInfo      Header;             // MEG data header
    ChannelInfo     *Channel;           // MEG channel info
    EpochInfo       *Epoch;             // MEG epoch info
    nifti_1_header  NiiHdr;             // NIFTI header
    unsigned char   *ExtHdr;            // NiiHdr[N] -- extended NIFTI header
    gsl_matrix      *Wgt;               // Wgt - VxM active SAM weights
    gsl_vector      *W;                 // W - Mx1 SAM weight vector
    //GIIATLAS        *atlas;             // gifti atlas info
    double          **x;                // x[TT][M] -- filtered multisensor data
    double          *In;                // In[TT] -- input time-series
    double          *Out;               // Out[TT] -- output time-series
    double          *da;                // double array
    double          time;               // trial time duration (what else?)
    double          Noise;              // mean-square noise
    double          Nse;                // voxel noise
    double          BWFreq;             // effective bandwidth
    double          percent;            // percent done
    double          mean;               // sample mean
    double          scale;              // bits --> Tesla conversion
    register double tmp;                // what else?
    float           **Image;            // Image[D][V] -- image voxels
    int             E;                  // number of good epochs
    int             M;                  // number of MEG or reference channels used
    int             N;                  // size of extended header
    int             T;                  // number of samples
    int             TT;                 // total number of samples E*T
    int             V;                  // number of voxels for active state weight file
    int             c;                  // input character & virtual channel
    int             pct;                // integer percentage
    int             m;                  // channel index
    int             n;                  // alternate channel index
    int             d;                  // SAM image index
    int             e;                  // epoch index
    int             f;                  // 1st primary sensor index
    int             t;                  // sample index
    int             tt;                 // alternate sample index
    int             i;                  // matrix index
    int             tid;                // thread id
    int             NumImg;             // number of SAM images
    int             eflg = FALSE;       // command-line error flag
    int             mflg = FALSE;       // analysis parameter file flag
    int             rflg = FALSE;       // dataset name flag
    int             gflg = FALSE;       // gifti format Atlas
    int             vflg = FALSE;       // verbose mode flag
    char            fpath[256];         // general path name
    char            *DSName = NULL;     // MEG dataset name
    char            DSpath[256];        // MEG dataset path
    char            SAMpath[256];       // SAM subdirectory path
    char            AtlasPath[256];     // atlas path
    char            WgtDir[256];        // path to weights
    char            *WtsName = NULL;    // used in weights file names
    char            *OutName = NULL;    // used in output file names
    char            Prefix[64];         // incomplete pathname
    char            FilterName[8];      // filter name -- IIR | FFT
    char            extension[4];       // extension
    FILE            *fp;                // data output file pointer

    // system-dependent variables
#if BTI
    static int      dflg = FALSE;       // pdf file name flag
    char            *PDFName;           // name of data file
#else
    unsigned char   **Bad;              // Bad[E][T] -- flags for bad samples
#endif

    // Register the parameters used by this program.

    reg_usage(argv[0], MINOR_REV, __DATE__);
    reg_std_parm();

    reg_parm("Marker");
    reg_parm("SegFile");
    reg_parm("CovType");
    reg_parm("CovBand");
    reg_parm("ImageBand");
    reg_parm("SmoothBand");
    reg_parm("FilterType");
    reg_parm("Notch");
    reg_parm("Hz");
    reg_parm("AtlasName");
    reg_parm("ImageFormat");
    reg_parm("ImageMetric");
    set_parm_arg_help("ImageMetric", "Hilbert");
    reg_parm("TimeStep");
    reg_parm("RemoveBaseline");
    reg_parm("ImageDirectory");
    reg_parm("PrefixLength");
    reg_parm("WtsName");
    reg_parm("OutName");

    // Parse the arguments, environment variables, and parameter files.

    do_parse_args(argc, argv);

    // Extract specified parameters.

    eflg = get_params(&Params);

    p = get_parm("DataSet");
    if (p->set) {
        rflg = TRUE;
        DSName = copy_string(Params.DataSetName);
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
    if(eflg || !rflg || !dflg || !mflg) {
#else
    if(eflg || !rflg || !mflg) {
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

    if (vflg) {
        printf("S A M   P O W E R   I M A G I N G   O F   S O U R C E   T I M E   S E R I E S\t%s\n", PRG_REV);
        printf("opening data file");
        fflush(stdout);
    }

    // get data with sensor structures
#if BTI
    if (GetMagnesInfo(DSName, PDFName, &Header, &Channel, &Epoch) == -1)
        Cleanup("can't read Magnes information");
    sprintf(DSpath, "%s/%s", Header.DsPath, Header.SetName);
#else
    GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);
    sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
#endif
    sprintf(SAMpath, "%s/SAM", DSpath);

    // set constants
    M = Header.NumPri;
    E = Header.NumEpochs;
    T = Header.MaxSamples;
    TT = E * T;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index

    // check to see if this is a continuous dataset
    if (E != 1)
        fprintf(stderr, " (reminder: multi-epoch data must be contiguous)");

    if (vflg) {
        printf(" - done\n");
        printf("'%s' uses %d epochs, %d channels, %d samples/channel/epoch\n", DSName, E, M, T);
    }

    // check imaging parameters
    if (Params.CovHP == -999. || Params.CovLP == -999.)
        Cleanup("CovBand parameters required");
    if (Params.ImageHP == -999. || Params.ImageLP == -999.)
        Cleanup("ImageBand parameters required");
    if (Params.ImageHP < Channel[f].AnalogHpFreq)
        Cleanup("ImageBand highpass cutoff frequency out of recorded range");
    if (Params.ImageLP > Channel[f].AnalogLpFreq)
        Cleanup("ImageBand lowpass cutoff frequency out of recorded range");
    if (Params.SmoothHP == -999. || Params.SmoothLP == -999.)
        Cleanup("SmoothBand specifications required");
    if (Params.SmoothHP < Channel[f].AnalogHpFreq)
        Cleanup("SmoothBand highpass cutoff frequency out of recorded range");
    if (Params.SmoothLP > Channel[f].AnalogLpFreq)
        Cleanup("SmoothBand lowpass cutoff frequency out of recorded range");
    if (Params.ImageMetric != HILBERT)
        Cleanup("ImageMetric must specify Hilbert");
    if (Params.TimeStep == -999.)
        Cleanup("TimeStep specification missing from parameter file");
    if (Params.CovType != GLOBAL_)
        Cleanup("CovType must be GLOBAL");
    if (Params.FilterType == IIR)
        strcpy(FilterName, "IIR");
    else
        strcpy(FilterName, "FFT");
    if (vflg) {
        printf(" - done\n");
        switch (Params.RemoveBaseline) {
            case NONE: printf("output Hilbert envelope of power images without baseline removal\n"); break;
            case BYVOXEL: printf("output Hilbert envelope of power images with mean removal by voxel\n"); break;
            case GLOBAL: printf("output Hilbert envelope of power images with global mean removed from each voxel\n"); break;
            default: break;
        }
        fflush(stdout);
    }

    // set prefix to NumPrefix characters (default is 8-character dataset hashcode)
    memset(Prefix, 0, sizeof(Prefix));
#if BTI
    memcpy(Prefix, "BIOMAG4D", 8);  // 4D datasets just have serial numbers
#else
    memcpy(Prefix, DSName, Params.NumPrefix);
#endif

    // check for an atlas

    p = get_parm("AtlasName");
    if (p->set) {
        GetMRIPath(AtlasPath, sizeof(AtlasPath), &Params, (char *)p->ptr, TRUE);
        Params.ImageFormat = ORIG;              // force ORIG if atlas is present

        // check for a GIFTI atlas, in which case the Wgts will be GIFTI

        i = strlen(AtlasPath);
        if (i > 4 && strcmp(&AtlasPath[i - 4], ".gii") == 0) {
            gflg = TRUE;
            // get the atlas info
            //atlas = GetGiiAtlas(AtlasPath);
        }
    }

    // read SAM weights
    if (vflg) {
        printf("reading Global SAM weights & noise estimate");
        fflush(stdout);
    }
    sprintf(WgtDir, "%s/%s,%-d-%-dHz", SAMpath, WtsName, (int)Params.CovHP, (int)Params.CovLP);
    GetNoise(WgtDir, "Global", &Noise);
    if (Params.ImageFormat == TLRC)
        sprintf(fpath, "%s/Global_at.nii", WgtDir);
    else
        sprintf(fpath, "%s/Global.nii", WgtDir);
    if (gflg) {
        Wgt = GetGIFTIWts(fpath);
    } else {
        Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
    }
    V = Wgt->size1;
    if (Wgt->size2 != M)
        Cleanup("number of channels in dataset does not match number in weight file");

    // if req'd, summarize epoch dimensions
    if (vflg) {
        printf(" - done\n");
        printf("dataset '%s' has %d epochs\n", DSName, E);
        printf("\tduration %6.3f seconds/epoch\n", (double)T/Header.SampleRate);
        fflush(stdout);
    }

    // set up 3d+time imaging constants)
    time = (double)TT / Header.SampleRate - 2.; // total duration (seconds) -- less 1 second at start & end
    NumImg = (int)rint(time / Params.TimeStep); // number of 3d images

    // allocate memory
    if (vflg) {
        printf("solving for %d SAM images\n", NumImg);
        printf("allocating memory");
        fflush(stdout);
    }
    W = gsl_vector_alloc(M);
    x = new_matrixE(TT, M, "x[][]");
    In = new_arrayE(double, TT, "In[]");
    Out = new_arrayE(double, TT, "Out[]");
    Image = new_arrayE(float *, NumImg, "Image[]");
    for (d=0; d<NumImg; d++)
        Image[d] = new_arrayE(float, V, "Image[][]");

    // determine filter type
    if (vflg) {
        printf(" - done\n");
        printf("making %6.3f to %6.3f Hz %s filter for weight bandpass", Params.ImageHP, Params.ImageLP, FilterName);
        fflush(stdout);
    }
    if (Params.FilterType == IIR)
        mkfilter(&DataIIR, SIG_ORDER, Header.SampleRate, Params.ImageHP, Params.ImageLP, &BWFreq);
    else
        Butterworth(&DataFFT, TT, Params.ImageHP, Params.ImageLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

    // generate IIR smoothing filter for Hilbert time-series
    if (vflg) {
        printf(" - done\n");
        printf("making %6.3f Hz lowpass %s filter for smoothing Hilbert envelope of power", Params.SmoothLP, FilterName);
        fflush(stdout);
    }
    if (Params.FilterType == IIR)
        mkfilter(&SmoothIIR, SIG_ORDER, Header.SampleRate, 0., Params.SmoothLP, &BWFreq);
    else
        Butterworth(&SmoothFFT, TT, 0., Params.SmoothLP, Header.SampleRate, LOWPASS, MAXORDER, &BWFreq);

    // compute notch filter
    if (Params.Notch) {
        if (Params.FilterType == FFT) {
            if (vflg) {
                printf(" - done\n");
                printf("making FFT notch filter");
                fflush(stdout);
            }
            mknotch(&NotchFFT, Params.CovHP, Params.CovLP, Header.SampleRate, TT, Params.Hz);
        } else {
            msg("\nreminder: when using IIR filters, data should be notch filtered for power mains prior to analysis\n");
        }
    }

    // read & filter MEG data
    if (vflg) {
        printf(" - done\n");
        printf("reading & filtering MEG data");
        fflush(stdout);
    }
    for (m=0; m<M; m++) {                           // read & process data in channel order
        n = Header.PsIndex[m];
        scale = Channel[n].Scale;
        for (e=0; e<E; e++) {                       // multi-epoch data is presumed to be contiguous!
#if BTI
            if (GetMagnesChan(DSpath, PDFName, &Header, Channel, e, m, In) == -1)
                Cleanup("'GetMagnesChan()' failed");
#else
            GetDsData(&Header, e, n, scale, In);    // get one epoch of data for this channel -- length T samples
#endif
            for (t=0; t<T; t++) {
                tt = e * T + t;
                Out[tt] = In[t];
            }
        }
        demean(Out, TT);
        da = In;        // filter result
        if (Params.FilterType == IIR) {
            bdiir(Out, In, TT, &DataIIR);
        } else {
            FFTfilter(Out, In, &DataFFT, TT);
            if (Params.Notch) {
                da = Out;
                FFTfilter(In, Out, &NotchFFT, TT);
            }
        }
        for (t=0; t<TT; t++) {
            x[t][m] = da[t];
        }

        // heartbeat
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }

    // for each latency step, compute one volumetric SAM image...
    if (vflg) {
        printf(" - done\n");
        printf("computing SAM images for each time step:\n");
        fflush(stdout);
    }

#ifdef _OPENMP
    tid = omp_get_thread_num();
msg("tid = %d\n", tid);
    // omp_get_num_threads
#endif

    // set up fftw plans for hilbert()

    fftwplan = new_fftwplan(TT, TRUE);

    // for each voxel...
    pct = -1;
    for (c=0; c<V; c++) {

        // test weight vector to determine if within head boundary
        gsl_matrix_get_row(W, Wgt, c);
        if (gsl_vector_isnull(W) == 0) {    // compute image only for non-zero weights

            // compute projected noise variance
            for (m=0, tmp=0.; m<M; m++)
                tmp += gsl_vector_get(W, m) * Noise * gsl_vector_get(W, m);
            Nse = sqrt(tmp);

            // compute virtual sensor time-series
            for (t=0; t<TT; t++) {
                for (m=0, tmp=0.; m<M; m++)
                    tmp += x[t][m] * gsl_vector_get(W, m);
                In[t] = tmp / Nse;          // In[] is the virtual sensor source time-series (S/N)
            }

            // compute Hilbert envelope of output power (S/N)
            hilbert(fftwplan, In, Out);     // Out[] is the phase-shifted VS time-series
            for (t=0; t<TT; t++)
                In[t] = (In[t] * In[t] + Out[t] * Out[t]);  // In[] is the unsmoothed Hilbert envelope of power

            // smooth Hilbert envelope & downsample
            if (Params.FilterType == IIR)
                bdiir(In, Out, TT, &SmoothIIR);
            else
                FFTfilter(In, Out, &SmoothFFT, TT);         // lowpass to smooth moment

            // remove filtering artifacts
            for (t=0; t<TT; t++)
                Out[t] = (Out[t] < 0.) ? 0. : Out[t];

            // downsample envelope of power
            for (d=0, time=1.; d<NumImg; d++, time+=Params.TimeStep) {
                t = (int)rint(time * Header.SampleRate);
                Image[d][c] = (float)Out[t];
            }

        } else {
            for (d=0; d<NumImg; d++)
                Image[d][c] = 0.;
        }

        // percent done (counting completed SAM images)
        if (vflg) {
            percent = 100. * (double)c / (double)V;
            if (rint(percent) > (double)pct) {
                pct = (int)rint(percent);
                printf("\r\t%d percent done", pct);
                fflush(stdout);
            }
        }
    }

    free_fftwplan(fftwplan);

    // if req'd, remove voxel baselines
    switch (Params.RemoveBaseline) {

        case NONE:      // no baseline removal
            break;

        case BYVOXEL:   // remove baseline for each voxel
            for (c=0; c<V; c++) {
                for (d=0, mean=0.; d<NumImg; d++)
                    mean += (double)Image[d][c];
                mean /= (double)NumImg;
                for (d=0; d<NumImg; d++)
                    Image[d][c] -= (float)mean;
            }
            break;

        case GLOBAL:    // remove global mean baseline from each voxel
            // find mean global power (S/N)
            for (d=i=0, mean=0.; d<NumImg; d++)
                for (c=0; c<V; c++)
                    if (Image[d][c] > 0.) { // non-zero voxels, only!
                        mean += (double)Image[d][c];
                        i++;                // count non-zero voxels!
                    }
            mean /= (double)(i * NumImg);

            // subtract mean
            for (d=0; d<NumImg; d++)
                for (c=0; c<V; c++)
                    if (Image[d][c] > 0.)
                        Image[d][c] -= (float)mean;
            break;

        default:
            break;
    }

    // create NIFTI file name
    if (!strncmp(Params.DirName, "NULL", 4))
        sprintf(fpath, "%s/%s,%s,%-d-%-dHz,PWR.nii", SAMpath, Prefix, OutName, (int)Params.ImageHP, (int)Params.ImageLP);
    else
        sprintf(fpath, "%s/%s,%s,%-d-%-dHz,PWR.nii", Params.DirName, Prefix, OutName, (int)Params.ImageHP, (int)Params.ImageLP);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open .nii file for write");

    // modify NIFTI header
    memset(extension, 0, 4);                                // for now, ignore the fact that we may have an extended header
    NiiHdr.dim[0] = 4;                                      // 3d+time
    NiiHdr.dim[4] = NumImg;                                 // number of 3d images
    NiiHdr.toffset = 0.;                                    // time offset shift
    NiiHdr.slice_duration = (float)Params.TimeStep;         // time step size
    NiiHdr.pixdim[4] = (float)Params.TimeStep;
    NiiHdr.slice_code = NIFTI_SLICE_SEQ_INC;
    NiiHdr.xyzt_units = (NIFTI_UNITS_MM | NIFTI_UNITS_SEC); // 3d+time units
    memset(NiiHdr.intent_name, 0, 16);
    memcpy(NiiHdr.intent_name, "Power", 5);

    // write header & extension
    if (fwrite((void *)&NiiHdr, sizeof(nifti_1_header), 1, fp) != 1)
        Cleanup("can't write nifti header to .nii file");
    if (fwrite((void *)extension, 4, 1, fp) != 1)
        Cleanup("can't write 'extension' to .nii file");
    if ((int)NiiHdr.vox_offset != 352 && ExtHdr != NULL) {
        N = (int)NiiHdr.vox_offset - 352;
        if (fwrite((void *)ExtHdr, N, 1, fp) != 1)
            Cleanup("can't write extended header to .nii file");
    }
    for (d=0; d<NumImg; d++)
        if (fwrite((void *)Image[d], sizeof(float), V, fp) != V)
            Cleanup("can't write data image file");
    fclose(fp);

    // that's all, folks!
    log_params(SAMpath);

    if (vflg) {
        msg(" - done\n'%s' done\n", Progname);
    }
    exit(0);
}
