// SAM4d -- synthetic aperture magnetometry imaging for source event related
//  time-locked signals. SAM4d reads all SAM weights for volume image, & computes
//  the source waveform for each voxel. Given specifications for active & control
//  time segments, SAM4d computes the mean power & variance time-series relative
//  to the specified event windows. The result is represented stored in a NIFTI image file
//  '.nii' file.
//
//  flags:
//      -r <dataset name>       -- MEG dataset name
//      -d <data file name>     -- MEG data file name (4D, only)
//      -m <parameter name>     -- name of parameter file
//      -v                      -- verbose mode
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

// include...
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <fcntl.h>
#include <gsl/gsl_matrix.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <rvelib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <filters.h>
#include <TimeSegs.h>
#include <nifti1.h>
#include <version.h>
#include "samutil.h"
#include "sam_param.h"
#include "sam_parse.h"

// internal constants
#define MINOR_REV   1
#define SIG_ORDER   4

typedef struct {
    int     Epoch;                  // epoch/trial
    int     TS;                     // starting sample
} SAM_SEG;

typedef struct {
    char        Name[64];           // marker file name
    gsl_matrix  *Wgt;               // Wgt - VxM SAM weights
    double      Noise;              // noise power
    SAM_HDR     WgtHdr;             // SAM weight file header
    SAM_SEG     *Segment;           // Segment[NumSegments] -- time-segment specs
    int         NumSegments;        // number of segments for this named marker
} IMAGE_STAT;

// INCOMPLETE: Need to work on DataSegments, integration, & time steps

int main(
    int     argc,
    char    **argv
) {
    IIRSPEC         DataIIR;        // IIR filter structure
    FFTSPEC         DataFFT;        // FFT data filter structure
    FFTSPEC         NotchFFT;       // notch FFT filter structure
    IMAGE_STAT      *Stats;         // Stats[N] -- structure for markers
    PARMINFO        Params;         // analysis parameters
    PARM            *p;             // pointer into active parameter list
    SAM_MARKS       *Marker;        // Marker[N] -- array of time-ordered markers
    HeaderInfo      Header;         // MEG data header
    ChannelInfo     *Channel;       // MEG channel info
    EpochInfo       *Epoch;         // MEG epoch info
    COV_HDR         CovHdr;         // noise covariance header
    nifti_1_header  NiiHdr;         // NIFTI header
    unsigned char   *ExtHdr;        // NiiHdr[N] -- extended NIFTI header
    gsl_matrix      *Cn;            // Cn - MxM noise covariance matrix
    gsl_vector      *W;             // W - Mx1 weight vector
    GIIATLAS        *atlas;         // gifti atlas info
    double          ***x;           // x[E][T][M] -- filtered multisensor data
    double          *In;            // In[T] -- input time-series
    double          *Out;           // Out[T] -- output time-series
    double          *Avg;           // Avg[NumImg] -- average
    double          *Var;           // Var[NumImg] -- variance
    double          *da;            // double array
    double          width;          // duration (sec)
    double          mark;           // time mark
    double          offset;         // offset from segment start to average window start
    double          BSoffset;       // baseline offset
    double          BWFreq;         // effective bandwidth
    double          BWratio;        // ratio of covariance bandwidth to image bandwidth (nominal, only)
    double          Noise;          // effective noise
    double          percent;        // percent done
    double          a1;             // sum of metric
    double          a2;             // sum of metric squared
    double          mean;           // mean baseline
    double          **sample;       // segment source moment
    float           **ImgM;         // ImgM[D][V] -- image mean
    float           **ImgV;         // ImgV[D][V] -- image variance
    double          tmp;            // what else?
    int             *ChanIndex;     // ChanIndex[M] -- index of channels used
    int             E;              // number of good epochs
    int             M;              // number of MEG or reference channels used
    int             N;              // just a number to take care of extended header read
    int             T;              // number of samples
    int             TS;             // data segment start
    int             TE;             // data segment end
    int             BS;             // baseline start sample
    int             BE;             // baseline end sample
    int             PS;             // polarity starting sample
    int             PE;             // polarity ending sample
    int             TSEG;           // number of samples in data segment
    int             V;              // number of voxels weight file
    int             v;              // voxel index
    int             c;              // virtual channel
    int             pct;            // integer percentage
    int             m;              // channel index
    int             d;              // image step index
    int             e;              // epoch index
    int             f;              // 1st primary channel index
    int             t;              // sample index
    int             tt;             // alternate sample index
    int             i;              // matrix index
    int             j;              // ...and another
    int             k;              // this is the last one -- for sure...
    int             n;              // condition index
    int             NumMarkers;     // number of markers in dataset -- e.g., 'MarkerFile.mrk'
    int             NumConditions;  // number of conditions in output image
    int             NumTot;         // total number of named markers
    int             NumImg;         // number of SAM images
    int             NumBox;         // number of boxcar integrator samples
    int             *vidx;          // array for shuffling voxel indices
    int             bflg = FALSE;   // baseline removal flag
    int             eflg = FALSE;   // command-line error flag
    int             mflg = FALSE;   // marker parameter file name flag
    int             nflg = FALSE;   // use noise covariance file
    int             rflg = FALSE;   // run name flag
    int             gflg = FALSE;   // gifti format Atlas
    int             pflg = FALSE;   // polarity test
    int             vflg = FALSE;   // verbose mode flag
    time_t          secs;           // random number seed
    char            *s;             // general string
    char            fpath[256];     // general path name
    char            *DSName = NULL; // MEG dataset name
    char            DSpath[256];    // MEG dataset path
    char            SAMpath[256];   // SAM subdirectory path
    char            AtlasPath[256]; // atlas path
    char            WgtDir[256];    // weight directory path name
    char            *CovName = NULL; // name used in covariance file names
    char            *WtsName = NULL; // name used in weights file names
    char            *OutName = NULL; // name used in output file names
    char            Prefix[64];     // incomplete pathname
    char            Name1[256];     // partial image name
    char            Name2[256];     // partial image name
    char            FilterName[8];  // filter name -- IIR | FFT
    char            extension[4];   // extension
    FILE            *fp;            // data output file pointer

    //system-dependent variables
#if BTI
    static int      dflg = FALSE;           // pdf file name flag
    char            *PDFName;               // name of data file
#else
    double          scale;                  // bits-->Tesla conversion
    unsigned char   **Bad;                  // Bad[E][T] -- flags for bad samples
#endif

    // Register the parameters used by this program.

    reg_usage(argv[0], MINOR_REV, __DATE__);
    reg_std_parm();

    reg_parm("Marker");
    reg_parm("DataSegment");
    reg_parm("SignSegment");
    reg_parm("Baseline");
    reg_parm("CovType");
    reg_parm("CovBand");
    reg_parm("ImageBand");
    reg_parm("SmoothBand");
    reg_parm("NoiseBand");
    reg_parm("TimeStep");
    reg_parm("BoxWidth");
    reg_parm("FilterType");
    reg_parm("Notch");
    reg_parm("Hz");
    reg_parm("AtlasName");
    reg_parm("ImageFormat");
    reg_parm("ImageMetric");
    set_parm_arg_help("ImageMetric", "Moment|Power|RankVectorEntropy");
    reg_parm("ImageDirectory");
    reg_parm("PrefixLength");
    reg_parm("CovName");
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
    CovName = OutName;
    p = get_parm("CovName");
    if (p->set) {
        CovName = (char *)p->ptr;
    }
    WtsName = CovName;
    p = get_parm("WtsName");
    if (p->set) {
        WtsName = (char *)p->ptr;
    }

    if (vflg) {
        printf("4 D   S A M   I M A G I N G\t%s\n", PRG_REV);
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

    // count data dimensions
    // set constants
    M = Header.NumPri;
    E = Header.NumEpochs;
    T = Header.MaxSamples;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index

    if (vflg) {
        printf(" - done\n");
        printf("'%s' uses %d epochs, %d channels, %d samples/channel/epoch\n", DSName, E, M, T);
    }

    // sanity check
    if (M == 0 || E == 0 || T == 0)
        Cleanup("no valid data");

    // set the random number seed
    time(&secs);
    srandom((unsigned int)secs);

    if (Params.NumMark < 1)
        Cleanup("There must be at least one marker.");      // @@@ maybe only 1

    // announce metric type
    if (Params.ImageMetric == -1)
        Cleanup("parameters must specify 'ImageMetric'");
    if (vflg) {
        printf(" - done\n");
        printf("Imaging functions of ");
        switch (Params.ImageMetric) {
            case MOMENT: printf("moment\n"); break;
            case POWER: printf("power\n"); break;
            case RV_ENTROPY: printf("rank vector entropy\n"); break;
            default: Cleanup("selected ImageMetric not supported"); break;
        }
    }

    // check imaging parameters
    if (Params.BaseStart != -999. && Params.BaseEnd != -999.)
        bflg = TRUE;
    if (Params.CovType == -1)
        Cleanup("'CovType' is mandatory");
    if (Params.CovHP == -999. || Params.CovLP == -999.)         // CovBand is always required!
        Cleanup("'CovBand' specifications required");
    if (Params.ImageHP == -999. || Params.ImageLP == -999.)     // ImageBand is always required!
        Cleanup("'ImageBand' specifications required");
    if (Params.ImageHP < Channel[f].AnalogHpFreq)
        Cleanup("ImageBand highpass cutoff frequency outside recorded range");
    if (Params.ImageLP > Channel[f].AnalogLpFreq)
        Cleanup("ImageBand lowpass cutoff frequency outside recorded range");
    if (Params.TimeStep < 0.)
        Cleanup("must specify 'TimeStep'");
    BWratio = (Params.ImageLP - Params.ImageHP) / (Params.CovLP - Params.CovHP);
    if (BWratio > 1.)
        Cleanup("covariance bandwidth must be >= image bandwidth");
    if (Params.FilterType == IIR)
        strcpy(FilterName, "IIR");
    else
        strcpy(FilterName, "FFT");

    // data segment extends the start & end times for each trial
    // allowing filtering and entropy computations on short time-series

    width = 0.;
    p = get_parm("BoxWidth");
    if (p->set) {
        width = *(double *)p->ptr;
    }
    p = get_parm("SmoothBand");
    if (p->set) {
        if (width != 0.) {
            fatalerr("Please specify only one of BoxWidth or SmoothBand.");
        }
        width = 1. / Params.SmoothLP;
    }
    NumBox = (int)rint(width * Header.SampleRate);
    if (NumBox <= 0) {
        fatalerr("Please specify either BoxWidth or SmoothBand.");
    }

        // @@@ This always uses marker 0!

    p = get_parm("DataSegment");
    if (!p->set) {
        // if DataSegment not specified, use the marker segment
        Params.DataSegStart = Params.Marker[0].MarkStart;
        Params.DataSegEnd = Params.Marker[0].MarkEnd;
    }

    // set up analysis constants
    width = Params.Marker[0].MarkEnd - Params.Marker[0].MarkStart;  // duration of average window (for all markers!)
    width -= NumBox / Header.SampleRate;                            // subtract one box
    NumImg = (int)rint(width / Params.TimeStep) + 1;                // number of 3d images
    width = Params.DataSegEnd - Params.DataSegStart;                // duration of segment analysis window
    TSEG = (int)rint(width * Header.SampleRate);                    // number of samples in segment analysis window
    offset = Params.Marker[0].MarkStart - Params.DataSegStart;
    if (offset < 0)
        Cleanup("response time window must begin after segment start");
    if (bflg) {
        TS = (int)rint(Params.DataSegStart * Header.SampleRate);    // start of segment
        TE = (int)rint(Params.DataSegEnd * Header.SampleRate);      // end of segment
        BS = (int)rint(Params.BaseStart * Header.SampleRate);       // baseline starting sample relative to trigger marker (can be negative!)
        BE = (int)rint(Params.BaseEnd * Header.SampleRate);         // baseline ending sample relative to trigger marker (can be negative!)
        if (BS < TS || BS > TE || BE < TS || BE > TE)
            Cleanup("baseline starting / ending samples out of range");
        BSoffset = Params.Marker[0].MarkStart;
    }

    // time window for the polarity test

    p = get_parm("SignSegment");
    if (p->set) {
        pflg = TRUE;
        BSoffset = Params.Marker[0].MarkStart;
        PS = (int)rint((((double *)p->ptr)[0] - BSoffset) * Header.SampleRate);
        PE = (int)rint((((double *)p->ptr)[1] - BSoffset) * Header.SampleRate);
        if (PS < TS || PS > TE || PE < TS || PE > TE)
            Cleanup("polarity segment starting & ending samples out of range");
    }

    // set prefix to NumPrefix characters (default is 8-character dataset hashcode)
    memset(Prefix, 0, sizeof(Prefix));
#if BTI
    memcpy(Prefix, "BIOMAG4D", 8);              // 4D datasets just have serial numbers
#else
    memcpy(Prefix, DSName, Params.NumPrefix);   // default: copy eight-letter NIMH hash code in dataset name
#endif

    // if req'd, initialize NOISE
    if (Params.NoiseHP != -999. && Params.NoiseLP != -999.) {
        if (Params.NoiseHP < Channel[f].AnalogHpFreq)
            Cleanup("NoiseBand highpass cutoff frequency outside recorded range");
        if (Params.NoiseLP > Channel[f].AnalogLpFreq)
            Cleanup("NoiseBand lowpass cutoff frequency outside recorded range");
        if (vflg) {
            printf("reading noise covariance file");
            fflush(stdout);
        }
        sprintf(fpath, "%s/%s,%-d-%-dHz/Noise.cov", SAMpath, CovName, (int)Params.CovHP, (int)Params.CovLP);
        Cn = gsl_matrix_alloc(M, M);
        GetCov(fpath, &CovHdr, &ChanIndex, Cn);
        nflg = TRUE;
        if (vflg)
            printf(" - done\n");
    }

    // if req'd, summarize epoch dimensions
    if (vflg) {
        printf("dataset '%s' has %d epochs\n", DSName, E);
        printf("\tduration %6.3f seconds/epoch\n", (double)T / Header.SampleRate);
    }

    // if req'd, read markers
    if (Params.NumMark > 0) {
#if BTI
        if (vflg) {
            printf("parsing trigger markers");
            fflush(stdout);
        }
        if (GetMagnesMarkers(DSpath, PDFName, &Header, Channel, &Marker, &NumMarkers) == -1)
            Cleanup("cannot get trigger markers");
#else
        if (vflg) {
            printf("reading 'MarkerFile.mrk'");
            fflush(stdout);
        }
        sprintf(fpath, "%s/MarkerFile.mrk", DSpath);
        GetMarkers(fpath, &Marker, &NumMarkers);
#endif
        if (vflg) {
            printf(" - found %d markers", NumMarkers);
            fflush(stdout);
        }

        // allocate one Stats structure per marker
        NumConditions = Params.NumMark;
        Stats = new_arrayE(IMAGE_STAT, NumConditions, "Stats[]");
        for (n=0; n<NumConditions; n++) {
            s = Params.Marker[n].MarkName2;
            if (s == NULL) {
                s = Params.Marker[n].MarkName;
            }
            strcpy(Stats[n].Name, s);
        }
    }

    // populate segment array
    if (vflg) {
        printf(" - done\n");
        printf("populating segments:\n");
    }

    // for each marker...
    for (n=0; n<NumConditions; n++) {
        if (vflg) {
            printf("\t'%s'\n", Stats[n].Name);
        }

        // count the total markers for each named marker
        s = Params.Marker[n].MarkName;
        for (i=NumTot=0; i<NumMarkers; i++)
            if (!strcmp(s, Marker[i].Name))
                NumTot++;
        if (NumTot == 0)
            Cleanup("no markers found with specified names");

        // allocate segment array
        Stats[n].Segment = new_arrayE(SAM_SEG, NumTot, "Stats[].Segment[]");

        // populate Stats for each condition
        for (i=k=0, e=(-1); i<NumMarkers; i++) {
            if (Marker[i].Epoch != e)
                e = Marker[i].Epoch;
            if (Epoch[e].Flag) {

                // populate segment windows
                if (!strcmp(s, Marker[i].Name)) {
                    mark = Marker[i].Latency;
                    TS = (int)rint((mark + Params.DataSegStart) * Header.SampleRate) + Header.PreTrigSamples;
                    TE = (int)rint((mark + Params.DataSegEnd) * Header.SampleRate) + Header.PreTrigSamples;
                    if (TS >= 0 && TE < T) {
                        Stats[n].Segment[k].Epoch = e;
                        Stats[n].Segment[k].TS = TS;
                        k++;
                    } else {
                        msg("\x1b[7mINFO:\x1b[m marker %s skipped, trial %d, t=%g\n", s, e+1, mark);
                    }
                }
            }
        }
        if (k == 0)
            Cleanup("there were no valid data segments for specified markers & time windows");
        Stats[n].NumSegments = k;
    }

    // check for an atlas

    p = get_parm("AtlasName");
    if (p->set) {
        GetMRIPath(AtlasPath, sizeof(AtlasPath), &Params, (char *)p->ptr, TRUE);
        Params.ImageFormat = ORIG;              // @@@ force ORIG if atlas is present

        // check for a GIFTI atlas, in which case the Wgts will be GIFTI

        i = strlen(AtlasPath);
        if (i > 4 && strcmp(&AtlasPath[i - 4], ".gii") == 0) {
            gflg = TRUE;
            // get the atlas info
            atlas = GetGiiAtlas(AtlasPath);
        }
    }

    // read weights according to CovType
    if (vflg) {
        printf("reading weights");
        fflush(stdout);
    }
    sprintf(WgtDir, "%s/%s,%-d-%-dHz", SAMpath, WtsName, (int)Params.CovHP, (int)Params.CovLP);
    switch (Params.CovType) {

        case GLOBAL_:       // read Global weights & Noise
            GetNoise(WgtDir, "Global", &Noise);
            if (Params.ImageFormat == TLRC)
                sprintf(fpath, "%s/Global_at.nii", WgtDir);
            else
                sprintf(fpath, "%s/Global.nii", WgtDir);
            for (n=0; n<NumConditions; n++) {
                if (gflg) {
                    Stats[n].Wgt = GetGIFTIWts(fpath);
                } else {
                    Stats[n].Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
                }
                V = Stats[n].Wgt->size1;
                if (Stats[n].Wgt->size2 != M)
                    Cleanup("number of channels in dataset does not match number in weight file");
                Stats[n].Noise = Noise;
            }
            break;

        case SUM_:          // read Sum weights
            GetNoise(WgtDir, "Sum", &Noise);
            if (Params.ImageFormat == TLRC)
                sprintf(fpath, "%s/Sum_at.nii", WgtDir);
            else
                sprintf(fpath, "%s/Sum.nii", WgtDir);
            for (n=0; n<NumConditions; n++) {
                if (gflg) {
                    Stats[n].Wgt = GetGIFTIWts(fpath);
                } else {
                    Stats[n].Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
                }
                V = Stats[n].Wgt->size1;
                if (Stats[n].Wgt->size2 != M)
                    Cleanup("number of channels in dataset does not match number in weight file");
                Stats[n].Noise = Noise;
            }
            break;

        case ALL_:      // read individual (All) weights
            for (n=0; n<NumConditions; n++) {
                GetNoise(WgtDir, Stats[n].Name, &Stats[n].Noise);
                if (Params.ImageFormat == TLRC)
                    sprintf(fpath, "%s/%s_at.nii", WgtDir, Stats[n].Name);
                else
                    sprintf(fpath, "%s/%s.nii", WgtDir, Stats[n].Name);
                if (gflg) {
                    Stats[n].Wgt = GetGIFTIWts(fpath);
                } else {
                    Stats[n].Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
                }
                V = Stats[n].Wgt->size1;
                if (Stats[n].Wgt->size2 != M)
                    Cleanup("number of channels in dataset does not match number in weight file");
            }
            break;

        default:
            Cleanup("unknown 'CovType'");
            break;
    }

    // allocate memory
    if (vflg) {
        printf(" - done\n");
        printf("allocating memory");
        fflush(stdout);
    }
    ImgM = new_arrayE(float *, NumImg, "ImgM[]");
    for (d=0; d<NumImg; d++)
        ImgM[d] = new_arrayE(float, V, "ImgM[][]");
    ImgV = new_arrayE(float *, NumImg, "ImgV[]");
    for (d=0; d<NumImg; d++)
        ImgV[d] = new_arrayE(float, V, "ImgV[][]");
    x = new_arrayE(double **, E, "x[]");
    for (e=0; e<E; e++) {
        x[e] = new_matrixE(T, M, "x[][]");
    }
    In = new_arrayE(double, T, "In[]");
    Out = new_arrayE(double, T, "Out[]");

    // compute filter
    if (vflg) {
        printf(" - done\n");
        printf("making %6.3f to %6.3f Hz %s filter for image bandpass", Params.ImageHP, Params.ImageLP, FilterName);
        fflush(stdout);
    }
    if (Params.FilterType == IIR)
        mkfilter(&DataIIR, SIG_ORDER, Header.SampleRate, Params.ImageHP, Params.ImageLP, &BWFreq);
    else
        Butterworth(&DataFFT, T, Params.ImageHP, Params.ImageLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

    // compute notch filter
    if (Params.Notch) {
        if (Params.FilterType == FFT) {
            if (vflg) {
                printf(" - done\n");
                printf("making FFT notch filter");
                fflush(stdout);
            }
            mknotch(&NotchFFT, Params.CovHP, Params.CovLP, Header.SampleRate, T, Params.Hz);
        } else {
            msg("\nreminder: when using IIR filters, data should be notch filtered for power mains prior to analysis\n");
        }
    }

    // summarize boxcar integrator & offset
    if (vflg) {
        printf(" - done\n");
        printf("boxcar integrator: %d samples (%f sec)\n", NumBox, (double)NumBox/Header.SampleRate);
        fflush(stdout);
    }

    // read & filter MEG data
    if (vflg) {
        printf(" - done\n");
        printf("reading & filtering MEG data");
        fflush(stdout);
    }
    for (e=0; e<E; e++) {
        for (m=0; m<M; m++) {
#if BTI
            if (GetMagnesChan(DSpath, PDFName, &Header, Channel, e, m, In) == -1)
                Cleanup("'GetMagnesChan()' failed");
#else
            n = Header.PsIndex[m];
            scale = Channel[n].Scale;
            GetDsData(&Header, e, n, scale, In);
#endif
            demean(In, T);
            da = Out;       // filter result
            if (Params.FilterType == IIR) {
                bdiir(In, Out, T, &DataIIR);
            } else {
                FFTfilter(In, Out, &DataFFT, T);
                if (Params.Notch) {
                    da = In;
                    FFTfilter(Out, In, &NotchFFT, T);
                }
            }
            for (t=0; t<T; t++) {
                x[e][t][m] = da[t];
            }
        }
        if (vflg) {         // if req'd - heartbeat
            printf(".");
            fflush(stdout);
        }
    }

    // run the voxel loop in random order, so the voxels with .ROI
    // set to false are spread out evenly across threads
    vidx = new_array(int, V);
    for (c = 0; c < V; c++) {
        vidx[c] = c;
    }
    shuffle(vidx, V);

    // for each named marker
    for (n=0; n<NumConditions; n++) {

        // announce condition
        if (vflg) {
            printf(" - done\n");
            printf("computing image voxels for '%s':\n", Stats[n].Name);
        }

        // loop on all coordinates in the ROI
        pct = -1;
        N = 0;
#pragma omp parallel private(W, da, Avg, Var, sample, i)
{
        W = gsl_vector_alloc(M);
        da = new_arrayE(double, TSEG, "da[]");
        Avg = new_arrayE(double, NumImg, "Avg[]");
        Var = new_arrayE(double, NumImg, "Var[]");
        // allocate metric time-series array per segment
        sample = new_arrayE(double *, Stats[n].NumSegments, "sample[]");
        for (i=0; i<Stats[n].NumSegments; i++)
            sample[i] = new_arrayE(double, NumImg, "sample[][]");

#pragma omp for private(c, m, tmp, j, k, a1, a2, mean, e, t, tt, TS, TE, Noise)
        for (v=0; v<V; v++) {
            c = vidx[v];

            // get voxel weight vector & determine if within head boundary
            gsl_matrix_get_row(W, Stats[n].Wgt, c);
            if (gsl_vector_isnull(W) == 0) {    // compute image only for non-zero weights

                // compute rms noise
                if (!nflg) {
                    for (m=0, tmp=0.; m<M; m++)
                        tmp += gsl_vector_get(W, m) * Stats[n].Noise * gsl_vector_get(W, m);
                    tmp *= BWratio;
                } else {
                    for (i=0, tmp=0.; i<M; i++)
                        for (j=0; j<M; j++)
                            tmp += gsl_vector_get(W, i) * gsl_matrix_get(Cn, i, j) * gsl_vector_get(W, j);
                    tmp *= (BWFreq / CovHdr.BWFreq);    // adjust noise power to image bandwidth
                }
                Noise = sqrt(tmp);

                // compute sam time-series for each segment
                for (k=0; k<Stats[n].NumSegments; k++) {

                    // set up segment initialization
                    e = Stats[n].Segment[k].Epoch;
                    tt = Stats[n].Segment[k].TS;

                    // compute SAM time series for segment
                    for (t=0; t<TSEG; t++, tt++) {
                        for (m=0, tmp=0.; m<M; m++)
                            tmp += x[e][tt][m] * gsl_vector_get(W, m);
                        da[t] = tmp / Noise;
                    }

                    // compute metric
                    switch (Params.ImageMetric) {
                        case MOMENT:        // event-related moment
                            for (i = 0; i < NumImg; i++) {
                                tmp = offset + i * Params.TimeStep;
                                t = (int)rint(tmp * Header.SampleRate);
                                TS = t;
                                TE = t + NumBox;
                                tmp = 0.;
                                for (t = TS; t < TE; t++) {
                                    tmp += da[t];
                                }
                                sample[k][i] = tmp / (TE - TS + 1);
                            }
                            break;
                        case POWER:         // power
                            for (i = 0; i < NumImg; i++) {
                                tmp = offset + i * Params.TimeStep;
                                t = (int)rint(tmp * Header.SampleRate);
                                TS = t;
                                TE = t + NumBox;
                                sample[k][i] = power(da, TS, TE);
                            }
                            break;
                        case RV_ENTROPY:    // rank vector entropy
#ifdef _OPENMP
                            fatalerr("OMP entropy unimplemented, recompile without -fopenmp");
#endif
                            // @@@ this is wrong anyway, sample[k] is not TSEG
                            fatalerr("entropy unimplemented");
                            Data2RVE(da, sample[k], TSEG, Params.Dims, Params.ImageLP, Params.Tau, Header.SampleRate);
                            break;
                        default:
                            Cleanup("ImageMetric not implemented");
                            break;
                    }

                }   // for(k...

                // compute average & variance
                for (i = 0; i < NumImg; i++) {
                    for (k=0, a1=a2=0.; k<Stats[n].NumSegments; k++) {
                        a1 += sample[k][i];
                        a2 += sample[k][i] * sample[k][i];
                    }
                    Avg[i] = (Stats[n].NumSegments > 2)? a1 / (double)Stats[n].NumSegments: 0.;
                    Var[i] = (Stats[n].NumSegments > 2)? (a2 - (a1 * a1 / (double)Stats[n].NumSegments)) / (double)(Stats[n].NumSegments - 1.): 0.;
                }

                // if req'd, remove baseline from Avg
                if (bflg) {
                    tmp = 0.;
                    j = 0;
                    for (i = 0; i < NumImg; i++) {
                        a1 = BSoffset + i * Params.TimeStep;
                        if (Params.BaseStart <= a1 && a1 <= Params.BaseEnd) {
                            tmp += Avg[i];
                            j++;
                        }
                    }
                    if (j == 0) {
                        fatalerr("No baseline time points in marker segment!");
                    }
                    mean = tmp / (double)j;
                    for (i = 0; i < NumImg; i++) {
                        Avg[i] -= mean;
                    }
                }

                // if req,d, correct polarity using atlas normals
                // @@@ lol that comment is funny, but we will implement that
                // @@@ meanwhile use signsegment
                if (pflg) {
                    for (i = 0; i < NumImg; i++) {
                        a1 = offset + i * Params.TimeStep;
                        j = (int)rint(a1 * Header.SampleRate);
                        if (PS <= j && j <= PE) {
                            tmp += Avg[i];
                        }
                    }
                    a2 = (tmp > 0.)? 1.: -1.;
                } else {
                    a2 = 1.;
                }

                // copy to output array

                for (i = 0; i < NumImg; i++) {
                    ImgM[i][c] = (float)(Avg[i] * a1);
                    ImgV[i][c] = (float)Var[i];
                }

            } else {    // zero weight vector
                for (i = 0; i < NumImg; i++) {
                    ImgM[i][c] = 0.;
                    ImgV[i][c] = 0.;
                }
            }

            // percent done
#pragma omp critical
            if (vflg) {
                percent = 100. * (double)N / (double)V;
                N++;
                if (rint(percent) > (double)pct) {
                    pct = (int)rint(percent);
                    printf("\r\t%d percent done", pct);
                    fflush(stdout);
                }
            }
        }       // for(c...

        free(da);
        gsl_vector_free(W);
        free(Avg);
        free(Var);
        for (i=0; i<Stats[n].NumSegments; i++)
            free(sample[i]);
        free(sample);
} // end parallel

        // write mean & variance images
        if (vflg) {
            printf("\nwriting 'Mean'");
            fflush(stdout);
        }
        if (strcmp(Params.DirName, "NULL") == 0)
            sprintf(Name1, "%s/%s,%s,%s,", SAMpath, Prefix, OutName, Stats[n].Name);
        else
            sprintf(Name1, "%s/%s,%s,%s,", Params.DirName, Prefix, OutName, Stats[n].Name);
        switch (Params.ImageMetric) {
            case MOMENT:        sprintf(Name2, "%s4D_MOM,", Name1); break;
            case POWER:         sprintf(Name2, "%s4D_PWR,", Name1); break;
            case RV_ENTROPY:    sprintf(Name2, "%s4D_RVE,", Name1); break;
            default:            Cleanup("ImageMetric not implemented"); break;
        }
        if (gflg) {
            /* Right hemisphere is stored first in ImgM! */

            sprintf(fpath, "%sMean.rh.gii", Name2);
            PutGIFTISurf(fpath, ImgM, NumImg, 0, atlas->ridx, atlas->nrh, atlas->meta, (float)Params.TimeStep);
            sprintf(fpath, "%sMean.lh.gii", Name2);
            PutGIFTISurf(fpath, ImgM, NumImg, atlas->nrh, atlas->lidx, atlas->nlh, atlas->meta, (float)Params.TimeStep);
        } else {
            sprintf(fpath, "%sMean.nii", Name2);
            fp = fileopen(fpath, "w");

            // modify NIFTI header
            memset(extension, 0, 4);                                    // for now, ignore the fact that we may have an extended header
            NiiHdr.dim[0] = 4;                                          // 3d+time
            NiiHdr.dim[4] = NumImg;                                     // number of 3d images
            NiiHdr.toffset = (float)Params.Marker[0].MarkStart;         // time offset shift
            NiiHdr.slice_duration = (float)Params.TimeStep;             // time step size --- unused?
            NiiHdr.pixdim[4] = (float)Params.TimeStep;                  // TR is really stored here
            NiiHdr.slice_code = NIFTI_SLICE_SEQ_INC;
            NiiHdr.xyzt_units = (NIFTI_UNITS_MM | NIFTI_UNITS_SEC);     // 3d+time units
            memset(NiiHdr.intent_name, 0, 16);
            memcpy(NiiHdr.intent_name, "Mean", 4);

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
            for (d=0; d<NumImg; d++)    // write time as BRIKs
                if (fwrite((void *)ImgM[d], sizeof(float), V, fp) != V)
                    Cleanup("can't write data image file");
            fclose(fp);
        }

        // write variance image
        if (vflg) {
            printf(" - done\n");
            printf("writing 'Variance'");
            fflush(stdout);
        }
        if (gflg) {
            /* Right hemisphere is stored first in ImgV! */

            sprintf(fpath, "%sVariance.rh.gii", Name2);
            PutGIFTISurf(fpath, ImgV, NumImg, 0, atlas->ridx, atlas->nrh, atlas->meta, (float)Params.TimeStep);
            sprintf(fpath, "%sVariance.lh.gii", Name2);
            PutGIFTISurf(fpath, ImgV, NumImg, atlas->nrh, atlas->lidx, atlas->nlh, atlas->meta, (float)Params.TimeStep);
        } else {
            sprintf(fpath, "%sVariance.nii", Name2);
            fp = fileopen(fpath, "w");

            // modify NIFTI header
            memset(extension, 0, 4);                                    // for now, ignore the fact that we may have an extended header
            NiiHdr.dim[0] = 4;                                          // 3d+time
            NiiHdr.dim[4] = NumImg;                                     // number of 3d images
            NiiHdr.toffset = (float)Params.Marker[0].MarkStart;         // time offset shift
            NiiHdr.slice_duration = (float)Params.TimeStep;             // time step size --- unused?
            NiiHdr.pixdim[4] = (float)Params.TimeStep;                  // TR is really stored here
            NiiHdr.slice_code = NIFTI_SLICE_SEQ_INC;
            NiiHdr.xyzt_units = (NIFTI_UNITS_MM | NIFTI_UNITS_SEC);     // 3d+time units
            memset(NiiHdr.intent_name, 0, 16);
            memcpy(NiiHdr.intent_name, "Variance", 8);

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
            for (d=0; d<NumImg; d++)    // write time as BRIKs
                if (fwrite((void *)ImgV[d], sizeof(float), V, fp) != V)
                    Cleanup("can't write data image file");
            fclose(fp);
        }
    }                   // for(n...

    log_params(SAMpath);

    if (vflg) {
        msg(" - done\n'%s' done\n", Progname);
    }
    if (nflg) {
        gsl_matrix_free(Cn);
    }
    exit(0);
}
