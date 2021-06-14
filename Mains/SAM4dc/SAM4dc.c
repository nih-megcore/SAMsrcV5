// SAM4dc -- synthetic aperture magnetometry imaging for source event related
//  time-locked signals. SAM4dc reads all SAM weights for volume image, & computes
//  the source waveform for each voxel. Given specifications for active & control
//  time segments, SAM4dc computes the mean power & variance time-series relative
//  to the specified event windows. The result is represented stored in a NIFTI
//  image file '.nii' file. It is presumed that the marker start and end times
//  are identical!
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
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <fcntl.h>
#include <gsl/gsl_matrix.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <rvelib.h>
#include <Params.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <samUtils.h>
#include <filters.h>
#include <TimeSegs.h>
#include <nifti1.h>
#include <version.h>

// internal constants
#define MINOR_REV   0
#define TRUE        1
#define FALSE       0
#define SIG_ORDER   4

typedef struct {
    int     Epoch;                  // epoch/trial
    int     TS;                     // starting sample
} SAM_SEG;

typedef struct {                    // image statistics
    char        Name[64];           // weight file name
    gsl_matrix  *Wgt;               // Wgt - VxM SAM weights
    double      Noise;              // noise power
    SAM_HDR     WgtHdr;             // SAM weight file header
    SAM_SEG     *Segment;           // Segment[TotSegments] -- dataset integration time-segment specs
    int         NumSegments;        // number of segments for this dataset
} IMAGE_STAT;


int main(
    int     argc,
    char    **argv)
{
    IIRSPEC         DataIIR;        // IIR data filter structure
    IIRSPEC         SmoothIIR;      // IIR smoothing filter structure
    FFTSPEC         DataFFT;        // FFT filter structure
    FFTSPEC         SmoothFFT;      // FFT filter for smoothing metric envelope
    FFTSPEC         NotchFFT;       // FFT notch filter
    IMAGE_STAT      *Stats;         // Stats[N] -- structure for conditions
    PARMINFO        Params;         // analysis parameters
    SAM_MARKS       *Marker;        // Marker[N] -- array of time-ordered markers
    HeaderInfo      Header;         // MEG data header
    ChannelInfo     *Channel;       // MEG channel info
    EpochInfo       *Epoch;         // MEG epoch info
    COV_HDR         CovHdr;         // noise covariance header
    nifti_1_header  NiiHdr;         // NIFTI header
    unsigned char   *ExtHdr;        // NiiHdr[N] -- extended NIFTI header
    gsl_matrix      *Cn;            // Cn - MxM noise covariance matrix
    gsl_vector      *W;             // W - Mx1 weight vector
    double          ***x;           // x[E][T][M] -- filtered multisensor data
    double          *avg;           // avg[TA] -- average (before downsampling)
    double          *var;           // var[TA] -- variance (before downsampling)
    double          *In;            // In[TT] -- input time-series
    double          *Out;           // Out[TT] -- output time-series
    double          mark;           // time mark
    double          time;           // time (seconds)
    double          BWFreq;         // effective bandwidth
    double          BWratio;        // nominal ratio of image bandwidth to noise bandwidth
    double          Noise;          // effective noise
    double          percent;        // percent done
    double          a1;             // sum of metric
    double          a2;             // sum of metric squared
    double          baseline;       // baseline of average
    double          tb;             // number of samples in baseline
    double          ns;             // number of segments
    float           **ImgM;         // ImgM[D][V] -- image mean
    float           **ImgV;         // ImgV[D][V] -- image variance
    register double tmp;            // what else?
    int             *ChanIndex;     // ChanIndex[M] -- index of channels used
    int             E;              // number of good epochs
    int             M;              // number of MEG or reference channels used
    int             N;              // just a number to take care of extended header read
    int             T;              // number of samples
    int             TT;             // E*T
    int             TA;             // number of samples in average window
    int             TB;             // number of samples in baseline window
    int             TS;             // average starting sample index
    int             TE;             // average ending sample index
    int             BS;             // baseline starting sample relative to average start
    int             BE;             // baseline ending sample relative to average start
    int             V;              // number of voxels weight file
    int             c;              // input character & virtual channel
    int             pct;            // integer percentage
    int             m;              // channel index
    int             d;              // image step index
    int             e;              // epoch index
    int             f;              // 1st primary sensor index
    int             t;              // sample index
    int             tt;             // alternate sample index
    int             i;              // matrix index
    int             j;              // ...and another
    int             k;              // this is the last one -- for sure...
    int             n;              // condition index
    int             NumTot;         // number of total segments
    int             NumMarkers;     // number of markers in 'MarkerFile.mrk'
    int             NumImg;         // number of SAM images
    int             LongIndex;      // position of arguments
    static int      bflg = FALSE;   // baseline flag
    static int      eflg = FALSE;   // command-line error flag
    static int      mflg = FALSE;   // marker parameter file name flag
    static int      nflg = FALSE;   // use noise covariance file
    static int      rflg = FALSE;   // run name flag
    static int      vflg = FALSE;   // verbose mode flag
    extern char     *optarg;
    extern int      opterr;
    char            fpath[256];     // general path name
    char            DSName[256];    // MEG dataset name
    char            DSpath[256];    // MEG dataset path
    char            SAMpath[256];   // SAM subdirectory path
    char            WgtDir[256];    // weight directory path name
    char            ParmPath[256];  // full path of parameter file
    char            ParmName[256];  // name of marker parameter file
    char            Prefix[64];     // incomplete pathname
    char            Name1[256];     // partial image name
    char            Name2[256];     // partial image name
    char            FilterName[8];  // filter name -- IIR | FFT
    char            extension[4];   // extension
    FILE            *fp;            // data output file pointer
    FILE            *np;            // noise file pointer

    //system-dependent variables
#if BTI
    static int      dflg = FALSE;           // pdf file name flag
    char            PDFName[256];           // name of data file
    char            ShortOpts[] = "r:d:m:v";
    static struct option    LongOpts[] = {  // command line options
        {"dataset", required_argument, NULL, 'r'},
        {"pdf_name", required_argument, NULL, 'd'},
        {"analysis_name", required_argument, NULL, 'm'},
        {"verbose", no_argument, NULL, 'v'},
        {NULL, no_argument, NULL, 0}
    };
#else
    double          scale;                  // bits-->Tesla conversion
    unsigned char   **Bad;                  // Bad[E][T] -- flags for bad samples
    char            ShortOpts[] = "r:m:v";
    static struct option    LongOpts[] = {  // command line options
        {"dataset", required_argument, NULL, 'r'},
        {"analysis_name", required_argument, NULL, 'm'},
        {"verbose", no_argument, NULL, 'v'},
        {NULL, no_argument, NULL, 0}
    };
#endif

    // parse command line parameters
    opterr = 1;             // enable error reporting
    while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
        switch (c) {
            case 'r':       // datset name
                sprintf(DSName, "%s", optarg);
                rflg = TRUE;
                break;
#if BTI
            case 'd':       // data file name
                sprintf(PDFName, "%s", optarg);
                dflg = TRUE;
                break;
#endif
            case 'm':       // parameter file name
                sprintf(ParmPath, "%s", optarg);
                mflg = TRUE;
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
        fprintf(stderr, "SAM4dc [arguments]\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
        fprintf(stderr, "\t-r <dataset name>\n");
        fprintf(stderr, "\t-d <data file name>\n");
#else
    if(eflg == TRUE || rflg == FALSE || mflg == FALSE) {
        fprintf(stderr, "SAM4dc [arguments]\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
        fprintf(stderr, "\t-r <dataset name>\n");
#endif
        fprintf(stderr, "\t-m <parameter file name>\n");
        fprintf(stderr, "\t-v -- verbose mode\n");
        exit(-1);
    }

    // get dataset information structures
    if(vflg == TRUE) {
        printf("4 D   S A M   I M A G I N G\t%s\n", PRG_REV);
        printf("opening data file");
        fflush(stdout);
    }

    // get data with sensor structures
#if BTI
    if(GetMagnesInfo(DSName, PDFName, &Header, &Channel, &Epoch) == -1)
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
    TT = E * T;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("'%s' uses %d epochs, %d channels, %d samples/channel/epoch\n", DSName, E, M, T);
        fflush(stdout);
    }

    // sanity check
    if(M == 0 || E == 0 || T == 0)
        Cleanup("no valid data");

    // parse analysis parameter specification file in current working directory
    ParseFileName(ParmPath, ParmName);
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("parsing '%s' parameter file", ParmName);
        fflush(stdout);
    }
    sprintf(fpath, "%s.param", ParmPath);
    GetParams(fpath, &Params);

    // announce metric type
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("Imaging functions of ");
        switch(Params.ImageMetric) {
            case MOMENT: printf(" moment\n");break;
            case POWER: printf("power\n"); break;
            case RV_ENTROPY: printf("rank vector entropy for %d dimensions\n", Params.Dims); break;
            case -1: Cleanup("parameters must specify 'ImageMetric'"); break;
            default: Cleanup("selected ImageMetric not supported"); break;
        }
        fflush(stdout);
    }

    // check imaging parameters
    if(Params.BaseStart != -999. && Params.BaseEnd != -999.)
        bflg = TRUE;
    if(Params.NumMark == -1)
        Cleanup("'NumMarkers' specification req'd");
    if(Params.CovType == -1)
        Cleanup("'CovType' is mandatory");
    if(Params.CovHP == -999. || Params.CovLP == -999.)          // CovBand is always required!
        Cleanup("'CovBand' specifications required");
    if(Params.ImageHP == -999. || Params.ImageLP == -999.)      // ImageBand is always required!
        Cleanup("'ImageBand' specifications required");
    if(Params.ImageHP < Channel[f].AnalogHpFreq)
        Cleanup("ImageBand highpass cutoff frequency outside recorded range");
    if(Params.ImageLP > Channel[f].AnalogLpFreq)
        Cleanup("ImageBand lowpass cutoff frequency outside recorded range");
    if(Params.SmoothHP == -999 || Params.SmoothLP == -999.)
        Cleanup("'SmoothBand' specifications required");
    if(Params.FilterType == IIR)
        strcpy(FilterName, "IIR");
    else
        strcpy(FilterName, "FFT");
    if(Params.TimeStep < 0.)
        Cleanup("must specify 'TimeStep'");
    if(Params.ImageMetric == RV_ENTROPY)
        if(Params.Dims < 3 || Params.Dims > 6)
            Cleanup("RVE dimensions must be 3 to 6");
    BWratio = (Params.ImageLP - Params.ImageHP) / (Params.CovLP - Params.CovHP);

    // set up 3d+time imaging constants
    time = Params.Marker[0].MarkEnd - Params.Marker[0].MarkStart;       // trial window duration (seconds)
    NumImg = (int)rint(time / Params.TimeStep);                         // number of 3d images

    // set prefix to NumPrefix characters (default is 8-character dataset hashcode)
    memset(Prefix, 0, 64);
#if BTI
    memcpy(Prefix, "BIOMAG4D", 8);              // 4D datasets just have serial numbers
#else
    memcpy(Prefix, DSName, Params.NumPrefix);   // default: copy eight-letter NIMH hash code in dataset name
#endif

    // use CovBand frequencies to select subdirectory
    if(vflg == TRUE) {
        printf("selecting filter parameters");
        fflush(stdout);
    }

    // if req'd, initialize NOISE
    if(Params.NoiseHP != -999. && Params.NoiseLP != -999.) {
        if(vflg == TRUE) {
            printf("reading noise covariance file");
            fflush(stdout);
        }
        sprintf(fpath, "%s/%s,%-d-%-dHz/Noise.cov", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
        Cn = gsl_matrix_alloc(M, M);
        GetCov(fpath, &CovHdr, &ChanIndex, Cn);
        nflg = TRUE;
    }

    // if req'd, summarize epoch dimensions
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("dataset '%s' has %d epochs\n", DSName, E);
        printf("\tduration %6.3f seconds/epoch\n", (double)T / Header.SampleRate);
        fflush(stdout);
    }

    // read markers
#if BTI
    if(vflg == TRUE) {
        printf("parsing trigger markers");
        fflush(stdout);
    }
    if(GetMagnesMarkers(DSpath, PDFName, &Header, Channel, &Marker, &NumMarkers) == -1)
        Cleanup("cannot get trigger markers");
#else
    if(vflg == TRUE) {
        printf("reading 'MarkerFile.mrk'");
        fflush(stdout);
    }
    sprintf(fpath, "%s/MarkerFile.mrk", DSpath);
    GetMarkers(fpath, &Marker, &NumMarkers);
#endif
    if(vflg == TRUE) {
        printf(" - found %d markers", NumMarkers);
        printf(" - done\n");
        fflush(stdout);
    }

    // compute limits
    time = Params.Marker[0].MarkEnd - Params.Marker[0].MarkStart;       // duration of average window (for all markers!)
    TA = (int)rint(time * Header.SampleRate);                           // number of samples in average window
    TS = (int)rint(Params.Marker[0].MarkStart * Header.SampleRate);     // starting sample, relative to marker
    TE = (int)rint(Params.Marker[0].MarkEnd * Header.SampleRate);       // ending sample, relative to marker
    if(bflg == TRUE) {
        time = Params.BaseEnd - Params.BaseStart;                       // duration of baseline window (seconds)
        TB = (int)rint(time * Header.SampleRate);                       // number of samples in baseline window
        tb = (double)TB;
        BS = (int)rint(Params.BaseStart * Header.SampleRate);           // baseline starting sample relative to trigger marker (can be negative!)
        BE = (int)rint(Params.BaseEnd * Header.SampleRate);             // baseline ending sample relative to trigger marker (can be negative!)
        if(BS < TS || BS > TE || BE < TS || BE > TE)
            Cleanup("baseline starting & ending samples out of range");
    }
    BS -= TS;                                                           // adjust baseline relative to start of average window

    // test for identical segment start & end times
    if(Params.NumMark > 1) {
        for(n=1; n<Params.NumMark; n++)
            if(Params.Marker[n].MarkStart != Params.Marker[0].MarkStart || Params.Marker[n].MarkEnd != Params.Marker[0].MarkEnd)
                Cleanup("marker time windows must be identical");
    }

    // populate segment array
    if(vflg == TRUE) {
        printf("populating segments:\n");
        fflush(stdout);
    }

    // allocate 'Stats' structure
    if((Stats = (IMAGE_STAT *)malloc((size_t)Params.NumMark * sizeof(IMAGE_STAT))) == NULL)
        allocfailed("Stats[]");

    // for each marker...
    for(n=0; n<Params.NumMark; n++) {

        // carry over marker name
        sprintf(Stats[n].Name, "%s", Params.Marker[n].MarkName);
        if(vflg == TRUE) {
            printf("\t'%s'\n", Stats[n].Name);
            fflush(stdout);
        }

        // count the total markers for each named marker
        for(i=NumTot=0; i<NumMarkers; i++)
            if(!strcmp(Params.Marker[n].MarkName, Marker[i].Name))
                NumTot++;
        if(NumTot == 0)
            Cleanup("no markers found with specified names");
        if(NumTot < 3)
            Cleanup("can't compute average or variance with < 3 markers");

        // allocate segment array
        if((Stats[n].Segment = (SAM_SEG *)malloc((size_t)NumTot * sizeof(SAM_SEG))) == NULL)
            allocfailed("Stats[].Segment[]");

        // populate Stats for each condition
        for(i=k=Stats[n].NumSegments=0, e=(-1); i<NumMarkers; i++) {
            if(Marker[i].Epoch != e)
                e = Marker[i].Epoch;
            if(Epoch[e].Flag == TRUE) {
                
                // populate segment windows
                if(!strcmp(Params.Marker[n].MarkName, Marker[i].Name)) {
                    mark = Marker[i].Latency;
                    Stats[n].Segment[k].Epoch = e;
                    TS = (int)rint((mark + Params.Marker[n].MarkStart) * Header.SampleRate) + Header.PreTrigSamples;
                    TE = (int)rint((mark + Params.Marker[n].MarkEnd) * Header.SampleRate) + Header.PreTrigSamples;
                    if(TS >= 0 && TE < T) {
                        Stats[n].Segment[k].TS = TS;
                        Stats[n].NumSegments++;
                        k++;
                    }
                }       // if(!strcmp...
            }           // if(Epoch[e].Flag...
        }               // for(i...
        if(Stats[n].NumSegments == 0)
            Cleanup("there were no valid data segments for specified markers & time windows");
    }                   // for(n...

    // read weights according to CovType
    if(vflg == TRUE) {
        printf("reading weights");
        fflush(stdout);
    }
    sprintf(WgtDir, "%s/%s,%-d-%-dHz", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
    switch(Params.CovType) {

        case GLOBAL_:       // read Global weights & Global_Noise
            sprintf(fpath, "%s/Global_Noise", WgtDir);
            if((np = fopen(fpath, "r")) == NULL)
                Cleanup("can't open noise file");
            if(fscanf(np, "%le", &Stats[GLOBAL_].Noise) != 1)
                Cleanup("can't read noise file");
            fclose(np);
            if(Params.ImageFormat == TLRC)
                sprintf(fpath, "%s/Global_at.nii", WgtDir);
            else
                sprintf(fpath, "%s/Global.nii", WgtDir);
            for(n=0; n<Params.NumMark; n++) {
                Stats[n].Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
                V = Stats[n].Wgt->size1;
                if(Stats[n].Wgt->size2 != M)
                    Cleanup("number of channels in dataset does not match number in weight file");
                Stats[n].Noise = Stats[GLOBAL_].Noise;
            }
            break;

        case SUM_:          // read Sum weights
            sprintf(fpath, "%s/Sum_Noise", WgtDir);
            if((np = fopen(fpath, "r")) == NULL)
                Cleanup("can't open noise file");
            if(fscanf(np, "%le", &Stats[SUM_].Noise) != 1)
                Cleanup("can't read noise file");
            fclose(np);
            if(Params.ImageFormat == TLRC)
                sprintf(fpath, "%s/Sum_at.nii", WgtDir);
            else
                sprintf(fpath, "%s/Sum.nii", WgtDir);
            for(n=0; n<Params.NumMark; n++) {
                Stats[n].Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
                V = Stats[n].Wgt->size1;
                if(Stats[n].Wgt->size2 != M)
                    Cleanup("number of channels in dataset does not match number in weight file");
                Stats[n].Noise = Stats[SUM_].Noise;
            }
            break;

        case ALL_:          // read weights for each individual marker
            for(n=0; n<Params.NumMark; n++) {
                sprintf(fpath, "%s/%s_Noise", WgtDir, Stats[n].Name);
                if((np = fopen(fpath, "r")) == NULL)
                    Cleanup("can't open noise file");
                if(fscanf(np, "%le", &Stats[n].Noise) != 1)
                    Cleanup("can't read noise file");
                fclose(np);
                if(Params.ImageFormat == TLRC)
                    sprintf(fpath, "%s/%s_at.nii", WgtDir, Stats[n].Name);
                else
                    sprintf(fpath, "%s/%s.nii", WgtDir, Stats[n].Name);
                Stats[n].Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
                V = Stats[n].Wgt->size1;
                if(Stats[n].Wgt->size2 != M)
                    Cleanup("number of channels in dataset does not match number in weight file");
            }
            break;

        default:
            Cleanup("unknown 'CovType'");
            break;
    }

    // allocate memory
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("allocating memory");
        fflush(stdout);
    };
    W = gsl_vector_alloc(M);
    if((ImgM = (float **)malloc((size_t)NumImg * sizeof(float *))) == NULL)
        allocfailed("ImgM[]");
    for(d=0; d<NumImg; d++)
        if((ImgM[d] = (float *)malloc((size_t)V * sizeof(float))) == NULL)
            allocfailed("ImgM[][]");
    if((ImgV = (float **)malloc((size_t)NumImg * sizeof(float *))) == NULL)
        allocfailed("ImgV[]");
    for(d=0; d<NumImg; d++)
        if((ImgV[d] = (float *)malloc((size_t)V * sizeof(float))) == NULL)
            allocfailed("ImgV[][]");
    if((x = (double ***)malloc((size_t)E * sizeof(double **))) == NULL)
        allocfailed("x[]");
    for(e=0; e<E; e++) {
        if((x[e] = (double **)malloc((size_t)T * sizeof(double *))) == NULL)
            allocfailed("x[][]");
        for(t=0; t<T; t++)
            if((x[e][t] = (double *)malloc((size_t)M * sizeof(double))) == NULL)
                allocfailed("x[][][]");
    }
    if((avg = (double *)malloc((size_t)TA * sizeof(double))) == NULL)
        allocfailed("avg[]");
    if((var = (double *)malloc((size_t)TA * sizeof(double))) == NULL)
        allocfailed("std[]");
    if((In = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
        allocfailed("In[]");
    if((Out = (double *)malloc((size_t)TT * sizeof(double))) == NULL)
        allocfailed("Out[]");

    // determine filter type
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("making %6.3f to %6.3f Hz %s filter for weight bandpass", Params.ImageHP, Params.ImageLP, FilterName);
        fflush(stdout);
    }
    if(Params.FilterType == IIR)
        mkfilter(&DataIIR, SIG_ORDER, Header.SampleRate, Params.ImageHP, Params.ImageLP, &BWFreq);
    else
        Butterworth(&DataFFT, T, Params.ImageHP, Params.ImageLP, Header.SampleRate, BANDPASS, MAXORDER, &BWFreq);

    // if req'd, generate IIR smoothing filter
    if(Params.ImageMetric != RV_ENTROPY) {
        if(vflg == TRUE) {
            printf(" - done\n");
            printf("making %6.3f Hz %s lowpass filter for smoothing metric", Params.SmoothLP, FilterName);
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
        if(vflg == TRUE) {                              // if req'd - heartbeat
            printf(".");
            fflush(stdout);
        }
    }

    // for each condition -- only 1 for GLOBAL & SUM, but multiple for ALL
    for(n=0; n<Params.NumMark; n++) {

        // announce condition
        if(vflg == TRUE) {
            printf(" - done\n");
            printf("computing image voxels for '%s':\n", Stats[n].Name);
            fflush(stdout);
        }

        // loop on all coordinates in the ROI
        for(c=0, pct=(-1); c<V; c++) {

            // get voxel weight vector & determine if within head boundary
            gsl_matrix_get_row(W, Stats[n].Wgt, c);
            if(gsl_vector_isnull(W) == 0) {     // compute image only for non-zero weights

                // compute rms noise
                if(nflg == FALSE) {
                    for(m=0, tmp=0.; m<M; m++)
                        tmp += gsl_vector_get(W, m) * Stats[n].Noise * gsl_vector_get(W, m);
                    tmp *= BWratio;
                } else {
                    for(i=0, tmp=0.; i<M; i++)
                        for(j=0; j<M; j++)
                            tmp += gsl_vector_get(W, i) * gsl_matrix_get(Cn, i, j) * gsl_vector_get(W, j);
                    tmp *= (BWFreq / CovHdr.BWFreq);    // adjust noise power to image bandwidth
                }
                Noise = sqrt(tmp);

                // compute voxel time-series for entire dataset
                for(e=tt=0; e<E; e++)
                    for(t=0; t<T; t++, tt++) {
                        for(m=0, tmp=0.; m<M; m++)
                            tmp += x[e][t][m] * gsl_vector_get(W, m);
                        In[tt] = tmp / Noise;
                    }

                // compute metric for entire dataset
                switch(Params.ImageMetric) {
                    case MOMENT:                // event-related moment
                        for(t=0; t<TT; t++)
                            Out[t] = In[t];
                        break;
                    case POWER:                 // Hilbert envelope of power
                        hilbert(TT, In, Out);
                        for(t=0; t<TT; t++)     // sum the in-phase & quadrature signal components
                            Out[t] = 0.5 * (In[t] * In[t] + Out[t] * Out[t]);
                        break;
                    case RV_ENTROPY:            // rank vector entropy
                        Data2RVE(In, Out, TT, Params.Dims, Params.ImageLP, Params.Tau, Header.SampleRate);
                        break;
                    default:
                        Cleanup("ImageMetric not implemented");
                        break;
                }

                // smooth metric time-series with a lowpass filter
                if(Params.ImageMetric != RV_ENTROPY) {
                    if(Params.FilterType == IIR)
                        bdiir(Out, In, TT, &SmoothIIR);     // lowpass to smooth metric
                    else
                        FFTfilter(Out, In, &SmoothFFT, TT); // lowpass to smooth metric
                } else {
                    for(t=0; t<TT; t++)
                        In[t] = 100. * Out[t];
                }

                // compute average (before downsampling)
                ns = (double)Stats[n].NumSegments;
                for(t=0; t<TA; t++) {
                    for(k=0, a1=a2=0.; k<Stats[n].NumSegments; k++) {
                        e = Stats[n].Segment[k].Epoch;
                        tt = t + e * T + Stats[n].Segment[k].TS;
                        a1 += In[tt];
                        a2 += In[tt] * In[tt];
                    }
                    avg[t] = a1 / ns;
                    var[t] = (a2 + a1 * a1 / ns) / (ns - 1.);
                }

                // if req'd, remove baseline
                if(bflg == TRUE) {
                    for(t=0, tt=BS, tmp=0.; t<TB; t++, tt++)
                        tmp += avg[tt];
                    baseline = tmp / tb;
                    for(t=0; t<TA; t++)
                        avg[t] -= baseline;
                }

                // downsample averge & variance
                for(d=0, time=0.; d<NumImg; d++, time+=Params.TimeStep) {
                    t = (int)rint(time * Header.SampleRate);
                    ImgM[d][c] = (float)avg[t];
                    ImgV[d][c] = (float)var[t];
                }
            } else {    // if(tmp...
                for(d=0; d<NumImg; d++) {
                    ImgM[d][c] = 0.;
                    ImgV[d][c] = 0.;
                }
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
        }           // for(c...

        // write mean & variance images
        if(vflg == TRUE) {
            printf("\nwriting 'Mean'");
            fflush(stdout);
        }
        if(!strncmp(Params.DirName, "NULL", 4))
            sprintf(Name1, "%s/%s,%s,%s,", SAMpath, Prefix, ParmName, Stats[n].Name);
        else
            sprintf(Name1, "%s/%s,%s,%s,", Params.DirName, Prefix, ParmName, Stats[n].Name);
        switch(Params.ImageMetric) {
            case MOMENT:        sprintf(Name2, "%s4DC_MOM,", Name1); break;
            case POWER:         sprintf(Name2, "%s4DC_PWR,", Name1); break;
            case RV_ENTROPY:    sprintf(Name2, "%s4DC_RVE,", Name1); break;
            default:            Cleanup("ImageMetric not implemented"); break;
        }
        sprintf(fpath, "%sMean.nii", Name2);
        if((fp = fopen(fpath, "w")) == NULL)
            Cleanup("can't open .nii file for write");

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
        if(fwrite((void *)&NiiHdr, sizeof(nifti_1_header), 1, fp) != 1)
            Cleanup("can't write nifti header to .nii file");
        if(fwrite((void *)extension, 4, 1, fp) != 1)
            Cleanup("can't write 'extension' to .nii file");
        if((int)NiiHdr.vox_offset != 352 && ExtHdr != NULL) {
            N = (int)NiiHdr.vox_offset - 352;
            if(fwrite((void *)ExtHdr, N, 1, fp) != 1)
                Cleanup("can't write extended header to .nii file");
        }
        for(d=0; d<NumImg; d++)     // write time as BRIKs
            if(fwrite((void *)ImgM[d], sizeof(float), V, fp) != V)
                Cleanup("can't write data image file");
        fclose(fp);

        // write variance image
        if(vflg == TRUE) {
            printf(" - done\n");
            printf("writing 'Variance'");
            fflush(stdout);
        }
        sprintf(fpath, "%sVariance.nii", Name2);
        if((fp = fopen(fpath, "w")) == NULL)
            Cleanup("can't open .nii file for write");

        // modify NIFTI header
        memset(extension, 0, 4);                                    // for now, ignore the fact that we may have an extended header
        NiiHdr.dim[0] = 4;                                          // 3d+time
        NiiHdr.dim[4] = NumImg;                                     // number of 3d images
        NiiHdr.toffset = (float)Params.Marker[0].MarkStart;         // time offset shift
        NiiHdr.slice_duration = (float)Params.TimeStep;             // time step size
        NiiHdr.pixdim[4] = (float)Params.TimeStep;                  // TR is really stored here
        NiiHdr.slice_code = NIFTI_SLICE_SEQ_INC;
        NiiHdr.xyzt_units = (NIFTI_UNITS_MM | NIFTI_UNITS_SEC);     // 3d+time units
        memset(NiiHdr.intent_name, 0, 16);
        memcpy(NiiHdr.intent_name, "Variance", 8);

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
        for(d=0; d<NumImg; d++)     // write time as BRIKs
            if(fwrite((void *)ImgV[d], sizeof(float), V, fp) != V)
                Cleanup("can't write data image file");
        fclose(fp);
    }               // for(n...

    gsl_matrix_free(Cn);
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("'SAM4dc' done\n");
        fflush(stdout);
    }
    exit(0);
}
