// SAMcov -- compute one or more covariance files
//
//  An analysis parameter file (*.param) is parsed to obtain the time window(s)
//  frequency passband, & markers fir the covariance analysis,
//  Also, allow for generation of a covariance for estimating noise (as before).
//
//  flags:
//      -r <run name>       -- MEG dataset name
//      -d <data file name> -- MEG data file name (4D, only)
//      -m                  -- read parameter specifications
//      -v                  -- verbose mode
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
#include <sys/time.h>
#include <errno.h>
#include <fcntl.h>
#include <ctype.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <siglib.h>
#include <samlib.h>
#include <filters.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <TimeSegs.h>
#include <version.h>
#include "samutil.h"
#include "sam_parse.h"

#define MINOR_REV   5
#define SIG_ORDER   4

typedef struct {
    double  **Cov;                  // there is a covariance matrix for each marker/segment definition
    double  Noise;                  // effective noise
    COV_SEG *Segment;               // Segment[TotSegments] -- covariance integration time-segment specs
    int     NumSegments;            // number of segments for this covariance file
    int     NumSamples;             // number of samples for this covariance file
    char    CovFileName[256];       // name of this covariance file
    int     Used;                   // flag indicating covariance in use
} COV_STATS;

int main(
    int     argc,
    char    **argv)
{
    IIRSPEC         CovIIR;         // IIR filter structure for computing data covariance
    IIRSPEC         OrientIIR;      // IIR filter structure for computing orientation covariance
    IIRSPEC         NoiseIIR;       // IIR filter structure for computing noise covariance
    FFTSPEC         CovFFT;         // FFT filter structure for computing data covariance
    FFTSPEC         OrientFFT;      // FFT filter structure for computing orientation covariance
    FFTSPEC         NoiseFFT;       // FFT filter structure for computing noise covariance
    FFTSPEC         NotchFFT;       // FFT notch filter structure
    COV_STATS       *Stats;         // separate statistics for each covariance file
    PARMINFO        Params;         // analysis parameters
    PARM            *p;             // pointer into active parameter list
    SAM_MARKS       *Marker;        // Marker[N] -- array of time-ordered markers
    HeaderInfo      Header;         // header information data
    ChannelInfo     *Channel;       // Channel[M] -- channel information data
    EpochInfo       *Epoch;         // Epoch[E] -- epoch information data
    COV_HDR         CovHeader;      // SAM covariance file header information
    gsl_matrix      *U;             // SVD matrix U
    gsl_matrix      *V;             // SVD matrix V
    gsl_vector      *S;             // singular values
    double          **A;            // A[M][M] -- covariance matrix of segment
    double          ***x;           // x[E][M][T] -- signal MEG data
    double          **seg;          // seg[M][T] -- segment time series (long enough for anything)
    double          *In;            // In[T] -- input time-series
    double          *Out;           // Out[T] -- output time-series
    double          *d;             // tmp pointer
    double          mark;           // time of marker
    double          CovBW;          // bandwidth of image filter
    double          OrientBW;       // bandwidth of "orientation" filter
    double          NoiseBW;        // bandwidth of "noise" filter
    double          BW;             // temporary bandwidth
    double          mean;           // mean value of segment
    double          Tolerance;      // span of singular values
    double          LSNoise;        // noise of least-significant singular value
    double          EffSamples;     // effective number of samples (normalized by number of primary sensors)
    double          min;            // minimum delta;
    double          delta;          // derivative
    double          time;           // the two-timer variable
    double          ttime;          // ...yup!
    double          tmp;            // integrator
    int             ndim = 3;       // number of skipped eigenvalues when determining minimum for noise
    int             f;              // 1st primary sensor channel index
    int             i;              // generic matrix index
    int             j;              // another generic matrix index
    int             k;              // marker index
    int             t;              // time index
    int             tt;             // alternate time index
    int             m;              // sensor index
    int             n;              // alternate sensor index
    int             e;              // epoch index
    int             E;              // even epoch count
    int             M;              // number of sensor channels used
    int             T;              // maximum number of samples in epoch
    int             TT;             // total samples
    int             NumCovs;        // number of covariance matrices used
    int             NumMarkers;     // number of markers in 'MarkerFile.mrk'
    int             TotSegments;    // number of time segments in final specs
    int             SumSegments;    // number of time segments in sum covariance
    static int      eflg = FALSE;   // error flag
    static int      mflg = FALSE;   // parameter file specification flag
    static int      nflg = FALSE;   // noise covariance filter flag
    static int      oflg = FALSE;   // orientation filter flag
    static int      pflg = FALSE;   // use pseudo-inverse and skip lowest eigenvalues of covariance
    static int      rflg = FALSE;   // dataset name flag
    static int      vflg = FALSE;   // verbose mode flag
    char            *DSName = NULL; // dataset name
    char            DSpath[256];    // full dataset path
    char            SAMpath[256];   // expanded SAM path
    char            COVpath[256];   // covariance file path
    char            fpath[256];     // general path name
    char            *CovName = NULL; // name for output files
    char            FilterName[8];  // filter name (IIR | FFT)
    char            *s;             // tmp char pointer
    FILE            *fp;            // summary file

    // system-dependent variables
#if BTI
    static int      dflg = FALSE;   // pdf file name flag
    char            *PDFName;       // name of data file
#else
    double          scale;          // bits-->Tesla conversion
    unsigned char   **Bad;          // Bad[E][T] -- flags for bad samples
#endif

    // Register the parameters used by this program.

    reg_usage(argv[0], MINOR_REV, __DATE__);
    reg_std_parm();

    reg_parm("Marker");
    reg_parm("SegFile");

    reg_parm("CovBand");
    reg_parm("ImageBand");
    reg_parm("NoiseBand");
    reg_parm("OrientBand");
    reg_parm("FilterType");
    reg_parm("DataSegment");
    reg_parm("Notch");
    reg_parm("Hz");
    reg_parm("Pinv");

    reg_parm("OutName");
    reg_parm("CovName");
    set_parm_help("CovName", "same as --OutName");

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
    p = get_parm("Pinv");
    pflg = p->set;
    if (pflg) {
        ndim = *(int *)p->ptr;
    }

#if BTI
    if (!rflg || !dflg || !mflg || eflg) {
#else
    if (!rflg || !mflg || eflg) {
#endif
        msg("dataset (-r) and parameter file (-m) are required\n");
        do_help();  // doesn't return
    }

    CovName = Params.ParmName;
    p = get_parm("OutName");
    if (p->set) {
        CovName = (char *)p->ptr;
    }
    p = get_parm("CovName");            // -C overrides -N if both are used ...
    if (p->set) {
        CovName = (char *)p->ptr;
    }

    // announce
    if (vflg) {
        printf("C O V A R I A N C E   F I L E   G E N E R A T O R\t%s\n", PRG_REV);
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
    GetFilePath(DSDIR, SAMpath, sizeof(SAMpath), &Params, "SAM", 0);

    // create SAM subdirectory
    if (mkdir(SAMpath, S_IRWXU | S_IRWXG | S_IRWXO) == -1)
        if (errno != EEXIST)
            Cleanup("can't create 'SAM' subdirectory");

    // count data dimensions
    E = Header.NumEpochs;
    T = Header.MaxSamples;
    TT = E * T;
    time = (double)T / Header.SampleRate;       // time/epoch
    ttime = (double)TT / Header.SampleRate;     // total time
    M = Header.NumPri;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index
    if (vflg) {
        printf(" - done\n");
        printf("'%s' uses %d epochs, %d channels, %d samples/channel/epoch (%6.1f seconds/epoch, %6.1f seconds total)\n", DSName, E, M, T, time, ttime);
    }

    // sanity check
    if (M == 0 || E == 0 || TT == 0)
        Cleanup("no valid data");

    NumCovs = MARKER1_;             // 4: global, orient, sum, and noise
    // if req'd, read markers
    if (Params.NumMark > 0) {
        NumCovs += Params.NumMark;
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
        GetFilePath(DSDIR, fpath, sizeof(fpath), &Params, "MarkerFile.mrk", 0);
        GetMarkers(fpath, &Marker, &NumMarkers);
#endif
        if (vflg) {
            printf(" - found %d markers", NumMarkers);
            fflush(stdout);
        }
    }

    // initialize Stats flags
    Stats = new_array(COV_STATS, NumCovs);
    for (n=0; n<NumCovs; n++)
        Stats[n].Used = FALSE;

    // always use CovBand frequencies
    if (vflg) {
        printf(" - done\n");
        printf("selecting filter parameters");
        fflush(stdout);
    }

    if (Params.CovHP == -999. || Params.CovLP == -999.)
        Cleanup("CovBand bandpass parameters required");
    if (Params.CovHP < Channel[f].AnalogHpFreq)
        Cleanup("CovBand highpass cutoff frequency out of recorded range");
    if (Params.CovLP > Channel[f].AnalogLpFreq)
        Cleanup("CovBand lowpass cutoff frequency out of recorded range");
    if (Params.ImageHP < Channel[f].AnalogHpFreq)
        Cleanup("ImageBand highpass cutoff frequency out of recorded range");
    if (Params.ImageLP > Channel[f].AnalogLpFreq)
        Cleanup("ImageBand lowpass cutoff frequency out of recorded range");
    if (Params.FilterType == IIR)
        strcpy(FilterName, "IIR");
    else
        strcpy(FilterName, "FFT");

    // generate covariance subdirectory prefix string by concatenating CovName & bandpass -- this will be the subdirectory name
    sprintf(fpath, "%s,%-d-%-dHz", CovName, (int)Params.CovHP, (int)Params.CovLP);

    // create covariance subdirectory
    sprintf(COVpath, "%s/%s", SAMpath, fpath);
    if (mkdir(COVpath, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1)
        if (errno != EEXIST)
            Cleanup("can't create covariance subdirectory");

    // default covariance matrices: GLOBAL & ORIENT
    sprintf(fpath, "%s/SAMcov-Segments", COVpath);
    fp = fileopen(fpath, "w");

    // initialize GLOBAL
    strcpy(Stats[GLOBAL_].CovFileName, "Global");
    Stats[GLOBAL_].Cov = new_matrixE(M, M, "Stats[GLOBAL_].Cov");
    Stats[GLOBAL_].NumSegments = E;
    Stats[GLOBAL_].Segment = new_arrayE(COV_SEG, E, "GlobalSegment[]");
    fprintf(fp, "'Global' segments:\n");
    for (e=0; e<E; e++) {
        Stats[GLOBAL_].Segment[e].Epoch = e;
        Stats[GLOBAL_].Segment[e].TS = 0;
        Stats[GLOBAL_].Segment[e].TE = T - 1;
        Stats[GLOBAL_].Segment[e].TT = T;
        fprintf(fp, "\t%d: e=%d, start=%d, end=%d\n", E, Stats[GLOBAL_].Segment[e].Epoch, Stats[GLOBAL_].Segment[e].TS, Stats[GLOBAL_].Segment[e].TE);
    }
    Stats[GLOBAL_].Used = TRUE;

    // initialize ORIENT
    strcpy(Stats[ORIENT_].CovFileName, "Orient");
    Stats[ORIENT_].Cov = new_matrixE(M, M, "Stats[ORIENT_].Cov");
    Stats[ORIENT_].NumSegments = E;
    Stats[ORIENT_].Segment = new_arrayE(COV_SEG, E, "OrientSegment[]");
    fprintf(fp, "'Orient' segments:\n");
    for (e=0; e<E; e++) {
        Stats[ORIENT_].Segment[e].Epoch = e;
        Stats[ORIENT_].Segment[e].TS = 0;
        Stats[ORIENT_].Segment[e].TE = T - 1;
        Stats[ORIENT_].Segment[e].TT = T;
        fprintf(fp, "\t%d: e=%d, start=%d, end=%d\n", E, Stats[ORIENT_].Segment[e].Epoch, Stats[ORIENT_].Segment[e].TS, Stats[ORIENT_].Segment[e].TE);
    }
    Stats[ORIENT_].Used = TRUE;

    // check for marker specifications
    j = 0;
    SumSegments = 0;
    if (Params.NumMark > 0) {                   // if there are markers...

        // initialize each marker
        for (i=0; i<Params.NumMark; i++) {      // pay attention! 'i' is the Params.Marker (conditions) index & 'j' is the Stats index
            j = MARKER1_ + i;                   // begin at marker1

            s = Params.Marker[i].MarkName;
            if (Params.Marker[i].MarkName2) {
                s = Params.Marker[i].MarkName2;
            }
            strcpy(Stats[j].CovFileName, s);
            Stats[j].Cov = new_matrixE(M, M, "Stats[Marker].Cov");

            fprintf(fp, "'%s' segments:\n", Params.Marker[i].MarkName);
            if (Params.Marker[i].SegFile) {
                GetFilePath(DSDIR, fpath, sizeof(fpath), &Params, Params.Marker[i].FileName, TRUE);
                Stats[j].Segment = read_segfile(fpath, &n, &Header);
                for (k = 0; k < n; k++) {
                    if (Stats[j].Segment[k].TS >= 0 && Stats[j].Segment[k].TE < T && Stats[j].Segment[k].TT > 0) {
                        fprintf(fp, "\t%d: e=%d, start=%d, end=%d\n", k, Stats[j].Segment[k].Epoch, Stats[j].Segment[k].TS, Stats[j].Segment[k].TE);
                    } else {
                        Cleanup("invalid segment %d in SegFile %s", k, Params.Marker[i].FileName);
                    }
                }
                Stats[j].NumSegments = n;
            } else {
                for (k=TotSegments=0; k<NumMarkers; k++)
                    if (!strcmp(Params.Marker[i].MarkName, Marker[k].Name))
                        TotSegments++;
                if (TotSegments == 0) {
                    fprintf(stderr, "\nmarker name '%s' not found\n", Params.Marker[i].MarkName);
                    Cleanup("there are no segments corresponding to this marker");
                }
                Stats[j].Segment = new_arrayE(COV_SEG, TotSegments, "Segment[]");

                if (Params.DataSegStart != -999. && Params.DataSegEnd != -999.) {   // note that DataSegment only changes covariance segment limits -- not image limits
                    Params.Marker[i].MarkStart = Params.DataSegStart;
                    Params.Marker[i].MarkEnd = Params.DataSegEnd;
                }

                // now search through all the markers for a match to the named marker
                for (n=k=0; n<NumMarkers; n++) {
                    if (!strcmp(Params.Marker[i].MarkName, Marker[n].Name)) {
                        e = Marker[n].Epoch;
                        mark = Marker[n].Latency;
                        if (Epoch[e].Flag) {
                            Stats[j].Segment[k].Epoch = e;
                            Stats[j].Segment[k].TS = (int)rint((mark + Params.Marker[i].MarkStart) * Header.SampleRate);
                            Stats[j].Segment[k].TS += Header.PreTrigSamples;
                            Stats[j].Segment[k].TE = (int)rint((mark + Params.Marker[i].MarkEnd) * Header.SampleRate);
                            Stats[j].Segment[k].TE += Header.PreTrigSamples;
                            Stats[j].Segment[k].TT = Stats[j].Segment[k].TE - Stats[j].Segment[k].TS;
                            if (Stats[j].Segment[k].TS >= 0 && Stats[j].Segment[k].TE < T && Stats[j].Segment[k].TT > 0) {
                                fprintf(fp, "\t%d: e=%d, start=%d, end=%d\n", k, Stats[j].Segment[k].Epoch, Stats[j].Segment[k].TS, Stats[j].Segment[k].TE);
                                k++;
                            }
                        }
                    }
                }
                if (k == 0) {
                    fprintf(stderr, "\nmarker '%s'\n", Params.Marker[i].MarkName);
                    Cleanup("no valid segments found within sample limits");
                }
                Stats[j].NumSegments = k;
            }
            if (Params.Marker[i].Sum) {
                SumSegments += Stats[j].NumSegments;
            }
            Stats[j].Used = TRUE;
        }
    }

    // check for marker specifications for sum covariance
    if (SumSegments > 0) {

        // initialize sum covariance
        strcpy(Stats[SUM_].CovFileName, "Sum");
        Stats[SUM_].Cov = new_matrixE(M, M, "Stats[SUM_].Cov");
        Stats[SUM_].Segment = new_arrayE(COV_SEG, SumSegments, "SumSegment[]");
        Stats[SUM_].NumSegments = SumSegments;

        // get the segments for each marker
        fprintf(fp, "'Sum' segments:\n");
        k = 0;
        for (i = 0; i < Params.NumMark; i++) {
            if (Params.Marker[i].Sum) {
                j = MARKER1_ + i;
                for (n = 0; n < Stats[j].NumSegments; n++) {
                    Stats[SUM_].Segment[k].Epoch = Stats[j].Segment[n].Epoch;
                    Stats[SUM_].Segment[k].TS = Stats[j].Segment[n].TS;
                    Stats[SUM_].Segment[k].TE = Stats[j].Segment[n].TE;
                    Stats[SUM_].Segment[k].TT = Stats[j].Segment[n].TT;
                    fprintf(fp, "\t%d: e=%d, start=%d, end=%d\n", k, Stats[SUM_].Segment[k].Epoch, Stats[SUM_].Segment[k].TS, Stats[SUM_].Segment[k].TE);
                    k++;
                }
            }
        }
        Stats[SUM_].Used = TRUE;
    }

    // if req'd, initialize NOISE
    if (Params.NoiseHP != -999. && Params.NoiseLP != -999.) {
        if (Params.NoiseHP < Channel[f].AnalogHpFreq)
            Cleanup("NoiseBand highpass cutoff frequency outside recorded range");
        if (Params.NoiseLP > Channel[f].AnalogLpFreq)
            Cleanup("NoiseBand lowpass cutoff frequency outside recorded range");
        strcpy(Stats[NOISE_].CovFileName, "Noise");
        Stats[NOISE_].Cov = new_matrixE(M, M, "Stats[NOISE_].Cov");
        Stats[NOISE_].NumSegments = E;
        Stats[NOISE_].Segment = new_arrayE(COV_SEG, E, "NoiseSegment[]");
        fprintf(fp, "'Noise' segments:\n");
        for (e=0; e<E; e++) {
            Stats[NOISE_].Segment[e].Epoch = e;
            Stats[NOISE_].Segment[e].TS = 0;
            Stats[NOISE_].Segment[e].TE = T - 1;
            Stats[NOISE_].Segment[e].TT = T;
            fprintf(fp, "\t%d: e=%d, start=%d, end=%d\n", E, Stats[NOISE_].Segment[e].Epoch, Stats[NOISE_].Segment[e].TS, Stats[NOISE_].Segment[e].TE);
        }
        Stats[NOISE_].Used = TRUE;
    }
    fclose(fp); // close segment list file

    // generate labels, & allocate arrays & matrices for output covariance files
    if (vflg) {
        printf(" - done\n");
        printf("covariance file names:\n");
        for (i=0; i<NumCovs; i++)
            if (Stats[i].Used)
                printf("\t'%s'\n", Stats[i].CovFileName);
    }

    // allocate memory for arrays
    if (vflg) {
        printf("allocating memory");
        fflush(stdout);
    }
    A = new_matrixE(M, M, "A[]");
    if ((x = (double ***)malloc((size_t)E * sizeof(double **))) == NULL)
        allocfailed("x[]");
    for (e=0; e<E; e++) {
        x[e] = new_matrixE(M, T, "x[][]");
    }
    In = new_arrayE(double, T, "In[]");
    Out = new_arrayE(double, T, "Out[]");
    seg = new_matrixE(M, T, "seg[][]");
    U = gsl_matrix_alloc(M, M);
    S = gsl_vector_alloc(M);
    V = gsl_matrix_alloc(M, M);

    // if req'd, determine filter type & generate filter coefficients
    if (Params.CovHP != -999. && Params.CovLP != -999.) {
        if (Params.CovHP < Channel[f].AnalogHpFreq)
            Cleanup("Covariance highpass cutoff out of recorded range");
        if (Params.CovLP > Channel[f].AnalogLpFreq)
            Cleanup("Covariance lowpass cutoff out of recorded range");
        if (vflg) {
            printf(" - done\n");
            printf("computing %6.3f to %6.3f Hz %s filter for covariance", Params.CovHP, Params.CovLP, FilterName);
            fflush(stdout);
        }
        if (Params.FilterType == IIR) {
            mkfilter(&CovIIR, SIG_ORDER, Header.SampleRate, Params.CovHP, Params.CovLP, &CovBW);
        } else {
            if (Params.CovHP > 0.)
                Butterworth(&CovFFT, T, Params.CovHP, Params.CovLP, Header.SampleRate, BANDPASS, MAXORDER, &CovBW);
            else
                Butterworth(&CovFFT, T, 0., Params.CovLP, Header.SampleRate, LOWPASS, MAXORDER, &CovBW);
        }
    }

    // test to see if additional filtering required for orientation
    if (Params.OrientHP != -999. && Params.OrientLP != -999. &&
        Params.OrientHP != Params.CovHP && Params.OrientLP != Params.CovLP) {
        if (Params.OrientHP < Channel[f].AnalogHpFreq)
            Cleanup("Orientation highpass cutoff out of recorded range");
        if (Params.OrientLP > Channel[f].AnalogLpFreq)
            Cleanup("Orientation lowpass cutoff out of recorded range");
        if (vflg) {
            printf(" - done\n");
            printf("computing %6.3f to %6.3f Hz %s filter for orientation covariance", Params.OrientHP, Params.OrientLP, FilterName);
            fflush(stdout);
        }
        if (Params.FilterType == IIR)
            mkfilter(&OrientIIR, SIG_ORDER, Header.SampleRate, Params.OrientHP, Params.OrientLP, &OrientBW);
        else
            Butterworth(&OrientFFT, T, Params.OrientHP, Params.OrientLP, Header.SampleRate, BANDPASS, MAXORDER, &OrientBW);
        oflg = TRUE;        // need to refilter for orientation covariance!
    } else {
        OrientBW = CovBW;
    }

    // if req'd, generate "noise" filter coefficients
    if (Params.NoiseHP > 0. && Params.NoiseLP > 0.) {
        if (vflg) {
            printf(" - done\n");
            printf("computing %8.3f to %8.3f Hz %s filter for measuring noise covariance", Params.NoiseHP, Params.NoiseLP, FilterName);
            fflush(stdout);
        }
        if (Params.FilterType == IIR)
            mkfilter(&NoiseIIR, SIG_ORDER, Header.SampleRate, Params.NoiseHP, Params.NoiseLP, &NoiseBW);
        else
            Butterworth(&NoiseFFT, T, Params.NoiseHP, Params.NoiseLP, Header.SampleRate, BANDPASS, MAXORDER, &NoiseBW);
        nflg = TRUE;
    }

    // compute notch filter
    if (Params.Notch) {
        if (Params.FilterType == FFT) {
            if (vflg) {
                printf(" - done\n");
                printf("making FFT notch filter at %g Hz", Params.Hz);
                fflush(stdout);
            }
            mknotch(&NotchFFT, Params.CovHP, Params.CovLP, Header.SampleRate, T, Params.Hz);
        } else {
            fprintf(stderr, "\nreminder: when using IIR filters, data should be notch filtered for power mains prior to analysis\n");
        }
    }

    // compute covariance for time-segments selected around each marker
    if (vflg) {
        printf(" - done\n");
        printf("reading & filtering data");
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
            d = Out;
            if (Params.FilterType == IIR) {
                bdiir(In, Out, T, &CovIIR);
            } else {
                FFTfilter(In, Out, &CovFFT, T);
                if (Params.Notch) {
                    d = In;
                    FFTfilter(Out, In, &NotchFFT, T);
                }
            }
            for (t=0; t<T; t++)
                x[e][m][t] = d[t];
        }

        // epoch heartbeat
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }

    // for each covariance matrix type (except ORIENT_ & NOISE_, which could use different filters)...
    if (vflg) {
        printf(" - done\n");
        printf("computing all marker covariance matrices:\n");
        fflush(stdout);
    }
    for (n=0; n<NumCovs; n++)
        if (Stats[n].Used && n != ORIENT_ && n != NOISE_) {

            // enumerate
            if (vflg) {
                printf("\t'%s'\n", Stats[n].CovFileName);
            }

            // zero matrices & sample counts
            for (i=0; i<M; i++)
                for (j=0; j<M; j++)
                    Stats[n].Cov[i][j] = 0.;
            Stats[n].NumSamples = 0;

            // for each segment...
            for (k=0; k<Stats[n].NumSegments; k++) {

                // compute covariance for this segment
                e = Stats[n].Segment[k].Epoch;

                // copy segment, compute mean, & subtract from segment
                for (m=0; m<M; m++) {
                    for (t=0, tt=Stats[n].Segment[k].TS, tmp=0.; tt<=Stats[n].Segment[k].TE; t++, tt++) {
                        seg[m][t] = x[e][m][tt];
                        tmp += seg[m][t];
                    }
                    mean = tmp / (double)Stats[n].Segment[k].TT;
                    for (t=0; t<(double)Stats[n].Segment[k].TT; t++)
                        seg[m][t] -= mean;
                }

                // integrate upper triangle
                for (i=0; i<M; i++)
                    for (j=i; j<M; j++) {
                        for (t=0, tmp=0.; t<Stats[n].Segment[k].TT; t++)
                            tmp += seg[i][t] * seg[j][t];
                        A[i][j] = tmp;      // unnormalized at this point
                    }

                // copy to lower triangle
                for (i=1; i<M; i++)
                    for (j=0; j<i; j++)
                        A[i][j] = A[j][i];

                // integrate segment into appropriate covariance matrix
                for (i=0; i<M; i++)
                    for (j=0; j<M; j++)
                        Stats[n].Cov[i][j] += A[i][j];

                // update sample count
                Stats[n].NumSamples += Stats[n].Segment[k].TT;
            }

            // normalize covariance matrices
            if (Stats[n].NumSamples > 0)
                for (i=0; i<M; i++)
                    for (j=0; j<M; j++)
                        Stats[n].Cov[i][j] /= (double)Stats[n].NumSamples;
            if (Stats[n].NumSamples < M)
                Cleanup("too few samples -- covariance is singular");
        }

    // compute orientation (Global) covariance matrix
    if (oflg) {             // need to refilter the data for the orientation covariance
        if (vflg) {
            printf("reading & filtering data");
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
                d = Out;
                if (Params.FilterType == IIR) {
                    bdiir(In, Out, T, &OrientIIR);
                } else {
                    FFTfilter(In, Out, &OrientFFT, T);
                    if (Params.Notch) {
                        d = In;
                        FFTfilter(Out, In, &NotchFFT, T);
                    }
                }
                for (t=0; t<T; t++)
                    x[e][m][t] = d[t];
            }

            // epoch heartbeat
            if (vflg) {
                printf(".");
                fflush(stdout);
            }
        }

        // zero orientation covariance matrix & sample count
        for (i=0; i<M; i++)
            for (j=0; j<M; j++)
                Stats[ORIENT_].Cov[i][j] = 0.;
        Stats[ORIENT_].NumSamples = 0;

        // for each segment..
        if (vflg) {
            printf(" - done\n");
            printf("computing orientation covariance matrix:\n");
            printf("\t'%s'\n", Stats[ORIENT_].CovFileName);
            fflush(stdout);
        }
        for (k=0; k<Stats[ORIENT_].NumSegments; k++) {

            // compute covariance for this segment
            e = Stats[ORIENT_].Segment[k].Epoch;

            // copy segment, compute mean, & subtract from segment
            for (m=0; m<M; m++) {
                for (t=0, tt=Stats[ORIENT_].Segment[k].TS, tmp=0.; tt<=Stats[ORIENT_].Segment[k].TE; t++, tt++) {
                    seg[m][t] = x[e][m][tt];
                    tmp += seg[m][t];
                }
                mean = tmp / (double)Stats[ORIENT_].Segment[k].TT;
                for (t=0; t<(double)Stats[ORIENT_].Segment[k].TT; t++)
                    seg[m][t] -= mean;
            }

            // integrate upper triangle
            for (i=0; i<M; i++)
                for (j=i; j<M; j++) {
                    for (t=0, tmp=0.; t<Stats[ORIENT_].Segment[k].TT; t++)
                        tmp += seg[i][t] * seg[j][t];
                    A[i][j] = tmp;      // unnormalized at this point
                }

            // copy to lower triangle
            for (i=1; i<M; i++)
                for (j=0; j<i; j++)
                    A[i][j] = A[j][i];

            // integrate segment into appropriate covariance matrix
            for (i=0; i<M; i++)
                for (j=0; j<M; j++)
                    Stats[ORIENT_].Cov[i][j] += A[i][j];

            // update sample count
            Stats[ORIENT_].NumSamples += Stats[ORIENT_].Segment[k].TT;
        }

        // normalize orientation covariance matrix
        if (Stats[ORIENT_].NumSamples > 0)
            for (i=0; i<M; i++)
                for (j=0; j<M; j++)
                    Stats[ORIENT_].Cov[i][j] /= (double)Stats[ORIENT_].NumSamples;
        if (Stats[ORIENT_].NumSamples < M)
            Cleanup("too few samples -- covariance is singular");

    } else {        // oflg == FALSE -- just copy GLOBAL covariance
        for (i=0; i<M; i++)
            for (j=0; j<M; j++)
                Stats[ORIENT_].Cov[i][j] = Stats[GLOBAL_].Cov[i][j];
        Stats[ORIENT_].NumSamples = Stats[GLOBAL_].NumSamples;
    }

    // if req,d, compute noise covariance matrix
    if (nflg) {
        if (vflg) {
            printf("reading & filtering data");
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
                d = Out;
                if (Params.FilterType == IIR) {
                    bdiir(In, Out, T, &NoiseIIR);
                } else {
                    FFTfilter(In, Out, &NoiseFFT, T);
                    if (Params.Notch) {
                        d = In;
                        FFTfilter(Out, In, &NotchFFT, T);
                    }
                }
                for (t=0; t<T; t++)
                    x[e][m][t] = d[t];
            }

            // epoch heartbeat
            if (vflg) {
                printf(".");
                fflush(stdout);
            }
        }

        // zero matrices & sample counts
        for (i=0; i<M; i++)
            for (j=0; j<M; j++)
                Stats[NOISE_].Cov[i][j] = 0.;
        Stats[NOISE_].NumSamples = 0;

        // for each segment...
        if (vflg) {
            printf(" - done\n");
            printf("computing noise covariance matrix:\n");
            printf("\t'%s'\n", Stats[NOISE_].CovFileName);
            fflush(stdout);
        }
        for (k=0; k<Stats[NOISE_].NumSegments; k++) {

            // compute covariance for this segment
            e = Stats[NOISE_].Segment[k].Epoch;

            // copy segment, compute mean, & subtract from segment
            for (m=0; m<M; m++) {
                for (t=0, tt=Stats[NOISE_].Segment[k].TS, tmp=0.; tt<=Stats[NOISE_].Segment[k].TE; t++, tt++) {
                    seg[m][t] = x[e][m][tt];
                    tmp += seg[m][t];
                }
                mean = tmp / (double)Stats[NOISE_].Segment[k].TT;
                for (t=0; t<(double)Stats[NOISE_].Segment[k].TT; t++)
                    seg[m][t] -= mean;
            }

            // integrate upper triangle
            for (i=0; i<M; i++)
                for (j=i; j<M; j++) {
                    for (t=0, tmp=0.; t<Stats[NOISE_].Segment[k].TT; t++)
                        tmp += seg[i][t] * seg[j][t];
                    A[i][j] = tmp;      // unnormalized at this point
                }

            // copy to lower triangle
            for (i=1; i<M; i++)
                for (j=0; j<i; j++)
                    A[i][j] = A[j][i];

            // integrate segment into appropriate covariance matrix
            for (i=0; i<M; i++)
                for (j=0; j<M; j++)
                    Stats[NOISE_].Cov[i][j] += A[i][j];

            // update sample count
            Stats[NOISE_].NumSamples += Stats[NOISE_].Segment[k].TT;
        }

        // normalize covariance matrix
        if (Stats[NOISE_].NumSamples > 0)
            for (i=0; i<M; i++)
                for (j=0; j<M; j++)
                    Stats[NOISE_].Cov[i][j] /= (double)Stats[NOISE_].NumSamples;
        if (Stats[NOISE_].NumSamples < M)
            Cleanup("too few samples -- covariance is singular");
    }

    // check for dead channels in global covariance matrix
    for (i=0; i<M; i++) {
        j = Header.PsIndex[i];
        tmp = 1.0e+15 * sqrt(Stats[GLOBAL_].Cov[i][i] / CovBW); // fT/rt-Hz
        if (tmp < 1.)    // less than 1 fT/rt-Hz based on variance indicates a dead channel
            fprintf(stderr, "\n***** Bad Channel '%s': total signal is %f fT/rt-Hz\n", Channel[j].ChannelName, tmp);
    }

    // covariance statistics
    if (vflg) {
        printf("covariance statistics:\n");
        fflush(stdout);
    }
    for (n=0; n<NumCovs; n++)
        if (Stats[n].Used) {

            // enumerate
            if (vflg) {
                printf("\t'%s'\n", Stats[n].CovFileName);
            }

            // report effective number of samples
            switch (n) {
                case ORIENT_: BW = OrientBW; break;
                case NOISE_: BW = NoiseBW; break;
                default: BW = CovBW; break;
            }
            EffSamples = (2. * BW * ((double)Stats[n].NumSamples) / Header.SampleRate) / (double)M;
            eflg = (EffSamples < 3.) ? TRUE : FALSE;
            if (vflg) {
                printf("\t\teffective samples is %8.3f x number of channels\n", EffSamples);
            }
            if (eflg) {
                fflush(stdout);
                msg("\07Warning! The effective number of samples is too low for reliable SAM analysis\n");
            }

            // perform an SVD to check for degeneracy of each covariance matrix
            for (i=0; i<M; i++)
                for (j=0; j<M; j++)
                    gsl_matrix_set(U, i, j, Stats[n].Cov[i][j]);
            if (gsl_linalg_SV_decomp_jacobi(U, V, S) == -1)
                cleanup("gsl_linalg_SV_decomp_jacobi() failed");
            Tolerance = gsl_vector_get(S, M-1) / gsl_vector_get(S, 0);

            // save singular values file
            sprintf(fpath, "%s/%s_Eigenvalues", COVpath, Stats[n].CovFileName);
            fp = fileopen(fpath, "w");
            fprintf(fp, "Covariance File: '%s'\n", Stats[n].CovFileName);
            fprintf(fp, "Number of Segments: %d\n", Stats[n].NumSegments);
            fprintf(fp, "Number of Samples: %d\n", Stats[n].NumSamples);
            fprintf(fp, "Singular Values:\n");
            for (i=0; i<M; i++)
                fprintf(fp, "%e\n", gsl_vector_get(S, i));
            fprintf(fp, "Tolerance=%e\n", Tolerance);
            fprintf(fp, "Condition=%e\n", 1. / Tolerance);
            fclose(fp);

            // look at derivatives - find k for minimum derivative
            for (k=M-ndim, min=HUGE; k>=2; k--) {
                delta = gsl_vector_get(S, k-2) - gsl_vector_get(S, k+2);
                if (delta < min) {
                    min = delta;
                    j = k;
                }
            }
            Stats[n].Noise = gsl_vector_get(S, j);
            LSNoise = gsl_vector_get(S, M-1);

            // warnings!
            if (Tolerance < 0. || fabs(1. / Tolerance) > 1.0e+15 || sqrt(LSNoise / BW) < 5.0e-16) {
                printf("\07Warning! '%s' appears to be degenerate\n", Stats[n].CovFileName);
                printf("\tcovariance matrix condition is %e, with least-significant singular value = %e\n", 1. / Tolerance, LSNoise);
                printf("\testimated measurement noise is %10.6f fT/rt-Hz\n", 1.0e+15 * sqrt(LSNoise / BW));
            } else {
                if (vflg) {
                    printf("\t\teffective noise floor at rank %d: %10.6f fT/rt-Hz\n", j, 1.0e+15 * sqrt(Stats[n].Noise / BW));
                    printf("\t\teffective noise floor at least-significant singular value: %10.6f fT/rt-Hz\n", 1.0e+15 * sqrt(LSNoise / BW));
                }
            }

            // save mean noise values file
            sprintf(fpath, "%s/%s_Noise", COVpath, Stats[n].CovFileName);
            fp = fileopen(fpath, "w");
            fprintf(fp, "%e\n", Stats[n].Noise);
            fclose(fp);
        }

    // free gsl matrices & vectors
    gsl_matrix_free(U);
    gsl_vector_free(S);
    gsl_matrix_free(V);

    // write req'd number of covariance files
    if (vflg) {
        printf("writing covariance files to disk");
        fflush(stdout);
    }
    for (n=0; n<NumCovs; n++) {
        memset(&CovHeader, 0, sizeof(CovHeader));
        if (Stats[n].Used && Stats[n].NumSegments > 0) {

            // set up generic covariance header & write the covariance matrix
            CovHeader.Version = SAM_REV;
            sprintf(CovHeader.SetName, "%s", Header.SetName);
            if (mflg)
                sprintf(CovHeader.SpecName, "%s", CovName);
            else
                sprintf(CovHeader.SpecName, "NONE");
            CovHeader.NumChans = M;
            switch (n) {
                case ORIENT_:
                    CovHeader.HPFreq = Params.OrientHP;
                    CovHeader.LPFreq = Params.OrientLP;
                    CovHeader.BWFreq = OrientBW;
                    break;
                case NOISE_:
                    CovHeader.HPFreq = Params.NoiseHP;
                    CovHeader.LPFreq = Params.NoiseLP;
                    CovHeader.BWFreq = NoiseBW;
                    break;
                default:
                    CovHeader.HPFreq = Params.CovHP;
                    CovHeader.LPFreq = Params.CovLP;
                    CovHeader.BWFreq = CovBW;
                    break;
            }
            CovHeader.CovType = n;
            CovHeader.NumSegments = Stats[n].NumSegments;
            CovHeader.NumSamples = Stats[n].NumSamples;
            CovHeader.Noise = Stats[n].Noise;
            sprintf(fpath, "%s/%s.cov", COVpath, Stats[n].CovFileName);
            PutCov(fpath, &CovHeader, Header.PsIndex, Stats[n].Cov);
        }
    }

    log_params(SAMpath);

    // say goodbye...
    if (vflg) {
        printf(" - done\n'%s' done\n", Progname);
    }

    exit(0);
}
