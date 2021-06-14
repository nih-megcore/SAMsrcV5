// SAM3d -- synthetic aperture magnetometry imaging for source event related
//  time-locked signals. SAM3d reads all SAM weights for volume image, & computes
//  the source waveform for each voxel. Given specifications for active & control
//  time segments, SAM3d computes the mean power & variance relative for the
//  specified event windows. The results are stored in NIFTI image ('.nii') files.
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
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
#define MINOR_REV   1               // new command line parser
#define SIG_ORDER   4

typedef struct {
    char        Name[64];           // marker name
    gsl_matrix  *Wgt;               // Wgt - VxM SAM weights
    double      Noise;              // noise power
    SAM_HDR     WgtHdr;             // SAM weight file header
    COV_SEG     *Segment;           // Segment[TotSegments] -- covariance integration time-segment specs
    int         NumSegments;        // number of segments for this condition
} IMAGE_STAT;

int main(
    int             argc,
    char            **argv
) {
    IIRSPEC         DataIIR;        // IIR filter structure
    FFTSPEC         DataFFT;        // FFT filter structure
    FFTSPEC         NotchFFT;       // notch FFT filter structure
    IMAGE_STAT      *Stats;         // Stats[N] -- structure for conditions (one for each marker or one if no markers)
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
    gsl_vector      *W;             // W - Mx1 beamformer coefficients
    GIIATLAS        *atlas;         // gifti atlas info
    double          ***x;           // x[E][T][M] -- filtered multisensor data
    double          *In;            // In[T] -- input time-series
    double          *Out;           // Out[T] -- output time-series
    double          *da;            // double array
    double          mark;           // time mark
    double          BWFreq;         // effective bandwidth
    double          BWratio;        // ratio of covariance bandwidth to image bandwidth (nominal, only)
    double          Src;            // source metric
    double          Noise;          // effective noise
    double          percent;        // percent done
    double          a1;             // sum of metric
    double          a2;             // sum of metric squared
    float           *ImgM;          // ImgM[V] -- image mean
    float           *ImgV;          // ImgV[V] -- image variance
    register double tmp;            // what else?
    int             *ChanIndex;     // ChanIndex[M] -- index of channels used
    int             E;              // number of good epochs
    int             M;              // number of MEG or reference channels used
    int             N;              // just a number
    int             T;              // number of samples
    int             TS;             // segment start
    int             TE;             // segment end
    int             TT;             // segment length
    int             V;              // number of voxels
    int             v;              // voxel index
    int             c;              // virtual channel
    int             pct;            // integer percentage
    int             m;              // channel index
    int             e;              // epoch index
    int             f;              // 1st primary sensor channel index
    int             t;              // sample index
    int             tt;             // alternate sample index
    int             i;              // matrix index
    int             j;              // ...and another
    int             k;              // this is the last one -- for sure...
    int             n;              // condition index
    int             NumTot;         // number of total segments
    int             NumMarkers;     // number of markers in 'MarkerFile.mrk'
    int             NumConditions;  // number of conditions in output image
    int             *vidx;          // array for shuffling voxel indices
    int             eflg = FALSE;   // command-line error flag
    int             mflg = FALSE;   // parameter file name flag
    int             nflg = FALSE;   // use noise covariance file
    int             rflg = FALSE;   // run name flag
    int             gflg = FALSE;   // gifti format Atlas
    int             vflg = FALSE;   // verbose mode flag
    time_t          secs;           // random number seed
    char            *s;             // general string
    char            fpath[256];     // general path name
    char            *DSName = NULL; // MEG dataset name
    char            DSpath[256];    // MEG dataset path
    char            SAMpath[256];   // SAM subdirectory path
    char            AtlasPath[256]; // atlas path
    char            WgtDir[256];    // weight directory path name
    char            *CovName = NULL; // input covariance file name
    char            *WtsName = NULL; // input weights file name
    char            *OutName = NULL; // output file name
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
    reg_parm("SegFile");
    reg_parm("CovType");
    reg_parm("CovBand");
    reg_parm("ImageBand");
    reg_parm("NoiseBand");
    reg_parm("FilterType");
    reg_parm("Notch");
    reg_parm("Hz");
    reg_parm("AtlasName");
    reg_parm("ImageFormat");
    reg_parm("ImageMetric");
    set_parm_arg_help("ImageMetric", "Power");
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
        printf("3 D   S A M   I M A G I N G\t%s\n", PRG_REV);
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
    f = Header.PsIndex[0];

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

    // announce metric type
    if (Params.ImageMetric == -1)
        Cleanup("parameters must specify 'ImageMetric'");
    if (vflg) {
        printf(" - done\n");
        printf("Imaging functions of ");
        switch (Params.ImageMetric) {
            case POWER: printf("power\n"); break;
            default: Cleanup("selected ImageMetric not supported"); break;
        }
    }

    // check imaging parameters
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
    BWratio = (Params.ImageLP - Params.ImageHP) / (Params.CovLP - Params.CovHP);
    if (BWratio > 1.)
        Cleanup("covariance bandwidth must be >= image bandwidth");
    if (Params.FilterType == IIR)
        strcpy(FilterName, "IIR");
    else
        strcpy(FilterName, "FFT");

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
    } else {        // if no markers, then allocate one Stats structure for Global measures
        if (vflg) {
            printf("no trigger markers");
            fflush(stdout);
        }
        NumConditions = 1;
        Stats = new_arrayE(IMAGE_STAT, NumConditions, "Stats[]");
        strcpy(Stats[0].Name, "Global");      // name it!
        Params.CovType = GLOBAL_;
    }

    // populate segment array
    if (vflg) {
        printf(" - done\n");
        printf("populating segments:\n");
    }

    if (Params.NumMark > 0) {       // if there are markers...

        // for each marker...
        for (n=0; n<NumConditions; n++) {
            if (vflg) {
                printf("\t'%s'\n", Stats[n].Name);
            }

            if (Params.Marker[n].SegFile) {
                GetFilePath(DSDIR, fpath, sizeof(fpath), &Params, Params.Marker[n].FileName, TRUE);
                Stats[n].Segment = read_segfile(fpath, &i, &Header);
                for (k = 0; k < i; k++) {
                    if (Stats[n].Segment[k].TS < 0 || Stats[n].Segment[k].TE >= T) {
                        Cleanup("invalid segment %d in SegFile %s", k, Params.Marker[n].FileName);
                    }
                }
                Stats[n].NumSegments = i;
            } else {
                // count marker matches
                s = Params.Marker[n].MarkName;
                for (i=NumTot=0; i<NumMarkers; i++)
                    if (!strcmp(s, Marker[i].Name))
                        NumTot++;
                if (NumTot == 0)
                    Cleanup("no markers found with specified names");

                // allocate segment array
                Stats[n].Segment = new_arrayE(COV_SEG, NumTot, "Stats[].Segment[]");

                // populate Stats for each condition
                for (i=k=0, e=(-1); i<NumMarkers; i++) {
                    if (Marker[i].Epoch != e)
                        e = Marker[i].Epoch;
                    if (Epoch[e].Flag) {

                        // populate segment windows
                        if (!strcmp(s, Marker[i].Name)) {
                            mark = Marker[i].Latency;
                            TS = (int)floor((mark + Params.Marker[n].MarkStart) * Header.SampleRate) + Header.PreTrigSamples;
                            TE = (int)floor((mark + Params.Marker[n].MarkEnd) * Header.SampleRate) + Header.PreTrigSamples;
                            TT = TE - TS + 1;
                            if (TS >= 0 && TE < T) {
                                Stats[n].Segment[k].Epoch = e;
                                Stats[n].Segment[k].TS = TS;
                                Stats[n].Segment[k].TE = TE;
                                Stats[n].Segment[k].TT = TT;
                                k++;
                            }
                        }
                    }
                }
                if (k == 0)
                    Cleanup("there were no valid data segments for specified markers & time windows");
                Stats[n].NumSegments = k;
            }
        }

    } else {

        // allocate segment array
        Stats[0].Segment = new_arrayE(COV_SEG, E, "Stats[].Segment[]");

        // for each epoch, populate segment windows...
        for (e=0; e<E; e++) {
            Stats[0].Segment[e].Epoch = e;
            Stats[0].Segment[e].TS = 0;
            Stats[0].Segment[e].TE = T - 1;
            Stats[0].Segment[e].TT = T;
        }
        Stats[0].NumSegments = E;
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

            if (Params.NumMark > 0) {
                for (n=0; n<NumConditions; n++) {   // all conditions use Global weights & noise
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
            } else {
                if (gflg) {
                    Stats[0].Wgt = GetGIFTIWts(fpath);
                } else {
                    Stats[0].Wgt = GetNIFTIWts(fpath, &NiiHdr, extension, &ExtHdr);
                }
                V = Stats[0].Wgt->size1;
                if (Stats[0].Wgt->size2 != M)
                    Cleanup("number of channels in dataset does not match number in weight file");
                Stats[0].Noise = Noise;
            }
            break;

        case SUM_:          // read Sum weights
            GetNoise(WgtDir, "Sum", &Noise);
            if (Params.ImageFormat == TLRC)
                sprintf(fpath, "%s/Sum_at.nii", WgtDir);
            else
                sprintf(fpath, "%s/Sum.nii", WgtDir);

            for (n=0; n<NumConditions; n++) {   // all conditions use sum weights
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

        case ALL_:          // read weights & noise for each individual marker
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
    ImgM = new_arrayE(float, V, "ImgM[]");
    ImgV = new_arrayE(float, V, "ImgV[]");
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

    // for each condition -- only 1 for GLOBAL & SUM, but multiple for ALL
    for (n=0; n<NumConditions; n++) {

        // announce condition
        if (vflg) {
            printf(" - done\n");
            printf("computing image voxels for '%s':\n", Stats[n].Name);
        }

        // loop on all coordinates in the ROI
        pct = -1;
        N = 0;
#pragma omp parallel private(W, da)
{
        W = gsl_vector_alloc(M);
        da = new_arrayE(double, T, "da[]");

#pragma omp for private(c, m, tmp, i, j, k, a1, a2, e, t, tt, Src, Noise)
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
                for (k=0, a1=a2=0.; k<Stats[n].NumSegments; k++) {

                    // set up segment initialization
                    e = Stats[n].Segment[k].Epoch;
                    tt = Stats[n].Segment[k].TS;

                    // compute SAM time series for segment
                    for (t=0; t<Stats[n].Segment[k].TT; t++, tt++) {
                        for (m=0, tmp=0.; m<M; m++)
                            tmp += x[e][tt][m] * gsl_vector_get(W, m);
                        da[t] = tmp / Noise;
                    }

                    // compute metric
                    switch (Params.ImageMetric) {
                        case POWER:
                            Src = power(da, 0, Stats[n].Segment[k].TT-1);
                            break;
                        default:
                            Cleanup("ImageMetric not implemented");
                            break;
                    }

                    // accumulate segment statistics
                    a1 += Src;
                    a2 += Src * Src;

                }   // for(k...

                // compute voxel values
                ImgM[c] = (Stats[n].NumSegments > 0)? a1 / (double)Stats[n].NumSegments: 0.;
                ImgV[c] = (Stats[n].NumSegments > 0)? (a2 - (a1 * a1 / (double)Stats[n].NumSegments)) / ((double)Stats[n].NumSegments - 1.): 0.;

            } else {        // if(tmp...
                ImgM[c] = 0.;
                ImgV[c] = 0.;
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
} // end parallel

        // write mean & variance images
        if (vflg) {
            printf("\nwriting 'Mean'");
            fflush(stdout);
        }
        if (strcmp(Params.DirName, "NULL") == 0) {
            sprintf(Name1, "%s/%s,%s,%s,", SAMpath, Prefix, OutName, Stats[n].Name);
        } else {
            // make sure the output directory exists
            if (mkdir(Params.DirName, 0755) == -1)
                if (errno != EEXIST)
                    Cleanup("can't create ImageDirectory");
            sprintf(Name1, "%s/%s,%s,%s,", Params.DirName, Prefix, OutName, Stats[n].Name);
        }
        switch (Params.ImageMetric) {
            case POWER:         sprintf(Name2, "%s3D_PWR,", Name1); break;
            default:            Cleanup("ImageMetric not implemented"); break;
        }
        if (gflg) {
            /* Right hemisphere is stored first in ImgM! */

            sprintf(fpath, "%sMean.rh.gii", Name2);
            PutGIFTISurf(fpath, &ImgM, 1, 0, atlas->ridx, atlas->nrh, atlas->meta, 0.);
            sprintf(fpath, "%sMean.lh.gii", Name2);
            PutGIFTISurf(fpath, &ImgM, 1, atlas->nrh, atlas->lidx, atlas->nlh, atlas->meta, 0.);
        } else {
            sprintf(fpath, "%sMean.nii", Name2);
            fp = fileopen(fpath, "w");

            // modify NIFTI header
            memset(extension, 0, 4);                    // for now, ignore the fact that we may have an extended header
            NiiHdr.dim[0] = 3;                          // 3d volume
            NiiHdr.dim[4] = 0;
            NiiHdr.toffset = 0.;                        // time offset shift
            NiiHdr.slice_duration = 1.;
            NiiHdr.slice_code = NIFTI_SLICE_UNKNOWN;
            NiiHdr.xyzt_units = NIFTI_UNITS_MM;         // units
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
            if (fwrite((void *)ImgM, sizeof(float), V, fp) != V)
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
            PutGIFTISurf(fpath, &ImgV, 1, 0, atlas->ridx, atlas->nrh, atlas->meta, 0.);
            sprintf(fpath, "%sVariance.lh.gii", Name2);
            PutGIFTISurf(fpath, &ImgV, 1, atlas->nrh, atlas->lidx, atlas->nlh, atlas->meta, 0.);
        } else {
            sprintf(fpath, "%sVariance.nii", Name2);
            fp = fileopen(fpath, "w");

            // modify NIFTI header
            memset(extension, 0, 4);                    // for now, ignore the fact that we may have an extended header
            NiiHdr.dim[0] = 3;                          // 3d volume
            NiiHdr.dim[4] = 0;
            NiiHdr.toffset = 0.;                        // time offset shift
            NiiHdr.slice_duration = 1.;
            NiiHdr.slice_code = NIFTI_SLICE_UNKNOWN;
            NiiHdr.xyzt_units = NIFTI_UNITS_MM;         // units
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
            if (fwrite((void *)ImgV, sizeof(float), V, fp) != V)
                Cleanup("can't write data image file");
            fclose(fp);
        }
    }           // for(n...

    log_params(SAMpath);

    if (vflg) {
        msg(" - done\n'%s' done\n", Progname);
    }
    if (nflg) {
        gsl_matrix_free(Cn);
    }
    exit(0);
}
