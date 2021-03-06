// SAM3dc -- synthetic aperture magnetometry imaging for source event related
//  time-locked signals. SAM3dc reads all SAM weights for volume image, & computes
//  the source waveform for each voxel. SAM3dc uses continuous functions for power
//  & RVE instead of segmented functions. Given specifications for active & control
//  time segments, SAM3dc computes the mean power & variance relative for the
//  specified event windows. The result is represented stored in a NIFTI image file
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
#include "sam_param.h"
#include "sam_parse.h"

// internal constants
#define MINOR_REV   4               // new command line parser
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
    int     argc,
    char    **argv)
{
    IIRSPEC         DataIIR;        // IIR filter structure
    FFTSPEC         DataFFT;        // FFT filter structure
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
    double          ***x;           // x[E][T][M] -- filtered multisensor data
    double          *In;            // In[TT] -- input time-series
    double          *Out;           // Out[TT] -- output time-series
    double          mark;           // time mark
    double          BWFreq;         // effective bandwidth
    double          BWratio;        // ratio of image bandwidth to covariance bandwidth
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
    int             N;              // just a number to take care of extended header read
    int             T;              // number of samples
    int             TT;             // E*T
    int             V;              // number of voxels weight file
    int             c;              // input character & virtual channel
    int             pct;            // integer percentage
    int             m;              // channel index
    int             e;              // epoch index
    int             f;              // 1st primary sensor channel index
    int             t;              // sample index
    int             ts;             // starting sample index
    int             te;             // ending sample index
    int             tt;             // alternate sample index
    int             i;              // matrix index
    int             j;              // ...and another
    int             k;              // this is the last one -- for sure...
    int             n;              // condition index
    int             NumTot;         // number of total segments
    int             NumMarkers;     // number of markers in 'MarkerFile.mrk'
    static int      eflg = FALSE;   // command-line error flag
    static int      mflg = FALSE;   // parameter file name flag
    static int      nflg = FALSE;   // use noise covariance file
    static int      rflg = FALSE;   // run name flag
    static int      vflg = FALSE;   // verbose mode flag
    char            *s;             // general string
    char            fpath[256];     // general path name
    char            *DSName = NULL; // MEG dataset name
    char            DSpath[256];    // MEG dataset path
    char            SAMpath[256];   // SAM subdirectory path
    char            WgtDir[256];    // weight directory path name
    char            *ParmName = NULL; // parameter file name
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

    reg_parm("ImageFormat");
    reg_parm("ImageMetric");
    set_parm_arg_help("ImageMetric", "Power|RankVectorEntropy");

    reg_parm("ImageDirectory");
    reg_parm("PrefixLength");
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

    ParmName = Params.ParmName;
    p = get_parm("OutName");
    if (p->set) {
        ParmName = (char *)p->ptr;
    }

    // get dataset information structures
    if (vflg) {
        printf("3 D C  S A M   I M A G I N G\t%s\n", PRG_REV);
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
    TT = E * T;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index
    if (vflg) {
        printf(" - done\n");
        printf("'%s' uses %d epochs, %d channels, %d samples/channel/epoch\n", DSName, E, M, T);
        fflush(stdout);
    }

    // sanity check
    if (M == 0 || E == 0 || T == 0)
        Cleanup("no valid data");

    // announce metric type
    if (Params.ImageMetric == -1)
        Cleanup("parameters must specify 'ImageMetric'");
    if (vflg) {
        printf(" - done\n");
        printf("Imaging functions of ");
        switch (Params.ImageMetric) {
            case POWER: printf("power\n"); break;
            case RV_ENTROPY: printf("rank vector entropy for %d dimensions\n", Params.Dims); break;
            default: Cleanup("selected ImageMetric not supported"); break;
        }
    }

    // check imaging parameters
    if (Params.NumMark < 1)
        Cleanup("There must be at least one marker.");
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
    if (Params.FilterType == IIR)
        strcpy(FilterName, "IIR");
    else
        strcpy(FilterName, "FFT");
#if 0
    if (Params.ImageMetric == RV_ENTROPY)
        if (Params.Dims < 4 || Params.Dims > 5)
            Cleanup("RVE dimensions must be 4 or 5");   // @@@
#endif

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
        sprintf(fpath, "%s/%s,%-d-%-dHz/Noise.cov", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
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
        fflush(stdout);
    }

    // populate segment array
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("populating segments & reading weights:\n");
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
                    Stats[n].Segment[k].TS = (int)floor((mark + Params.Marker[n].MarkStart) * Header.SampleRate) + Header.PreTrigSamples;
                    Stats[n].Segment[k].TE = (int)floor((mark + Params.Marker[n].MarkEnd) * Header.SampleRate) + Header.PreTrigSamples;
                    Stats[n].Segment[k].TT = Stats[0].Segment[k].TE - Stats[0].Segment[k].TS; // + 1;
                    if(Stats[n].Segment[k].TS >= 0 && Stats[0].Segment[k].TE < T) {
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
    }
    W = gsl_vector_alloc(M);
    if((ImgM = (float *)malloc((size_t)V * sizeof(float))) == NULL)
        allocfailed("ImgM[]");
    if((ImgV = (float *)malloc((size_t)V * sizeof(float))) == NULL)
        allocfailed("ImgV[]");
    if((x = (double ***)malloc((size_t)E * sizeof(double **))) == NULL)
        allocfailed("x[]");
    for(e=0; e<E; e++) {
        if((x[e] = (double **)malloc((size_t)T * sizeof(double *))) == NULL)
            allocfailed("x[][]");
        for(t=0; t<T; t++)
            if((x[e][t] = (double *)malloc((size_t)M * sizeof(double))) == NULL)
                allocfailed("x[][][]");
    }
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
            if(Params.FilterType == IIR)
                bdiir(In, Out, T, &DataIIR);
            else
                FFTfilter(In, Out, &DataFFT, T);
            for(t=0; t<T; t++)
                x[e][t][m] = Out[t];
        }
        if(vflg == TRUE) {          // if req'd - heartbeat
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
                for(e=0; e<E; e++)
                    for(t=0, tt=e*T; t<T; t++, tt++) {
                        for(m=0, tmp=0.; m<M; m++)
                            tmp += x[e][t][m] * gsl_vector_get(W, m);
                        In[tt] = tmp / Noise;
                    }

                // compute metric
                switch(Params.ImageMetric) {
                    case POWER:
                        hilbert(TT, In, Out);
                        for(t=0; t<TT; t++)
                            Out[t] = 0.5 * (In[t] * In[t] + Out[t] * Out[t]);   // power S/N
                        break;
                    case RV_ENTROPY:
                        Data2RVE(In, Out, TT, Params.Dims, Params.ImageLP, Params.Tau, Header.SampleRate); break;
                        break;
                    default:
                        Cleanup("ImageMetric not implemented");
                        break;
                }

                // integrate source over segment
                for(k=0, a1=a2=0.; k<Stats[n].NumSegments; k++) {

                    // set up segment initialization
                    e = Stats[n].Segment[k].Epoch;
                    ts = e * T + Stats[n].Segment[k].TS;
                    te = e * T + Stats[n].Segment[k].TE;
                    for(t=ts, Src=0.; t<=te; t++)
                        Src += Out[t];
                    Src /= (double)Stats[n].Segment[k].TT;

                    // accumulate segment statistics
                    a1 += Src;
                    a2 += Src * Src;

                }   // for(k...

                // compute voxel values
                ImgM[c] = (Stats[n].NumSegments > 2)? a1 / (double)Stats[n].NumSegments: 0.;
                ImgV[c] = (Stats[n].NumSegments > 2)? (a2 - (a1 * a1 / (double)Stats[n].NumSegments)) / ((double)Stats[n].NumSegments - 1.): 0.;

            } else {        // if(tmp...
                ImgM[c] = 0.;
                ImgV[c] = 0.;
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
        }       // for(c...

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
            case POWER:         sprintf(Name2, "%s3DC_PWR,", Name1); break;
            case RV_ENTROPY:    sprintf(Name2, "%s3DC_RVE,", Name1); break;
            default:            Cleanup("ImageMetric not implemented"); break;
        }
        sprintf(fpath, "%sMean.nii", Name2);
        if((fp = fopen(fpath, "w")) == NULL)
            Cleanup("can't open .nii file for write");

        // modify NIFTI header
        memset(extension, 0, 4);                    // for now, ignore the fact that we may have an extended header
        NiiHdr.dim[0] = 3;                          // 3d volume
        NiiHdr.dim[4] = 0;
        NiiHdr.toffset = 0.;                        // time offset shift
        NiiHdr.slice_duration = 1.;
        NiiHdr.slice_code = NIFTI_SLICE_UNKNOWN;
        NiiHdr.xyzt_units = NIFTI_UNITS_MM;                     // units
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
        if(fwrite((void *)ImgM, sizeof(float), V, fp) != V)
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
        memset(extension, 0, 4);                    // for now, ignore the fact that we may have an extended header
        NiiHdr.dim[0] = 3;                          // 3d volume
        NiiHdr.dim[4] = 0;
        NiiHdr.toffset = 0.;                        // time offset shift
        NiiHdr.slice_duration = 1.;
        NiiHdr.slice_code = NIFTI_SLICE_UNKNOWN;
        NiiHdr.xyzt_units = NIFTI_UNITS_MM;                     // units
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
        if(fwrite((void *)ImgV, sizeof(float), V, fp) != V)
            Cleanup("can't write data image file");
        fclose(fp);
    }           // for(n...

    if(vflg == TRUE) {
        printf("'SAM3dc' done\n");
        fflush(stdout);
    }
    exit(0);
}
