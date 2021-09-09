// OPMsim -- simulate solution of cortical source activity &
//  independent source resolution, given an array of OPM
//  sensors & an array of test dipole locations. Thus simulation
//  employs the forward solution for a dipole in a homogeneous
//  spherical conductor (Sarvas).
//
//  Note: The average skull thickness for adult males is 6.5 mm
//      & for adult females is 7.1 mm. CSF will double this to 
//      about 15mm.  
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>
#include <siglib.h>
#include <samlib.h>
#include <DataFiles.h>
#include <samutil.h>
#include <filters.h>

#define HEAD_RADIUS     0.09        // 9 cm head radius
#define DURATION        300.0       // 5 minute duration
#define SAMPLERATE      350.0       // 350 Hz sample rate
#define BANDWIDTH       175.0       // 175 Hz OPM bandwidth
#define NOTCHFREQ       60.0        // 60 Hz & harmonic notch filter
#define NOISE           1.0e-14     // 10 fT/√Hz nominal OPM noise density
#define DIPMOMENT       2.5e-08     // default 25 nA-m rms dipole moment
#define OUTERRAD        0.08        // 8 cm outer radius
#define INNERRAD        0.06        // 6 cm inner radius
#define IOTA            1.0e-09     // small distance for loop completion
#define SIM_VERSION     "v4.1"


int
main(
    int             argc,
    char            **argv
)
{
    // variables
    HeaderInfo  DataHeader;     // measured data header
    ChannelInfo *DataChan;      // measured channel info
    HeaderInfo  SimHeader;      // simulation data header
    ChannelInfo *SimChan;       // SimChan[M] -- channel information data 
    ChannelInfo *OrigChan;   // SimChanOrig[M] -- channel information data
    EpochInfo   *Epoch;         // MEG epoch info
    FFTSPEC     DataFFT;        // data filter
    FFTSPEC     NotchFFT;       // notch filter
    FIELD       Q6;             // 6-parameter test dipole
    FIELD       Q5;             // 5-parameter test dipole
    VOXELINFO   *Dipole;        // Dipole[D] -- array of test dipoles
    double      rot;      	// angle between North pole and sensor orientation vector
    gsl_matrix  *A;             // copy of C (for LU inversion)
    gsl_matrix  *C;             // covariance matrix - MxM
    gsl_matrix  *Cinv;          // inverse covariance matrix - MxM
    gsl_vector  *Bp;            // forward solution - Mx1
    gsl_vector  *W;             // beamformer weights - Mx1
    gsl_matrix  *V;             // for SVD
    gsl_vector  *Sing;          // singular values
    gsl_vector  *work;          // for SVD
    gsl_permutation *p;         // LU decomposition permutation
    double      rand_rad3;      // Random radius between INNERRAD^3 and OUTERRAD^3
    double      **data;         // data[M][T] -- MEG data time-series
    double      *In;            // In[T] -- input epoch time-series
    double      *Out;           // Out[T] -- output epoch time-series
    double      **dip;          // dip[D][T] -- dipole source time series
    double      **src;          // src[D][T] -- SAM time-series
    double      **forward;      // forward[D][M] -- dipole forward solution matrix
    double      *Brain;         // Brain[T] -- simulated brain noise time-series
    double      *BrainScaled;   // Brain[T] -- scaled simulated brain noise time-series
    double      *Gain;          // Gain[D] -- recovered waverform gain -- relative to simulated waveform
    double      **distance;     // distance -- distances between sensors or dipoles
    double      **correlation;  // correlation[D][D] -- correlation between dipole & source time-series
    double      **crosstalk;    // crosstalk[M][M] -- crosstalk matrix
    double      Vc[3];          // cartesian vector
    double      Vs[3];          // spherical vector
    double      Vb[3];          // baseline vector
    double      Vx[3];          // x' tangential vector
    double      Vy[3];          // y' tangential vector
    double      Voff[3];        // orientation offset vector (spherical coords)
    double      Voffc[3];       // orientation offset vector (cartesian coords)
    double      u[3];       	// cross product of North pole and sensor orientation vector  
    double      Vnew[3];        // new sensor orientation vector 
    double      MaxSensor;      // maximum z-coordinate of sensor array
    double      MaxDipole;      // maximum z-coordinate of dipole array
    double      Dx, Dy;         // test dipole coordinates
    double      interval;       // dipole spacing
    double      NoiseDensity = NOISE;   // OPM sensor noise density (T/√Hz)
    double      OPMNoise;       // rms OPM sensor noise
    double      BrainMax;       // maximum brain noise (w/o bandwidth)
    double      BrainMin;       // minimum brain noise (w/o bandwidth)
    double      MaxBrainNoise;  // rms maximum brain noise
    double      MinBrainNoise;  // rms minimum brain noise
    //double      BrainScale;   // peakiness of Levy distribution
    double      BrainPow;       // peakiness of Levy distribution
    double      BrainNoise;     // just more noise, as far as I'm concerned...
    double      HPfreq = 0.;    // band-reject highpass frequency
    double      LPfreq = 0.;    // band-reject lowpass frequency
    double      BandWidth = BANDWIDTH;  // filter or simulation bandwidth
    double      Noise2;         // sensor noise power
    double      GainError;      // rms gain error
    double      Src;            // source power
    double      Nse;            // noise power
    double      span;           // span of eigensystem solution
    double      sx;             // sum of x
    double      sx2;            // sum of x^2
    double      sy;             // sum of y
    double      sxy;            // sum of x*y
    double      isr;            // independent source resolution
    double      r, r2, rvar;
    double      min;            // minimum distance between sensors
    double      rms;            // rms
    double      angle;          // angle (degrees)
    double      MomVec;         // tangential moment vector (radians)
    double      Moment = DIPMOMENT; // default dipole moment
    double      ArcOffset = 0.; // sensor vector offset
    double      Muddle = 0.;    // muddle the positions of the sensors
    double      CorticalDepth = 0.; // depth that sources are placed into the noise shell 
    double      percent;        // percent done
    double      edge;           // source array edge
    double      gap = 0.;       // array offset
    double      xtalk;          // crosstalk to nearest neighbor
    double      scale;          // crosstalk scale factor
    double      range;          // image range (edge + 10 mm)
    double      tmp;            // accumulator of sums
    int         NumNoise;       // number of noise sources
    int         FilterType;     // filter
    int         D;              // number of dipoles
    int         M;              // number of primary sensors
    int         R;              // number of reference sensors
    int         S;              // total sensors M+R
    int         T;              // number of samples
    int         c;              // coil index
    int         d;              // dipole index
    int         dd;             // alternate dipole index
    int         i;              // useful index
    int         j;              // 'nuther
    int         m;              // channel index
    int         mm;             // alternate channel index
    int         t;              // sample index
    int         v;              // vector index
    int         signum;         // sign of LU decomposition
    int         pct;            // integer percentage done
    int         LongIndex;      // position of arguments
    int         Modulus = 1;    // modulus for number of brain noise time series
    int         seed = 0;       // seed for random number generation for noise dipoles
    int         aflg = FALSE;   // arc-offset flag
    int         bflg = FALSE;   // recorded brain noise flag
    int         fflg = FALSE;   // fixed source grid flag
    int         nflg = FALSE;   // optionally naming flag
    int         pflg = FALSE;   // gain error (percentage) flag
    int         rflg = FALSE;   // random moment vector flag
    int         sflg = FALSE;   // simulation file name flag
    int         vflg = FALSE;   // verbose flag
    int         xflg = FALSE;   // crosstalk flag
    int         dflg = FALSE;   // flag for specified seed
    int         eflg = FALSE;   // error flag
    int		mflg = FALSE;   // muddle the sensor positions
    char        ArrayName[32];  // name of OPM sensor array file
    char        OutName[64];    // optional output file label
    char        SetName[64];    // dataset file namme
    char        SimName[64];    // simulation file name
    char        Prefix[128];     // filename prefix
    char        fpath[128];     // file path name
    char        Line[256];      // line buffer
    char        Label[32];      // parameter label
    char        State[8];       // state -- TRUE | FALSE
    FILE        *fp;            // generic file pointer
    FILE        *fp2;           // generic file pointer number 2
    FILE        *flog;          // log file pointer
    extern char *optarg;
    extern int  opterr;
    char        ShortOpts[] = "b:p:s:e:n:m:t:a:o:g:c:u:frx:d:z:v";
    static struct option    LongOpts[] = {                  // command line options
        {"OPM_noise_file", required_argument, NULL, 'b'},   // read recorder brain noise
        {"parameter_file", required_argument, NULL, 'p'},   // read brain noise simulation file
        {"sensor_name", required_argument, NULL, 's'},      // name of file with sensor coordinates
        {"array_edge", required_argument, NULL, 'e'},       // source array edge (mm)
        {"number sources", required_argument, NULL, 'n'},   // number of dipole sources
        {"dipole_moment", required_argument, NULL, 'm'},    // dipole moment (nA-m)
        {"sensoe_noise", required_argument, NULL, 't'},     // optional sensor noise (can be over-ridden by parameter file)
        {"gap", required_argument, NULL, 'g'},              // gap between sensors and scalp
        {"output_name", required_argument, NULL, 'z'},      // output file prefix
        {"cortical_depth", required_argument, NULL, 'c'},   // cortical depth
        {"pct_gain_error", required_argument, NULL, 'a'},   // percent gain error
        {"vector_offset", required_argument, NULL, 'o'},    // sensor vector [origin] offset
        {"fixed_dipoles", no_argument, NULL, 'f'},          // fixed dipole grid
        {"random", no_argument, NULL, 'r'},                 // random source moment vector
        {"crosstalk", required_argument, NULL, 'x'},        // include crosstalk
        {"verbose", no_argument, NULL, 'v'},
        {"seed", required_argument, NULL, 'd'},
        {"muddle", required_argument, NULL, 'u'},	    // muddle the sensor positions
        {NULL, no_argument, NULL, 0}
    };

    // parse command line parameters
    opterr = 1;                             // enable error reporting
    while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
        switch (c) {
            case 'b':       // read recorded OPM brain noise
                sprintf(SetName, "%s", optarg);
                bflg = TRUE;
                break;
            case 'p':       // read simulation parameters
                sprintf(SimName, "%s", optarg);
                sflg = TRUE;
                break;
            case 's':           // sensor file name
                sprintf(ArrayName, "%s", optarg);
                sflg = TRUE;
                break;
            case 'z':          // output file prefix
                sprintf(OutName, "%s", optarg);
                nflg = TRUE;
                break;
            case 'e':       // area edge
                edge = 0.001 * atof(optarg);
                break;
            case 'n':           // number of sources
                D = atoi(optarg);
                break;
            case 'm':       // dipole moment (optional)
                Moment = 1.0e-09 * atof(optarg);
                break;
            case 't':       // sensor noise (optional)
                NoiseDensity = 1.0e-15 * atof(optarg);
                break;
            case 'g':       // additional gap between sensors and scalp (optional)
                gap = .001*atof(optarg);
                break;
            case 'a':       // percent gain error
                GainError = 0.01 * atof(optarg);
                pflg = TRUE;
                break;
            case 'o':       // vector offset degrees X & Y
                ArcOffset = atof(optarg);
                if (ArcOffset > 0)
                    aflg = TRUE;
                break;
            case 'c':       // Cortical Depth 
                CorticalDepth = atof(optarg);
                CorticalDepth = .001*CorticalDepth;
                break;
            case 'u':       // muddle the positions 
                Muddle = atof(optarg);
                Muddle = .001*Muddle;
                if (Muddle > 0)
                    mflg = TRUE;
                break;
            case 'f':       // fixed dipole grid -- otherwise random positions
                fflg = TRUE;
                break;
            case 'r':       // generate random source moment vectors -- else parallel moments
                rflg = TRUE;
                break;
            case 'x':       // include crosstalk matrix
                xtalk = atof(optarg);
                xflg = TRUE;
                break;
            case 'v':       // verbose output
                vflg = TRUE;
                break;
            case 'd':       // use a specified seed instead of system time
                seed = atoi(optarg);
                dflg = TRUE;
                break;
            default:
                eflg = TRUE;
            break;
        }
    }

    if(eflg || !sflg) {
        fprintf(stderr, "OPMsim [arguments]\t%s\n", SIM_VERSION);
        fprintf(stderr, "\t-b <OPM brain noise dataset> (optional)\n");
        fprintf(stderr, "\t-p <simulation parameter file name>\n");
        fprintf(stderr, "\t-t <sensor noise density fT/√Hz (optional)\n");
        fprintf(stderr, "\t-v -- verbose\n");
        exit(-1);
    }

    // announce program
    if (vflg) {
        printf("OPMsim version %s\n", SIM_VERSION);
        fflush(stdout);
    }

    // read simulation parameters
    if (vflg) {
        printf("parsing brain noise simulation parameters");
        fflush(stdout);
    }
    
    sprintf(fpath, "%s", SimName);
    if ((fp = fopen(fpath, "r")) == NULL)
        Cleanup("can't open simulation parameter file");
    while (fgets(Line, 128, fp) != NULL) {              // scan parameter values
        if (strlen(Line) == 1)                          // skip empty lines
            continue;
        if (Line[0] == '#')                             // skip comment lines
            continue;
        if (sscanf(Line, "%s", Label) != 1)             // read the label field
            Cleanup("can't read label field");
        if (!strncmp("NumNoise", Label, 8)) {           // read 'NumNoise'
            if (sscanf(Line, "%s%d", Label, &NumNoise) != 2)
                Cleanup("can't read 'NumNoise'");
        } else if (!strncmp("ArrayName", Label, 9)) {   // read 'ArrayName'
            if (sscanf(Line, "%s%s", Label, ArrayName) != 2)
                Cleanup("can't read 'ArrayName'");
        } else if (!strncmp("SourceEdge", Label, 10)) { // read 'SourceEdge'
            if (sscanf(Line, "%s%lf", Label, &edge) != 2)
                Cleanup("can't read 'SourceEdge'");
            edge *= 0.001;                              // convert mm to meters
        } else if (!strncmp("NumDipoles", Label, 10)) { // read 'NumDipoles'
            if (sscanf(Line, "%s%d", Label, &D) != 2)
                Cleanup("can't read 'ArrayName'");
        } else if (!strncmp("DipMoment", Label, 9)) {   // read 'DipMoment'
            if (sscanf(Line, "%s%lf", Label, &Moment) != 2)
                Cleanup("can't read 'DipMoment'");
            Moment *= 1.0e-09;                          // convert nA-m to A-m
        } else if (!strncmp("OPMOffset", Label, 9)) {   // read 'OPMOffset'
            if (sscanf(Line, "%s%lf", Label, &gap) != 2)
                Cleanup("can't read 'OPMOffset'");
            gap *= 0.001;                         // convert mm to meters
        } else if (!strncmp("CorticalDepth", Label, 9)) {   // read 'OPMOffset'
            if (sscanf(Line, "%s%lf", Label, &CorticalDepth) != 2)
                Cleanup("can't read 'CorticalDepth'");
        } else if (!strncmp("FixedGrid", Label, 9)) {   // read 'FixedGrid'
            if (sscanf(Line, "%s%s", Label, State) != 2)
                Cleanup("can't read 'FixedGrid'");
            if (!strncmp("TRUE", State, 4))
                fflg = TRUE;
            else
                fflg = FALSE;
        } else if (!strncmp("RandomOri", Label, 9)) {   // read 'RandomOri'
            if (sscanf(Line, "%s%s", Label, State) != 2)
                Cleanup("can't read 'RandomOri'");
            if (!strncmp("TRUE", State, 4))
                rflg = TRUE;
            else
                rflg = FALSE;
        } else if (!strncmp("GainError", Label, 9)) {   // read 'GainError'
            if (sscanf(Line, "%s%lf", Label, &GainError) != 2)
                Cleanup("can't read 'GainError'");
            GainError *= 0.01;                          // convert percent to fraction
            pflg = TRUE;
        } else if (!strncmp("CrossTalk", Label, 9)) {   // read 'CrossTalk'
            if (sscanf(Line, "%s%lf", Label, &xtalk) != 2)
                Cleanup("can't read 'CrossTalk'");
            xtalk *= 0.01;                              // convert percent to fraction
            if (xtalk !=0) 
                xflg = TRUE;
        } else if (!strncmp("ArcOffset", Label, 9)) {   // read 'ArcOffset'
            if (sscanf(Line, "%s%lf", Label, &ArcOffset) != 2)
                Cleanup("can't read 'ArcOffset'");
            if (ArcOffset > 0)
                aflg = TRUE;
        } else if (!strncmp("OPMNoise", Label, 8)) {    // read 'OPMNoise' (fT/√Hz)
            if (sscanf(Line, "%s%lf", Label, &NoiseDensity) != 2)
                Cleanup("can't read 'NoiseDensity'");
            NoiseDensity *= 1.0e-15;
        } else if (!strncmp("BrainMax", Label, 8)) {    // read 'BrainMax'
            if (sscanf(Line, "%s%le", Label, &BrainMax) != 2)
                Cleanup("can't read 'BrainMax'");
        } else if (!strncmp("BrainMin", Label, 8)) {    // read 'BrainMin'
            if (sscanf(Line, "%s%le", Label, &BrainMin) != 2)
                Cleanup("can't read 'BrainMin'");
        } else if (!strncmp("Modulus", Label, 7)) {     // read 'Modulus'
            if (sscanf(Line, "%s%d", Label, &Modulus) != 2)
                Cleanup("can't read 'Modulus'");
        //} else if (!strncmp("BrainScale", Label, 10)) { // read 'BrainScale'
        //    if (sscanf(Line, "%s%lf", Label, &BrainScale) != 2)
        //        Cleanup("can't read 'BrainScale'");
        } else if (!strncmp("BrainPow", Label, 10)) { // read 'BrainPow'
            if (sscanf(Line, "%s%lf", Label, &BrainPow) != 2)
                Cleanup("can't read 'BrainPow'");
        } else if (!strncmp("Filter", Label, 10)) { // read 'Filter'
            if (sscanf(Line, "%s%lf%lf", Label, &HPfreq, &LPfreq) != 3)
                Cleanup("can't read 'Filter' frequencies");
        } else {
            fprintf(stderr, "\nBad label: '%s'\n", Label);
            fflush(stdout);
        }
    }
    fclose(fp);

    // form prefix for either simulation or brain noise measurement
    if (bflg) {
        sprintf(Prefix, "%s_%s_OPMMEG", SimName, SetName);
        nflg = FALSE;
    } else if (nflg) {
        sprintf(Prefix, "%s_%s_OPMSIM", SimName, OutName);
    } else {
        sprintf(Prefix, "%s_OPMSIM", SimName);
    }

    // open log file
    sprintf(fpath, "%s_Log.txt", Prefix);
    if((flog = fopen(fpath, "w")) == NULL)
        Cleanup("can't open log for write");

    // if reading brain noise, open dataset
    if (bflg) {

        // open OPM dataset
        if (vflg) {
            printf(" - done\n");
            printf("opening OPM brain noise dataset");
            fflush(stdout);
        }
        GetDsInfo(SetName, &DataHeader, &DataChan, &Epoch, NULL, TRUE);
        M = DataHeader.NumPri;          // number of primary sensors in recorded data
        T = DataHeader.MaxSamples;
        if (DataHeader.NumEpochs > 1)
            Cleanup("single epoch dataset req'd");

        // create highpass & notch filters
        if (vflg) {
            printf(" - done\n");
            printf("making %g to %g Hz bandpass & %g Hz FFT notch filters", HPfreq, LPfreq, NOTCHFREQ);
            fflush(stdout);
        }
        if (HPfreq == 0. && LPfreq == 0.)
            Cleanup("'Filter' req'd");
        if (HPfreq == 0. && LPfreq > 0.)
            FilterType = LOWPASS;
        else if (HPfreq > 0. && LPfreq > 0.)
            FilterType = BANDPASS;
        else if (HPfreq > 0. && LPfreq ==0.)
            FilterType = HIGHPASS;
        Butterworth(&DataFFT, T, HPfreq, LPfreq, DataHeader.SampleRate, FilterType, MAXORDER, &BandWidth);
        mknotch(&NotchFFT, HPfreq, LPfreq, DataHeader.SampleRate, T, NOTCHFREQ);

        // allocate memory
        data = new_matrixE(M, T, "data[M][T]");
        In = new_arrayE(double, T, "In[T]");
        Out = new_arrayE(double, T, "Out[T]");

        // read & filter data
        if (vflg) {
            printf(" - done\n");
            printf("reading data");
            fflush(stdout);
        }
        for (m=0; m<M; m++) {
            mm = DataHeader.PsIndex[m];
            scale = DataChan[mm].Scale;
            GetDsData(&DataHeader, 0, mm, scale, In);
            demean(In, T);
            FFTfilter(In, Out, &DataFFT, T);
            FFTfilter(Out, In, &NotchFFT, T);
            for (t=0; t<T; t++)
                data[m][t] = In[t];

            // channel heartbeat
            if(vflg == TRUE) {
                printf(".");
                fflush(stdout);
            }
        }
        free(In);
        free(Out);

        // compute total measured noise
        for (m=0, tmp=0.; m<M; m++)
            for (t=0; t<T; t++)
                tmp += data[m][t] * data[m][t];
        rms = sqrt(tmp / (double)(M * T));
        fprintf(flog, "measured brain noise: %6.1f fT rms\n", 1.0e+15*rms);
        fprintf(flog, "measured brain noise density: %5.1f fT/rt-Hz\n", 1.0e+15*rms/sqrt(BandWidth));
        fflush(flog);

    } // end reading real brain noise

    // read sensor file & allocate 'Channel'
    if (vflg) {
        printf(" - done\n");
        printf("generating 'SimChan' array from sensor file: '%s'", ArrayName);
        fflush(stdout);
    }
    if((fp = fopen(ArrayName, "r")) == NULL)
        Cleanup("can't open sensor file for read");
    if(fscanf(fp, "%d%d", &m, &R) != 2)         // list of channels begins with primary sensors, then reference sensors
        Cleanup("can't read number of sensors");

    // sanity check
    if (bflg && m != M) {
        fprintf(stderr, "\nnumber of sensors (%d) in OPM recorded data does not agree with that in sensor file (%d)", m, M);
        Cleanup("can't continue");
    }
    M = m;
    S = M + R;

    // create sensors -- sensors are in the order of: ADC, primary sensors, reference sensors (same as OPM datasets)
    if (vflg) {
        printf(" - done\n");
        printf("populating 'SimChan' for %d primary & %d reference OPM sensors", M, R);
        fflush(stdout);
    }
    SimChan = new_arrayE(ChannelInfo, S+1, "SimChan[S]");
    OrigChan = new_arrayE(ChannelInfo, S+1, "OrigChan[S]");

    // make 1st channel the ADC
    SimChan[0].ChannelType = TYPE_EXTERNAL;
    SimChan[0].Geom.NoGeom = TRUE;
    sprintf(SimChan[0].ChannelName, "ADC");
    SimChan[0].ChannelNumber = 0;
    SimChan[0].Units = UNIT_VOLTS;
    SimChan[0].Scale = 1.;

    if (dflg) {
        printf("\nRandom number generator initialized with Seed %d \n",seed);
        srand48(seed);   // seed random number generator with user specified seed
    } else {
        srand48((long)time((time_t *)0));   // seed random number generator with system clock
    }

    // offset sensors (for simulating effect of gap between sensors & scalp)
    if (vflg && gap > 0.) {
        printf(" - done\n");
        printf("adding %4.1f mm offset from scalp", 1000.*gap);
        fflush(stdout);
    }

    // read primary & reference sensors
    for (m=1, MaxSensor=0.; m<(S+1); m++) {             // 1st channel is ADC

        // fill in 'Channel' elements
        if (m < M+1) {                                  // primary sensors are assigned radial orientations

            SimChan[m].ChannelType = TYPE_MEG;
            if (fscanf(fp, "%s%lf%lf%lf",               // read OPM primary sensor position, only (mm) & name
              SimChan[m].ChannelName, &Vc[X_], &Vc[Y_], &Vc[Z_]) != 4) {
                printf("%s%lf%lf%lf",SimChan[m].ChannelName, Vc[X_], Vc[Y_], Vc[Z_]);
                Cleanup("can't read primary sensor position");
            }
          
#ifdef EBUG
            printf("\nm=%d, name=%s\n", m, SimChan[m].ChannelName);
            printf("\t%f %f %f\n", Vc[X_], Vc[Y_], Vc[Z_]);
            fflush(stdout);
#endif
            for (v=X_, r2=0.; v<=Z_; v++)
                r2 += Vc[v] * Vc[v];
            r = sqrt(r2);
            for (v=X_; v<=Z_; v++)                      // Vb is an outward pointing unit vector from array origin (assuming common origin)
                Vb[v] = Vc[v] / r;
            for (v=X_; v<=Z_; v++) {
                SimChan[m].Geom.OPMSensor.origin.p[v] = 0.001 * Vc[v];
                SimChan[m].Geom.OPMSensor.origin.v[v] = Vb[v];
                SimChan[m].Geom.OPMSensor.LocalSphere[v] = 0.;
            }

            // add gap, if requested
            if (gap > 0.)
                SimChan[m].Geom.OPMSensor.origin.p[Z_] += gap;


            // find maximum z-coordinate of sensor array
            if (SimChan[m].Geom.OPMSensor.origin.p[Z_] > MaxSensor)
                MaxSensor = SimChan[m].Geom.OPMSensor.origin.p[Z_];

        } else {                                        // create references sensors

            SimChan[m].ChannelType = TYPE_REF_MAG;
            if (fscanf(fp, "%s%lf%lf%lf%lf%lf%lf",      // read OPM reference sensor parameters & name
              SimChan[m].ChannelName, &Vc[X_], &Vc[Y_], &Vc[Z_], &Vb[X_], &Vb[Y_], &Vb[Z_]) != 7)
                Cleanup("can't read reference sensor parameters");
#ifdef EBUG
            fflush(stdout);
#endif
        }
        SimChan[m].Flag = TRUE;
        SimChan[m].ChannelNumber = m;
        SimChan[m].Units = UNIT_TESLA;
        SimChan[m].Scale = 1.;
        SimChan[m].AnalogHpFreq = 0.;
        SimChan[m].AnalogLpFreq = BandWidth;
        SimChan[m].HeadR2 = HEAD_RADIUS * HEAD_RADIUS;
        SimChan[m].BalOrder = 0;
        SimChan[m].Balance = NULL;
        SimChan[m].NumBal = 0;
    }
 
    fclose(fp);

    // before perturbing angles, copy the SimChan to an original version
    memcpy(OrigChan, SimChan, sizeof(ChannelInfo) * (S+1)); 
    

    // optionally, offset sensor vector by adding a fixed angular offset in a random direction in the tangential plane. 
    
    if (aflg & ArcOffset > 0) {

       sprintf(fpath, "%s_Orientations.txt", Prefix);
       if ((fp2 = fopen(fpath,"w")) == NULL)
          Cleanup("can't open orientation angle file");

        if (vflg) {
            printf("perturbing vectors by up to %4.1f degrees\n", ArcOffset);
            fflush(stdout);
        }

        // read primary & reference sensors
        for (m=1; m<(S+1); m++) {             // 1st channel is ADC

            // fill in 'Channel' elements
            if (m < M+1) {                                  // primary sensors are assigned radial orientations

            for (v=X_; v<=Z_; v++) {
                Vb[v] = SimChan[m].Geom.OPMSensor.origin.v[v];
            }

	    Voff[THETA_] = drand48() * 2 * M_PI;
            Voff[PHI_] = ArcOffset / 180 * M_PI;
            Voff[RHO_] = 1;

            // convert to cartesian coordinates
            StoC(Voff, Voffc);

	    // find cross product of [0,0,1] with orientation vector to get rotation axis
            u[X_] = -1*Vb[Y_];  // Vb is the original orientation
            u[Y_] = 1*Vb[X_];
            u[Z_] = 0;   	    // don't really need, but...
            for (v=X_, r2=0.; v<=Z_; v++)   // normalize u
                r2 += u[v] * u[v];
            r = sqrt(r2);
            for (v=X_; v<=Z_; v++)
                u[v] = u[v] / r;
	    // product of rotation axis with orientation vector gives rotation angle
            rot = acos(Vb[Z_]);
            // multiply the rotation matrix by the cartesian vector Voffc to get Vnew
            Vnew[X_] = (cos(rot) + u[X_]*u[X_]*(1 - cos(rot))) * Voffc[X_] + (u[X_]*u[Y_]*(1 - cos(rot))) * Voffc[Y_] + (u[Y_] * sin(rot)) * Voffc[Z_];
            Vnew[Y_] = (u[Y_]*u[X_]*(1 - cos(rot))) * Voffc[X_] + (cos(rot) + u[Y_]*u[Y_]*(1 - cos(rot))) * Voffc[Y_] + (-u[X_] * sin(rot)) * Voffc[Z_];
            Vnew[Z_] = (-u[Y_]*sin(rot)) * Voffc[X_] + (u[X_] * sin(rot)) * Voffc[Y_] + cos(rot) * Voffc[Z_];
            for (v=X_; v<=Z_; v++)
                SimChan[m].Geom.OPMSensor.origin.v[v] = Vnew[v];
                fprintf(fp2, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",Vc[X_],Vc[Y_],Vc[Z_],Vb[X_],Vb[Y_],Vb[Z_],Vnew[X_],Vnew[Y_],Vnew[Z_]);
            }
       }
       fclose(fp2);
    }

    // populate simulation header
    if (vflg) {
        printf("\npopulating HeaderInfo & initializing constants");
        fflush(stdout);
    }
 
    SimHeader.Version = 1;
    SimHeader.NumChannels = S;
    SimHeader.NumAcquired = S;
    SimHeader.NumPri = M;
    SimHeader.NumRef = R;
    SimHeader.SqIndex = new_arrayE(int, S, "SqIndex[S]");
    for (m=0; m<S; m++)
        SimHeader.SqIndex[m] = m + 1;
    SimHeader.PsIndex = new_arrayE(int, M, "PsIndex[M]");
    for (m=0; m<M; m++)
        SimHeader.PsIndex[m] = m + 1;
    SimHeader.RsIndex = new_arrayE(int, R, "RsIndex[R]");
    for (m=0; m<R; m++)
        SimHeader.RsIndex[m] = m + M;
    SimHeader.NumEpochs = 1;
    SimHeader.MaxSamples = (int)(SAMPLERATE * DURATION);
    SimHeader.SampleRate = SAMPLERATE;
    SimHeader.PowerFreq = 60.;
    SimHeader.DataFormat = 0;

    if (dflg) {
       srand48(seed); // reset the seed here so the random positions are consistent across runs 
    }
    
    // optionally, muddle the positions of the sensors as well
    if (mflg & Muddle > 0) {

       if (vflg) {
           printf("Muddling the sensor positions by up to %f mm, then writing ideal and muddled positions",1000*Muddle);
           fflush(stdout);
       }

       sprintf(fpath, "%s_PositionOffsets.txt", Prefix);
       if ((fp2 = fopen(fpath,"w")) == NULL)
          Cleanup("can't open PositionOffset file");
        
        // read primary & reference sensors
        for (m=1; m<(S+1); m++) {             // 1st channel is ADC

            // fill in 'Channel' elements
            if (m < M+1) {                                  // primary sensors are assigned radial orientations

                for (v=X_; v<=Z_; v++) { 
                    Vc[v] = SimChan[m].Geom.OPMSensor.origin.p[v];
		    Voffc[v] = Vc[v] + (Muddle * (2 *drand48() - 1));  
                    SimChan[m].Geom.OPMSensor.origin.p[v] = Voffc[v]; 
                } 
                fprintf(fp2, "%f\t%f\t%f\t%f\t%f\t%f\n",Vc[X_],Vc[Y_],Vc[Z_],Voffc[X_],Voffc[Y_],Voffc[Z_]);

            }
        }
        fclose(fp2);

    } else { 
    
    // output sensor positions to file
    
        if (vflg) {
            printf(" - done\n");
            printf("writing sensor positions");
            fflush(stdout);
        }
        sprintf(fpath, "%s_Sensor%d.txt", Prefix, M);
        if((fp = fopen(fpath, "w")) == NULL)
            Cleanup("can't open 'Sensors' file for write");
        for (m=0; m<M; m++) {
            mm = SimHeader.PsIndex[m];
            fprintf(fp, "%7.3f\t%7.3f\t%7.3f\n",
              1000.*SimChan[mm].Geom.OPMSensor.origin.p[X_], 1000.*SimChan[mm].Geom.OPMSensor.origin.p[Y_], 1000.*SimChan[mm].Geom.OPMSensor.origin.p[Z_]);
        }
        fclose(fp);
    }

    // compute distances between sensors
    distance = new_matrixE(M, M, "distance[M][M]");
    for (m=0, min=HUGE; m<M; m++) {
        i = SimHeader.PsIndex[m];
        for (mm=0; mm<M; mm++) {
            j = SimHeader.PsIndex[mm];
            for (v=X_, r2=0.; v<=Z_; v++) {
                tmp = SimChan[i].Geom.OPMSensor.origin.p[v] - SimChan[j].Geom.OPMSensor.origin.p[v];
                r2 += tmp * tmp;
            }
            distance[m][mm] = sqrt(r2);
            if (distance[m][mm] > 0. && distance[m][mm] < min)
                min = distance[m][mm];
        }
    }

    // create crosstalk matrix
    if (xflg) {
        if (vflg) {
            printf(" - done\n");
            printf("creating crosstalk matrix -- %4.1f percent to nearest neighbor sensor", xtalk);
            fflush(stdout);
        }
        fprintf(flog, "crosstalk matrix: %4.1f percent to nearest sensor\n", xtalk);
        fflush(flog);
        scale = 0.01 * xtalk * pow(min, 3.);
        crosstalk = new_matrixE(M, M, "crosstalk[M][M]");
        for (m=0; m<M; m++)
            for (mm=0; mm<M; mm++) {
                if (m == mm)
                    crosstalk[m][mm] = 1.;
                else
                    crosstalk[m][mm] = scale / pow(distance[m][mm], 3.);
            }
        free_matrix(distance);
    }

    // simulate brain + sensor noise
    if (!bflg) {
   
        //initialize some constants
        T = SimHeader.MaxSamples;
        OPMNoise = NoiseDensity * sqrt(BandWidth);
        MaxBrainNoise = BrainMax * sqrt(BandWidth);
        MinBrainNoise = BrainMin * sqrt(BandWidth);
        Noise2 = OPMNoise * OPMNoise;       // noise power

        // allocate space for data (excluding the reference channels)
        if (vflg) {
            printf(" - done\n");
            printf("allocating data array");
            fflush(stdout);
        }
        data = new_matrixE(M, T, "data[M][T]");

        // simulate sensor noise
        if (vflg) {
            printf(" - done\n");
            printf("generating simulated sensor noise time-series");
            fflush(stdout);
        }
        if (dflg) {
 	    srand48(seed); // reset the seed here so the noise model is the same every time 
        }
        fprintf(flog, "simulating sensor noise: %5.1f fT/√Hz (%7.1f fT rms)\n", 1.0e+15*NoiseDensity, 1.0e+15*OPMNoise); 
        for (m=0; m<M; m++)         // for each sensor...
            for (t=0; t<T; t++)     // generate Gaussian random time-series with zero mean & standard deviation of OPMNoise
                data[m][t] = RandomGauss(0., OPMNoise);

#ifdef HISTO
        // file BrainNoise for histogram plot
        sprintf(fpath, "%s_BrainNoise.txt", Prefix);
        if ((fp = fopen(fpath, "w")) == NULL)
            Cleanup("can't open 'BrainNoise' for write");
#endif

        // simulate brain noise in 20 mm thick hemispheric shell
        if (vflg) {
            printf(" - done\n");
            printf("simulating brain noise time-series\n");
            fflush(stdout);
        }
        Bp = gsl_vector_alloc(M);
        Brain = new_arrayE(double, T, "Brain[T]");
        BrainScaled = new_arrayE(double, T, "BrainScaled[T]");

        if (dflg) {
 	    srand48(seed); // reset the seed here so the noise model is the same as that produced by SAMsimulate
        }

        for (i=0, pct=(-1); i<NumNoise; i++) {

            // generate a random unit dipole within the shell, in spherical coordinates

            Q5.p[THETA_] = 2. * M_PI * drand48();       // 0 <= theta <= 2*pi
            Q5.p[PHI_] = M_PI;                          // initially set phi in the Southern hemisphere
            while(Q5.p[PHI_] > M_PI / 2. )              // keep picking phi values until you get one in the Northern hemisphere
              Q5.p[PHI_] = acos(2. * drand48() -1.);    // northern hemisphere: 0 <= phi <= pi/2
            rand_rad3 = (pow(OUTERRAD,3) - pow(INNERRAD,3)) * drand48() + pow(INNERRAD,3);
            Q5.p[RHO_] = cbrt(rand_rad3);               // area of a sphere scales like r^3 - need equal dipoles per unit area!

            Q5.v[PSI_] = 2. * M_PI * drand48();             // 0 <= psi <= 2*pi -- tangential moment vector
            Q5.v[Q_] = 1.;                                  // |Q| = 1.0
            P5toP6(&Q5, &Q6);                               // convert this to 6-parameter cartesian dipole

            // simulate source field using Sarvas equation for unit (1.0 A-m dipole)
            for(m=0; m<M; m++) {
                mm = SimHeader.PsIndex[m];                  // make sure we are selecting primary sensors
                gsl_vector_set(Bp, m, FwdOPM(&Q6, &SimChan[mm].Geom.OPMSensor));
            }

            // for each time sample compute brain noise time-series -- this is the tricky part!
            //do
            //    BrainNoise = MaxBrainNoise * Levy((20. * drand48() + 1.1e-12), 1.0e-12, BrainScale);
            //while (BrainNoise < MinBrainNoise);
   
            BrainNoise = (pow((double)i / (double)NumNoise, BrainPow) * (MaxBrainNoise-MinBrainNoise)) + MinBrainNoise;

            // generate Gaussian random time-series for this noise dipole
            if (i % Modulus == 0)
                for (t=0; t<T; t++)
                    Brain[t] = RandomGauss(0., 1);
#ifdef HISTO
                fprintf(fp, "%d\t%e\n", i, BrainNoise);
#endif
            for (t=0; t<T; t++)
		        BrainScaled[t] = BrainNoise * Brain[t];

            // sum this dipole MEG signal
            for (m=0; m<M; m++)
                for(t=0; t<T; t++)
                    data[m][t] += BrainScaled[t] * gsl_vector_get(Bp, m);

            // show progress (so I can be bored to tears)
            if (vflg) {
                percent = 100. * (double)i / (double)NumNoise;
                if (rint(percent) > (double)pct) {
                    pct = (int)rint(percent);
                    printf("\r\t%d percent done", pct);
                    fflush(stdout);
                }
            }
        }
#ifdef HISTO
        fclose(fp);
#endif
        free(Brain);
        free(BrainScaled);
        gsl_vector_free(Bp);

        // compute rms sensor noise plus brain noise
        for (m=0, tmp=0.; m<M; m++)
            for (t=0; t<T; t++)
                tmp += data[m][t] * data[m][t];
        rms = sqrt(tmp / (double)(M * T));
        fprintf(flog, "brain noise (%d sources): %6.1f fT rms\n", NumNoise, 1.0e+15*rms);
        fprintf(flog, "brain noise density: %5.1f fT/rt-Hz\n", 1.0e+15*rms/sqrt(BandWidth));
        fflush(flog);
    }

    // generate test dipoles
    if (vflg) {
        printf(" - done\n");
        printf("simulating %d", D);
        if (rflg)
            printf(" random");
        else
            printf(" parallel");
        printf(" dipoles in %4.1f x %4.1f area", 1000.*edge, 1000.*edge);
        fflush(stdout);
    }
    fprintf(flog, "simulated %d", D);
    if (rflg)
        fprintf(flog, " random");
    else
        fprintf(flog, " parallel");
    fprintf(flog, " dipoles in %4.1f x %4.1f area\n", 1000.*edge, 1000.*edge);
    fprintf(flog, "dipole moment: %4.1f nA-m\n", 1.0e+09*Moment);
    if (dflg) {
        srand48(seed); // reset the seed here so the test dipoles end up in the same place 
        //printf("\nre-initializing the seed to %d\n",seed);
    }
    fflush(flog);
    Dipole = new_arrayE(VOXELINFO, D, "Dipole[D]");
    if (fflg) {
        interval = edge / (sqrt((double)D) - 1.);
        Dx = -edge / 2.;
        Dy = -edge / 2.;
    }
    for (d=0, MaxDipole=0.; d<D; d++) {
        if(fflg) {
            Vc[X_] = Dx;
            Vc[Y_] = Dy;
            Vc[Z_] = OUTERRAD - CorticalDepth;
        } else {
            Vc[X_] = drand48() * edge - edge / 2.;
            Vc[Y_] = drand48() * edge - edge / 2.;
            Vc[Z_] = OUTERRAD - CorticalDepth;;
        }
        CtoS(Vc, Vs);                   // convert dipole coordinates to spherical coordinates
        Vs[RHO_] = OUTERRAD - CorticalDepth;
        StoC(Vs, Vc);
        for (v=X_, r2=0.; v<=Z_; v++)   // compute unit radial vector
        r2 += Vc[v] * Vc[v];
            r = sqrt(r2);
        for (v=X_; v<=Z_; v++)
            Vb[v] = Vc[v] / r;
        OrthoPlane(Vb, Vx, Vy);         // compute tangential unit vectors Vx & Vy
        if (rflg)
            MomVec = 2. * M_PI * drand48();
        else
            MomVec = M_PI_4;
        for (v=X_; v<=Z_; v++) {
            Dipole[d].p[v] = Vc[v];
            Dipole[d].v[v] = Vx[v] * cos(MomVec) + Vy[v] * sin(MomVec);
        }
        Dipole[d].Solve = TRUE;
        Dipole[d].ROI = TRUE;
        if (fflg) {
            Dx += interval;
            if(Dx > edge / 2.) {
                Dx = -edge / 2.;
                Dy += interval;
            }
        }

        // find maximum z-coordinate of test dipole array
        if (Dipole[d].p[Z_] > MaxDipole)
            MaxDipole = Dipole[d].p[Z_];
    }

    // output sensor & dipole positions to file
    if (vflg) {
        printf(" - done\n");
        printf("gap between sensor & dipole arrays is %4.1f mm\n", 1000.*(MaxSensor-MaxDipole));
        printf("writing dipole positions");
        fflush(stdout);
    }
    fprintf(flog, "gap between sensor & dipole arrays is %4.1f mm\n", 1000.*(MaxSensor-MaxDipole));
    fflush(flog);
    sprintf(fpath, "%s_Dipoles%d.txt", Prefix, D);
    if((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'Dipoles' file for write");
    for (d=0; d<D; d++)
        fprintf(fp, "%f\t%f\t%f\n", 1000.*Dipole[d].p[X_], 1000.*Dipole[d].p[Y_], 1000.*Dipole[d].p[Z_]);
    fclose(fp);

    // compute distances between dipole source pairs
    distance = new_matrixE(D, D, "distance[D][D]");
    for (d=0; d<D; d++)
        for (dd=0; dd<D; dd++) {
            for (v=X_, sx2=0.; v<=Z_; v++) {
                tmp = Dipole[d].p[v] - Dipole[dd].p[v];
                sx2 += tmp * tmp;
            }
            distance[d][dd] = sqrt(sx2);
        }

    if (dflg) {
       srand48(seed); // reset the seed here so the noise model is the same every time 
    }
    // simulate dipole time courses
    if (vflg) {
        printf(" - done\n");
        printf("simulating dipole waveforms");
        fflush(stdout);
    }
    dip = new_matrixE(D, T, "dip[D][T]");
    forward = new_matrixE(D, M, "forward[D][M]");
    for (d=0; d<D; d++) {

        // convert dipole parameters to FIELD
        for (v=X_; v<=Z_; v++) {
            Q6.p[v] = Dipole[d].p[v];
            Q6.v[v] = Dipole[d].v[v];
        }

        // generate random dipole time-series
        for (t=0; t<T; t++)
            dip[d][t] = RandomGauss(0., Moment);    // moment is A-m

        // simulate source using Sarvas equation for unit (1.0 A-m dipole)
        for (m=0; m<M; m++) {
            mm = SimHeader.PsIndex[m];
            forward[d][m] = FwdOPM(&Q6, &SimChan[mm].Geom.OPMSensor);
        }

        // sum the MEG signal for this dipole after scaling by moment -- including crosstalk
        if (xflg) {
            for (m=0; m<M; m++)
                for (mm=0; mm<M; mm++)
                    for (t=0; t<T; t++)
                        data[m][t] += dip[d][t] * crosstalk[m][mm] * forward[d][mm];
        } else {
            for (m=0; m<M; m++)
                for (t=0; t<T; t++)
                    data[m][t] += dip[d][t] * forward[d][m];
        }

        // heartbeat
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }
    if (dflg) {
       srand48(seed); // reset the seed here so the noise model is the same every time 
    }
    // if req'd, introduce random gain errors
    if (pflg) {
        if (vflg) {
            printf(" - done\n");
            printf("introducing +/-%5.3f%% random gain errors", 100.*GainError);
            fflush(stdout);
        }
        fprintf(flog, "+/-%5.3f%% random gain errors", 100.*GainError);
        fflush(flog);
        for (m=0; m<M; m++) {
            tmp = 1. + (2. * drand48() - 1.) * GainError;
            for (t=0; t<T; t++)
                data[m][t] *= tmp;
        }
    }

    // compute spatial angle matrix
    if (vflg) {
        printf(" - done\n");
        printf("computing spatial angle matrix");
        fflush(stdout);
    }
    sprintf(fpath, "%s_SpatialAngles.txt", Prefix);
    if((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open spatial angles file for write");
    for (d=0; d<D; d++) {           // make dipole forward solutions unit vectors
        for (m=0, r2=0.; m<M; m++)
            r2 += forward[d][m] * forward[d][m];
        r = sqrt(r2);
        for (m=0; m<M; m++)
            forward[d][m] /= r;
    }
    for (d=0; d<D; d++)             // compute dot product (& angle) between dipole forward solutions
        for (dd=0; dd<D; dd++) {
            for (m=0, sxy=0.; m<M; m++)
                sxy += forward[d][m] * forward[dd][m];
            angle = 180. / M_PI * acos(fabs(sxy));
            fprintf(fp, "angle[%d][%d]=%7.3f\n", d, dd, angle);
            fflush(fp);
        }
    fclose(fp);
    free_matrix(forward);

    // compute covariance
    if (vflg) {
        printf(" - done\n");
	printf("computing covariance matrix & its inverse");
	fflush(stdout);
    } 
    C = gsl_matrix_alloc(M, M);
    Cinv = gsl_matrix_alloc(M, M);
    A = gsl_matrix_alloc(M, M);
    p = gsl_permutation_alloc(M);
	for (i=0; i<M; i++)
		for (j=0; j<M; j++) {
			for (t=0, tmp=0.; t<T; t++)
			    tmp += data[i][t] * data[j][t];
            gsl_matrix_set(C, i, j, tmp / (double)T);
		}
    gsl_matrix_memcpy(A, C);    // make a copy of C
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_invert(A, p, Cinv);
    gsl_permutation_free(p);

    // compute eigenvalue spectrum using SVD
    if (vflg) {
	    printf(" - done\n");
	    printf("computing eigenvalue spectrum of covariance matrix");
	    fflush(stdout);
	}
    gsl_matrix_memcpy(A, C);    // make another copy of C
    V = gsl_matrix_alloc(M, M);
    Sing = gsl_vector_alloc(M);
    work = gsl_vector_alloc(M);
    gsl_linalg_SV_decomp(A, V, Sing, work);
    sprintf(fpath, "%s_Eigenvalues.txt", Prefix);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'eigenvalues' for write");
    gsl_vector_fprintf(fp, Sing, "%e");
//    for (i=0; i<M; i++)
//        fprintf(fp, "%e\n", gsl_vector_get(Sing, i));
    fclose(fp);
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(Sing);
    gsl_vector_free(work);

    // compute optimal weights for each dipole source location
    if (vflg) {
        printf(" - done\n");
        printf("computing weights & source correlation");
        fflush(stdout);
    }
    src = new_matrixE(D, T, "src[D][T]");
    Gain = new_arrayE(double, D, "Gain[D]");
    Bp = gsl_vector_alloc(M);
    W = gsl_vector_alloc(M);

    for (d=0; d<D; d++) {

        // set voxel coordinate & solve for optimal target vector
	// be sure to use the original unperturbed version of the channel array
        DIPSolve2D(&Dipole[d], Bp, Cinv, &SimHeader, OrigChan, 0, &span);
        SAMsolve(C, Cinv, Bp, W, Noise2, &Src, &Nse);

        // compute source time-series estimate
        for (t=0, sx=sx2=0; t<T; t++) {
            for (m=0, tmp=0.; m<M; m++)
                tmp += gsl_vector_get(W, m) * data[m][t];
            src[d][t] = tmp;
            sx += tmp;
            sx2 += tmp * tmp;
        }
        Gain[d] = sqrt((sx2 - sx * sx / (double)T) / (double)T) / Moment;   // ratio of rms recovered signal to rms simulated signal
    }

    // compute map of source strength (z^2) on 1 mm grid extending 10 mm beyond dipole grid
    if (vflg) {
        printf(" - done\n");
        printf("writing SAM image file");
        fflush(stdout);
    }
    range = edge / 2. + 0.005;
    tmp = 1000. * edge + 11.; tmp *= tmp;
    sprintf(fpath, "%s_Image.txt", Prefix);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open OPM_Image file for write");
    Vc[Z_] = OUTERRAD - CorticalDepth;;
    for (Vc[X_]=(-range), i=0, pct=(-1); Vc[X_]<=(range+IOTA); Vc[X_]+=0.001) {
        for (Vc[Y_]=(-range); Vc[Y_]<=(range+IOTA); Vc[Y_]+=0.001, i++) {
            CtoS(Vc, Vs);           // convert to spherical coordinates
            Vs[RHO_] = OUTERRAD - CorticalDepth;    // set image radius
            StoC(Vs, Vb);           // convert to cartesian coordinates
            for (v=X_; v<=Z_; v++) {
                Dipole[0].p[v] = Vb[v];
                Dipole[0].v[v] = 0.;
            }
            Dipole[0].Solve = TRUE;
            Dipole[0].ROI = TRUE;

            // compute source z^2 & save to file
            DIPSolve2D(&Dipole[0], Bp, Cinv, &SimHeader, OrigChan, 0, &span);
            SAMsolve(C, Cinv, Bp, W, Noise2, &Src, &Nse);
            fprintf(fp, "%f\t%f\t%f\t%f\t%f\n", 1000.*Vb[X_], 1000.*Vb[Y_], 1000.*Vb[Z_], 1.0e+18*Src, Src/Nse);
            fflush(fp);

            // show progress
            if (vflg) {
                percent = 100. * (double)i / tmp;
                if(rint(percent) > (double)pct) {
                    pct = (int)rint(percent);
                    printf("\r\t%d percent done", pct);
                    fflush(stdout);
                }
            }
        }
    }
    fclose(fp);

    // free memory
    gsl_matrix_free(C);
    gsl_matrix_free(Cinv);
    gsl_vector_free(Bp);
    gsl_vector_free(W);
    free(Dipole);

    // open & write gain file
    sprintf(fpath, "%s_Gain.txt", Prefix);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'OPM_Gain.txt' file for write");
    for (d=0; d<D; d++)
        fprintf(fp, "%8.6f\n", Gain[d]);
    fclose(fp);
    free(Gain);

    // open matrix file
    sprintf(fpath, "%s_Matrix.txt", Prefix);
    if((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'OPM_Matrix.txt' file for write");

    // compute correlation coefficient of SAM estimate with original dipole time-series
    if (vflg) {
        printf(" - done\n");
        printf("computing correlation of beamformer dipole estimates with original dipoles\n");
        fflush(stdout);
    }
    correlation = new_matrixE(D, D, "correlation[D][D]");
    for (d=0, pct=(-1); d<D; d++) {       // dipole index
        for (dd=0; dd<D; dd++) {        // source index
            correlation[d][dd] = fabs(gsl_stats_correlation(dip[d], 1, src[dd], 1, T));

            // tabulate results
            fprintf(fp, "%11.8f\n", correlation[d][dd]);
            fflush(stdout);

            // show progress
            if (vflg) {
                percent = 100. * (double)d / (double)D;
                if(rint(percent) > (double)pct) {
                    pct = (int)rint(percent);
                    printf("\r\t%d percent done", pct);
                    fflush(stdout);
                }
            }
        }
    }
    fclose(fp);

    // open & write isr file
    sprintf(fpath, "%s_ISR.txt", Prefix);
    if((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'OPM_ISR.txt' file for write");
    for (d=0; d<D; d++) {
        for (dd=0; dd<D; dd++) {
            fprintf(fp, "%9.6f\t%9.6f\n", 1000.*distance[d][dd], fabs(correlation[d][dd]));
        }
     }

    fclose(fp);

    free_matrix(src);
    free_matrix(dip);
    //
    // compute independent source resolution (mean diagonal - mean off-diagonal correlations)
    if (vflg) {
        printf(" - done\n");
        printf("computing independent source resolution");
        fflush(stdout);
    }
    for (d=0, sx=sy=0.; d<D; d++)
        for (dd=0; dd<D; dd++) {
            if (d == dd)
                sx += correlation[d][dd];
            else
                sy += correlation[d][dd];
        }
    isr = (sx / (double)D) - (sy / (double)(D * (D - 1)));
    rvar = 100. * (1. - isr * isr);
    if (vflg) {
        printf(" - done\n");
        printf("ISR = %f\n", isr);
        printf("variance = %7.3f percent\n", rvar);
    }
    free_matrix(correlation);

    // that's all, folks!
    fprintf(flog, "ISR = %f\n", isr);
    fprintf(flog, "residual variance = %f percent\n", rvar);
    fprintf(flog, "gap = %4.1f mm\n", 1000.*(MaxSensor-MaxDipole));
    fclose(flog);
    if (vflg) {
        printf(" - done\n\n");
        printf("OPM simulation done\n");
        fflush(stdout);
    }
    exit(0);
}
