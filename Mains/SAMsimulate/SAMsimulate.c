// SAMsimulate -- MEG simulation utility
//
//  Superposes a single dipole source of specified coordinate & magnitude
//  on recorded brain noise (task-free dataset), using a single sphere
//  model with the origin at the CTF head origin. Subsequently, this
//  routine applies beamformer imaging to the source & displays its
//  extent in three line transections through the simulated sources.
//
//  This version uses the local sphere parameters to move the simulated
//  brain noise & dipole source, accordingly.
//
//	flags:
//		-r <dataset name>			-- MEG dataset name
//      -s <dipole spec file>       -- file with simulated dipole parameters
//      -n                          -- simulated brain & sensor noise
//		-v							-- verbose mode
//
//		Author: Stephen E. Robinson
//				MEG Core Facility
//				NIMH
//


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <fcntl.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <filters.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <samutil.h>
#include <voxel.h>
#include <version.h>

#define MINOR_REV   6
#define NOISE       5.0e-15     // sensor noise density 5 fT/√Hz
#define IOTA        1.0e-12     // small distance increment to assure loop completion
#define OUTERRAD    0.08        // 8 cm outer radius
#define INNERRAD    0.06        // 6 cm inner radius
#define ORDER       3           // 3rd-order integration over coil area
#define RANGE       0.01        // +/- 1 cm
#define STEP        1.0e-04     // 0.1 mm steps

#define HISTO       TRUE

int	main(
    int     argc,
    char    **argv
)

{
    HeaderInfo  Header;         // MEG data header
    ChannelInfo *Channel;       // MEG channel info
    EpochInfo   *Epoch;         // MEG epoch info
    FFTSPEC     RejectFFT;      // band-reject filter
    FIELD       Q5;             // 5-parameter dipole for brain noise simulation
    FIELD       Q6;             // test dipole
    VOXELINFO   Dipole;         // test dipole
    gsl_matrix  *A;             // copy of C (for LU inversion)
    gsl_matrix  *C;             // covariance matrix - MxM
    gsl_matrix  *Cinv;          // inverse covariance matrix - MxM
    gsl_vector  *Bp;            // forward solution - Mx1
    gsl_vector  *W;             // beamformer weights - Mx1
    gsl_matrix  *V;             // for SVD
    gsl_vector  *S;             // singular values
    gsl_vector  *work;          // for SVD
    gsl_permutation *p;         // LU decomposition permutation
    double	rand_rad3;      // Random radius between INNERRAD^3 and OUTERRAD^3
    double      **x;            // x[M][TT] -- multisensor data time-series
    double      *In;            // In[T] -- input epoch time-series
    double      *Out;           // Out[T] -- output epoch time-series
    double      *dip;           // dip[TT] -- Gaussian random time-series
    double      *Brain;         // Brain[TT] -- simulated brain noise
    double      *BrainScaled;   // BrainScaled[TT] -- scaled simulated brain noise
    double      Vc[3];          // Cartesian vector
    double      Vd[3];          // dipole coordinates
    double      Vx[3];          // tangential vector x'
    double      Vy[3];          // tangential vector y'
    double      Vo[3];          // local sphere origin
    double 	Vnew[3];
    double 	Vold[3];
    double      Fwd;            // forward solution
    double      Moment;         // dipole moment (A-m)
    double      MomVec;         // tangential moment vector (rad)
    double      BandWidth;      // dataset bandwidth
    double      HPfreq;         // band-reject highpass frequency
    double      LPfreq;         // band-reject lowpass frequency
    double      CTFNoise;       // rms noise
    double      BrainMax;
    double      BrainMin;
    double      MaxBrainNoise;  // rms maximum brain noise
    double      MinBrainNoise;  // rms minimum brain noise
    //double      BrainScale;     // peakiness of Levy distribution
    double      BrainPow;     // change random dipole distribution
    double      BrainNoise;     // just more noise, as far as I'm concerned...
    double      Noise2;         // noise variance
    double      rms;            // root mean square
    double      Src;            // projected source power
    double      Nse;            // projected noise power
    double      scale;          // T/bit
    double      span;           // eigenvalue span
    double      r;				// distance
    double      r2;				// r^2
    double      min;            // minimum distance between sensors & dipole
    double      min_test;       // minimum distance between sensors & dipole for simulated data
    double      percent;        // percent done
    double      tmp;            // what else?
    double	CorrDist=0.;	// correlation distance
    int         NumNoise;       // number of brain noise dipole sources
    int         Modulus;        // repeat for brain noise time-series
    int         Seed;           // allow user to specify a seed for drand48()
    int         E;				// total number of epochs
    int         M;              // number of MEG or reference channels used
    int         T;              // number of samples
    int         TT;             // total samples
    int         c;              // vertex index
    int         e;              // epoch index
    int         f;              // 1st primary sensor index
    int         i;              // i-index
    int         j;              // j-index
    int         m;              // channel index
    int         mm;             // alternate channel index
    int         zmax_index;     // sensor with maximum z 
    double      zmax;           // sensor with maximum z 
    int         t;              // sample index
    int         tt;             // alternate sample index
    int         v;              // voxel index
    int         signum;         // sign of decomposition
    int         pct;            // integer percentage done
    int         LongIndex;      // position of arguments
    int         eflg = FALSE;   // command-line error flag
    int         bflg = FALSE;   // band reject flag
    int         nflg = FALSE;   // simulate all noise flag
    //int         oflg = FALSE;   // local sphere origin flag
    int         rflg = FALSE;   // dataset name flag
    int         sflg = FALSE;   // simulation parameter flag
    int         dflg = FALSE;   // seed parameter flag
    int         vflg = FALSE;   // verbose mode flag
    extern char *optarg;
    extern int  opterr;
    char        fpath[256];     // general path name
    char        DSName[256];    // MEG dataset name
    char        SimName[256];   // parameter file name
    char        DipName[256];   // dipole parameter file name
    char        Prefix[256];    // prefix for output file name
    char        DataPrefix[4];    // prefix for output file name
    char        Line[256];
    char        Label[64];
    unsigned char **Bad;        // Bad[E][T] -- flags for bad samples
    char        ShortOpts[] = "r:s:nv";
    static struct option    LongOpts[] = {  // command line options
        {"dataset_name", required_argument, NULL, 'r'},
        {"simulation_name", required_argument, NULL, 's'},
        {"sensor_noise", no_argument, NULL, 'n'},
        {"verbose", no_argument, NULL, 'v'},
        {NULL, no_argument, NULL, 0}
    };
    FILE        *fp;            // generic file pointer

    // parse command line parameters
    opterr = 1;             // enable error reporting
    while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
        switch(c) {
            case 'r':       // run name
                sprintf(DSName, "%s", optarg);
                rflg = TRUE;
                break;
            case 's':       // simulation parameter file name
                sprintf(SimName, "%s", optarg);
                sflg = TRUE;
                break;
		    case 'n':       // simulate brain & sensor noise
		        nflg = TRUE;
		        break;
			case 'v':       // verbose mode
				vflg = TRUE;
				break;
			default:
				eflg = TRUE;
				break;
		}
	}
	if(eflg || !rflg || !sflg) {
		fprintf(stderr, "SAMsim\t-r <run name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "\t-s <simulation parameter file name>\n");
		fprintf(stderr, "\t-n -- simulated brain & sensor noise\n");
		fprintf(stderr, "\t-v -- verbose mode\n");
		exit(-1);
	}

	// get dataset information structures
	if(vflg == TRUE) {
		printf("S A M   D A T A   S I M U L A T O R\t%s\n", PRG_REV);
		printf("opening data file");
		fflush(stdout);
	}

	// get data with sensor structures
	GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);

    // extract constants
    M = Header.NumPri;
    E = Header.NumEpochs;
    T = Header.MaxSamples;
    TT = E * T;
    f = Header.PsIndex[0];              // 1st MEG primary sensor channel index


    // read simulation parameters
    if (vflg) { 
        printf(" - done\n");
        printf("reading brain noise simulation parameters");
        fflush(stdout);
    }
    sprintf(fpath, "%s", SimName);
    if ((fp = fopen(fpath, "r")) == NULL)
        Cleanup("can't open simulation parameter file");
    while (fgets(Line, 128, fp) != NULL) {                  // scan parameter values
        if (strlen(Line) == 1)                              // skip empty lines
            continue;
        if (Line[0] == '#')                                 // skip comment lines	
            continue;
        if (sscanf(Line, "%s", Label) != 1)                 // read the label field
            cleanup("can't read label field");
        if (!strncmp("DipoleFile", Label, 10)) {            // read 'DipoleFile'
            if (sscanf(Line, "%s%s", Label, DipName) != 2)
                Cleanup("can't read 'DipName'");
        } else if (!strncmp("NumNoise", Label, 8)) {        // read 'NumNoise'
            if (sscanf(Line, "%s%d", Label, &NumNoise) != 2)
                Cleanup("can't read 'NumNoise'");
        } else if (!strncmp("BrainMax", Label, 8)) {        // read 'BrainMax'
            if (sscanf(Line, "%s%le", Label, &BrainMax) != 2)
                Cleanup("can't read 'BrainMax'");
        } else if (!strncmp("BrainMin", Label, 8)) {        // read 'BrainMin'
            if (sscanf(Line, "%s%le", Label, &BrainMin) != 2)
                Cleanup("can't read 'BrainMin'");
        //} else if (!strncmp("BrainScale", Label, 10)) {     // read 'BrainScale'
        //    if (sscanf(Line, "%s%lf", Label, &BrainScale) != 2)
        //        Cleanup("can't read 'BrainScale'");
        } else if (!strncmp("BrainPow", Label, 10)) {     // read 'BrainPow'
            if (sscanf(Line, "%s%lf", Label, &BrainPow) != 2)
                Cleanup("can't read 'BrainPow'");
        } else if (!strncmp("Modulus", Label, 7)) {         // read 'Modulus'
            if (sscanf(Line, "%s%d", Label, &Modulus) != 2)
                Cleanup("can't read 'Modulus'");
        } else if (!strncmp("CorrDist", Label, 8)) { 	    // read 'CorrDist'
            if (sscanf(Line, "%s%lf", Label, &CorrDist) != 2)
                Cleanup("can't read 'CorrDist");
            CorrDist *= 0.001;      // convert mm to m
        } else if (!strncmp("BandReject", Label, 10)) {     // read 'BandReject'
            if (sscanf(Line, "%s%lf%lf", Label, &HPfreq, &LPfreq) != 3)
                Cleanup("can't read 'BandReject' frequencies");
            bflg = TRUE;
        } else if (!strncmp("Seed", Label, 4)) {           // read 'Seed'
            if (sscanf(Line, "%s%d", Label, &Seed) !=2)
              Cleanup("can't read 'Seed'");
            dflg = TRUE;
        } else {
            fprintf(stderr, "\nBad label: '%s'\n", Label);
            fflush(stdout);
        }
    }
    fclose(fp);

    // set prefix
    //if (nflg)
    //    sprintf(Prefix, "%s,SIM", SimName);
    //else
    //    sprintf(Prefix, "%s,MEG", SimName);
    //
    // generate prefix
    //
    memset(DataPrefix, 0, sizeof(DataPrefix));
    memcpy(DataPrefix, Header.SetName, 3);
    //printf("\nDataPrefix: %s DipName: %s SimName: %s Header.setname: %s\n",DataPrefix, DipName, SimName, Header.SetName);
    if (!nflg)
       sprintf(Prefix, "%s_%s_%s_MEG", DataPrefix, DipName, SimName);
    else
       sprintf(Prefix, "%s_%s_%s_SIM", DataPrefix, DipName, SimName);
    //

    // set sensor integration points & local sphere from the dataset

    if (vflg) {
        printf(" - done\n");
        printf("setting sensor coil integration points & local sphere");
        fflush(stdout);
    }
    for (m=0; m<Header.NumChannels; m++)
        if (Channel[m].ChannelType == TYPE_REF_MAG || Channel[m].ChannelType == TYPE_REF_GRAD || Channel[m].ChannelType == TYPE_MEG) {
            Channel[m].BalOrder = 0;
            Channel[m].NumBal = 0;
            FwdSensIntPnt(&Channel[m].Geom.MEGSensor, ORDER);
        }

    // Read dipole parameters

    if (vflg) { 
        printf(" - done\n");
        printf("reading dipole parameters");
        fflush(stdout);
    }
    sprintf(fpath, "%s", DipName);
    if ((fp = fopen(fpath, "r")) == NULL)
        Cleanup("can't open dipole specification file");
    if (fscanf(fp, "%le%le%le%le", &Vd[X_], &Vd[Y_], &Vd[Z_], &Moment) != 4)
        Cleanup("can't read dipole parameters");
    for (v=X_; v<=Z_; v++)
        Vd[v] *= 0.01;      // convert spec'd cm to m
    Moment *= 1.0e-09;      // convert spec'd moment to nA-m (rms)
    fclose(fp);
    
    // compute distance from nearest primary sensor to dipole
    // also find the sensor with max z

    for (m=0, min=HUGE, zmax=0; m<M; m++) {
        mm = Header.PsIndex[m];
        for (v=X_, r2=0.; v<=Z_; v++) {
            tmp = Channel[mm].Geom.MEGSensor.origin.p[v] - Vd[v];
            r2 += tmp * tmp;
        }
        r = sqrt(r2);
        if (r < min) 
            min = r;
        if (Channel[mm].Geom.MEGSensor.origin.p[Z_] > zmax) { 
            zmax = Channel[mm].Geom.MEGSensor.origin.p[Z_];
            zmax_index = mm;
        }
    }
    printf("\nMEG dataset has minimum distance from nearest sensor to simulated dipole: %7.3f mm", 1000.*min);

    // if we are simulating brain noise...
 
    if (nflg) { 

       // set the test dipole directly below the sensor with maximum z value, the same distance away
       // as the minimum distance from the sensors to the test dipole from the real brain data ("min")

       printf("\nFor simulated data, setting test dipole %3.3f mm from sensor %d at maximum z=%3.3f mm\n", 1000*min, zmax_index, 1000*zmax);  

       Vd[X_] = Channel[zmax_index].Geom.MEGSensor.origin.p[X_];
       Vd[Y_] = Channel[zmax_index].Geom.MEGSensor.origin.p[Y_];
       Vd[Z_] = zmax - min; 

       printf("Dipole Coords Vd[x]=%3.3f Vd[y]=%3.3f Vd[z]=%3.3f\n", Vd[X_],Vd[Y_],Vd[Z_]);
 
       // we also have to shift the local sphere origin so that the sphere of noise sources will also be shifted accordingly
      
       Vo[X_] = Channel[zmax_index].Geom.MEGSensor.origin.p[X_];
       Vo[Y_] = Channel[zmax_index].Geom.MEGSensor.origin.p[Y_];
       Vo[Z_] = zmax - (((OUTERRAD+INNERRAD)/2) + min); 

       printf("Shifting noise sphere origin accordingly: Vo[x]=%3.3f Vo[y]=%3.3f Vo[z]=%3.3f\n", Vo[X_],Vo[Y_],Vo[Z_]); 

       // figure out the new minimum distance - because of the particular geometry, a neighbor sensor may now be slightly closer

       for (m=0, min_test=HUGE; m<M; m++) {
           mm = Header.PsIndex[m];
           for (v=X_, r2=0.; v<=Z_; v++) {
              tmp = Channel[mm].Geom.MEGSensor.origin.p[v] - Vd[v];
              r2 += tmp * tmp;
           }
           r = sqrt(r2);
           if (r < min_test) {
               min_test = r;
           }
        } 
        printf("New minimum distance to any sensor is %3.3f mm\n",1000*min_test);
    }

    if (!nflg)
       if (vflg) {
           printf(" - done\n");
           printf("distance from nearest sensor to simulated dipole: %7.3f mm", 1000.*min);
           fflush(stdout);
       }
    
    // make band-reject filter
    if (vflg && bflg) {
        printf(" - done\n");
        printf("making %g to %g Hz band-reject filter", HPfreq, LPfreq);
        fflush(stdout);
    }
    if (bflg)
        Butterworth(&RejectFFT, T, HPfreq, LPfreq, Header.SampleRate, BANDREJECT, MAXORDER, &BandWidth);
    else
        BandWidth = Channel[f].AnalogLpFreq - Channel[f].AnalogHpFreq;

    // set bandwidth-dependent parameters
    CTFNoise = NOISE * sqrt(BandWidth);
    Noise2 = CTFNoise * CTFNoise;       // noise power
    MaxBrainNoise = BrainMax * sqrt(BandWidth);
    MinBrainNoise = BrainMin * sqrt(BandWidth);

    // allocate matrices & vectors
    if (vflg) {
        printf(" - done\n");
        printf("allocating matrices & vectors");
        fflush(stdout);
    }
    x = new_matrixE(M, TT, "x[M][TT]");             // allocate memory for meaured MEG signals
    In = new_arrayE(double, T, "In[T]");            // allocate memory for single epoch
    Out = new_arrayE(double, T, "Out[T]");          // allocate memory for single epoch
    dip = new_arrayE(double, TT, "dip[TT]");        // allocate memory for simulated moment time-series
    Brain = new_arrayE(double, TT, "Brain[TT]");    // allocate memory for simulated brain noise
    BrainScaled = new_arrayE(double, TT, "BrainScaled[TT]");    // allocate memory for simulated brain noise
    C = gsl_matrix_alloc(M, M);
    Cinv = gsl_matrix_alloc(M, M);
    Bp = gsl_vector_alloc(M);
    W = gsl_vector_alloc(M);
    A = gsl_matrix_alloc(M, M);
    p = gsl_permutation_alloc(M);

    // either read MEG data or simulate brain noise
    if (!nflg) {

        // read MEG data
        if (vflg) {
            printf(" - done\n");
            printf("reading MEG data");
            fflush(stdout);
        }
        for (m=0; m<M; m++) {                               // read & process data in channel order
            mm = Header.PsIndex[m];
            scale = Channel[mm].Scale;
            for (e=tt=0; e<E; e++) {                        // assemble multi-epoch contiquous data for this channel
                GetDsData(&Header, e, mm, scale, In);       // get one epoch of data for this channel -- length T samples
                if (!bflg) {
                    for (t=0; t<T; t++, tt++)
                        x[m][tt] = In[t];                      // 'x' is contiguous
                } else {
                    FFTfilter(In, Out, &RejectFFT, T);
                    for (t=0; t<T; t++, tt++)
                        x[m][tt] = Out[t];                      // 'x' is contiguous
                }
            }
            //demean(x[m], T); // pretty sure this is wrong
            demean(x[m],TT);

            // channel heartbeat
		    if (vflg) {
			    printf(".");
			        fflush(stdout);
		    }
	    }

	} else {        // simulate!

        // simulate sensor & brain noise
        if (vflg) {
            printf(" - done\n");
            printf("simulating sensor noise");
            fflush(stdout);
        }
        for (m=0; m<M; m++)
            for (t=0; t<T; t++)
                x[m][t] = RandomGauss(0., CTFNoise);
#ifdef HISTO
        // file BrainNoise for histogram plot
        if (nflg) {
            sprintf(fpath, "%s_BrainNoise.txt", Prefix);
            if ((fp = fopen(fpath, "w")) == NULL)
                Cleanup("can't open 'BrainNoise' for write");
        }
#endif

        // simulate brain noise in shell
        if (vflg) {
            printf(" - done\n");
            printf("simulating brain noise time-series (please be patient)\n");
            //if (oflg)
            printf("Brain noise shell extends up to max z of %3.3f mm with dipole at %3.3f mm\n", 1000.*(OUTERRAD+Vo[Z_]), 1000.*Vd[Z_]);
            fflush(stdout);
        }

        if (dflg) { 
           srand48(Seed);
           printf("Using random seed %d for noise dipole position and orientation\n", Seed );
        } else { 
           srand48((long)time((time_t *)0));   // seed random number generator
           printf("Using clock time as seed, first rand = %f\n", drand48());
        }
  
        if (CorrDist > 0.)	// disable modulus if correlation distance is set
	   Modulus = 1.;
        for (v=X_; v<=Z_; v++)
            Vold[v] = 0.;

        printf("new rand: %f\n",drand48());
        for (i=c=0, pct=(-1); i<NumNoise; i++) {

            // generate a random unit dipole within the shell, in spherical coordinates -- shells are relative to CTF head frame origin

            Q5.p[THETA_] = 2. * M_PI * drand48();       // 0 <= theta <= 2*pi
            Q5.p[PHI_] = M_PI;				// initially set phi in the Southern hemisphere 
            while(Q5.p[PHI_] > M_PI / 2. ) 		// keep picking phi values until you get one in the Northern hemisphere 
              Q5.p[PHI_] = acos(2. * drand48() -1.);    // northern hemisphere: 0 <= phi <= pi/2
            rand_rad3 = (pow(OUTERRAD,3) - pow(INNERRAD,3)) * drand48() + pow(INNERRAD,3);
            Q5.p[RHO_] = cbrt(rand_rad3);	        // area of a sphere scales like r^3 - need equal dipoles per unit area!

            Q5.v[PSI_] = 2. * M_PI * drand48();         // 0 <= psi <= 2*pi
            Q5.v[Q_] = 1.;                              // |Q| = 1.0
            P5toP6(&Q5, &Q6);                           // convert this to 6-parameter cartesian dipole

            // don't forget to shift the noise source sphere
            for (v=X_; v<=Z_; v++)
                Q6.p[v] += Vo[v];

   	    // compute distance between this brain noise source & the previous one
            for (v=X_, r2=0.; v<=Z_; v++) {
                Vnew[v] = Q6.p[v];
                tmp = Vnew[v] - Vold[v];
                r2 += tmp * tmp;
                Vold[v] = Vnew[v];
            }
            r = sqrt(r2);

            // simulate source field using Sarvas equation for unit (1.0 A-m dipole)
            for(m=0; m<M; m++) {
                mm = Header.PsIndex[m];                     // make sure we are selecting primary sensors
                gsl_vector_set(Bp, m, FwdSensor(&Q6, &Channel[mm].Geom.MEGSensor, ORDER));
            }

            // for each time sample compute brain noise time-series -- this is the tricky part!

            BrainNoise = (pow((double)i / (double)NumNoise, BrainPow) * (MaxBrainNoise-MinBrainNoise)) + MinBrainNoise;

            //do
            //    BrainNoise = MaxBrainNoise * Levy((20. * drand48() + 1.1e-12), 1.0e-12, BrainScale);
            ////while (BrainNoise < MinBrainNoise);
#ifdef HISTO
            if (nflg)
                fprintf(fp, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i, BrainNoise,Q6.p[0],Q6.p[1],Q6.p[2],Q6.v[0],Q6.v[1],Q6.v[2]);
#endif
            if (r > CorrDist && i % Modulus == 0) {  // only generate a new Brain[t] timeseries if Modulus = 0
                for (t=0; t<TT; t++)
                    Brain[t] = RandomGauss(0., 1);
	    } else {
                c++;
            }

            // scale the existing Brain[t] by the noise value for this dipole - this may be a new random time series or 
            // a repeated one, depending on the modulus 
           
            for (t=0; t<TT; t++)
                 BrainScaled[t] = BrainNoise * Brain[t]; 

            // sum this dipole MEG signal
            for (m=0; m<M; m++)
                for(t=0; t<TT; t++)
                    x[m][t] += BrainScaled[t] * gsl_vector_get(Bp, m);

            // show progress (so I can be bored to tears)
            if (vflg) {
                percent = 100. * (double)i / (double)NumNoise;
                if(rint(percent) > (double)pct) {
                    pct = (int)rint(percent);
                    printf("\r\t%d percent done", pct);
                    fflush(stdout);
                }
            }
        }
    }
#ifdef HISTO
    if (nflg)
        fclose(fp);
    free(Brain);
    free(BrainScaled);
    printf("\nnumber of correlated sources: %d\n",c);
#endif
	// summarize noise statistics
	sprintf(fpath, "%s_Noise.txt", Prefix);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open noise statistics for write");

    // compute rms sensor noise plus brain noise
    for (m=0, tmp=0.; m<M; m++)
        for (t=0; t<TT; t++)
            tmp += x[m][t] * x[m][t];
    rms = sqrt(tmp / (double)(M * TT));
    if (nflg) {
        fprintf(fp, "simulated brain noise (%d sources): %6.1f fT rms\n", NumNoise, 1.0e+15*rms);
        printf("simulated brain noise (%d sources): %6.1f fT rms\n", NumNoise, 1.0e+15*rms);
    } else {
        fprintf(fp, "measured brain noise: %6.1f fT rms\n", 1.0e+15*rms);
        printf("measured brain noise: %6.1f fT rms\n", 1.0e+15*rms);
    }
    fprintf(fp, "brain noise density: %5.1f fT/√Hz\n", 1.0e+15*rms/sqrt(BandWidth));
    printf("brain noise density: %5.1f fT/√Hz\n", 1.0e+15*rms/sqrt(BandWidth));
    fclose(fp);

    // compute unit radial vector from sphere origin to dipole
    if (vflg) {
        printf(" - done\n");
        printf("setting tangential test dipole moment vector");
        fflush(stdout);
    }
    for (v=X_, r2=0.; v<=Z_; v++)   // compute unit radial vector
        r2 += Vd[v] * Vd[v];
    r = sqrt(r2);
    for (v=X_; v<=Z_; v++)
        Vc[v] = Vd[v] / r;
    OrthoPlane(Vc, Vx, Vy);         // compute tangential unit vectors Vx & Vy
    MomVec = 2. * M_PI * drand48(); // random vector from 0 to 2 pi
//    MomVec = M_PI / 4.;             // 45 degree fixed orientation
    for (v=X_; v<=Z_; v++) {
        Dipole.p[v] = Vd[v];
        Dipole.v[v] = Vx[v] * cos(MomVec) + Vy[v] * sin(MomVec);
    }
    Dipole.Solve = TRUE;
    Dipole.ROI = TRUE;
#ifdef EBUG
    if (vflg) {
        printf(" - done\n");
        printf("%s: %f %f %f %f %f %f\n", SimName,
          1000.*Dipole.p[X_], 1000.*Dipole.p[Y_], 1000.*Dipole.p[Z_], Dipole.v[X_], Dipole.v[Y_], Dipole.v[Z_]);
        fflush(stdout);
    }
#endif
    // convert dipole parameters to FIELD
    for (v=X_; v<=Z_; v++) {
        Q6.p[v] = Dipole.p[v];
        Q6.v[v] = Dipole.v[v];
    }

    // compute Gaussian random moment time series for dipole
    if (vflg) {
        printf(" - done\n");
	    printf("computing Gaussian random moment time-series");
	    fflush(stdout);
	}
    for (t=0; t<TT; t++)
        dip[t] = RandomGauss(0., Moment);   // Moment is A-m rms

    // compute forward solution for dipole in a sphere (Sarvas) & sum random waveform to data
    if (vflg) {
	    printf(" - done\n");
	    printf("summing simulated dipole to MEG data");
	    fflush(stdout);
	}
	for (m=0; m<M; m++) {
	    mm = Header.PsIndex[m];
	    Fwd = FwdSensor(&Q6, &Channel[mm].Geom.MEGSensor, ORDER);
        for (t=0; t<TT; t++)
            x[m][t] += Fwd * dip[t];
	}

	// compute covariance
	if (vflg) {
	    printf(" - done\n");
	    printf("computing covariance matrix & its inverse");
	    fflush(stdout);
	}
	for (i=0; i<M; i++)
		for (j=0; j<M; j++) {
			for (t=0, tmp=0.; t<TT; t++)
			    tmp += x[i][t] * x[j][t];
            gsl_matrix_set(C, i, j, tmp / (double)TT);
		}
    gsl_matrix_memcpy(A, C);    // make a copy of C
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_invert(A, p, Cinv);

    // compute eigenvalue spectrum using SVD
    if (vflg) {
	    printf(" - done\n");
	    printf("computing eigenvalue spectrum of covariance matrix");
	    fflush(stdout);
	}
    gsl_matrix_memcpy(A, C);    // make another copy of C
    V = gsl_matrix_alloc(M, M);
    S = gsl_vector_alloc(M);
    work = gsl_vector_alloc(M);
    gsl_linalg_SV_decomp(A, V, S, work);
    sprintf(fpath, "%s_Eigenvalues.txt", Prefix);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'eigenvalues' for write");
    for (i=0; i<M; i++)
        fprintf(fp, "%e\n", gsl_vector_get(S, i));
    fclose(fp);
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
        
    // scan SAM scan along x-axis through true source location
    if (vflg) {
	    printf(" - done\n");
	    printf("computing x-axis scan");
	    fflush(stdout);
	}
    sprintf(fpath, "%s_Xscan.txt", Prefix);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'Xscan' for write");
    for (r=(-RANGE); r<=(RANGE+IOTA); r+=STEP) {
        Dipole.p[X_] = r + Vd[X_];
        Dipole.p[Y_] = Vd[Y_];
        Dipole.p[Z_] = Vd[Z_];

        // compute source z^2 & save to file
        DIPSolve2D(&Dipole, Bp, Cinv, &Header, Channel, ORDER, &span);
        SAMsolve(C, Cinv, Bp, W, Noise2, &Src, &Nse);
        fprintf(fp, "%f\t%f\t%f\n", 1000.*r, 1.0e+09*sqrt(Src), Src/Nse);
        fflush(fp);

        // heartbeat
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }
    fclose(fp);

    // scan SAM output along y-axis through true source location
    if (vflg) {
	    printf(" - done\n");
	    printf("computing y-axis scan");
	    fflush(stdout);
	}
    sprintf(fpath, "%s_Yscan.txt", Prefix);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'Yscan' for write");
    for (r=(-RANGE); r<=(RANGE+IOTA); r+=STEP) {
        Dipole.p[X_] = Vd[X_];
        Dipole.p[Y_] = r + Vd[Y_];
        Dipole.p[Z_] = Vd[Z_];

        // compute source z^2 & save to file
        DIPSolve2D(&Dipole, Bp, Cinv, &Header, Channel, ORDER, &span);
        SAMsolve(C, Cinv, Bp, W, Noise2, &Src, &Nse);
        fprintf(fp, "%f\t%f\t%f\n", 1000.*r, 1.0e+09*sqrt(Src), Src/Nse);
        fflush(fp);

        // heartbeat
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }
    fclose(fp);

    // scan SAM output along z-axis through true source location
    if (vflg) {
	    printf(" - done\n");
	    printf("computing z-axis scan");
	    fflush(stdout);
	}
    sprintf(fpath, "%s_Zscan.txt", Prefix);
    if ((fp = fopen(fpath, "w")) == NULL)
        Cleanup("can't open 'Zscan' for write");
    for (r=(-RANGE); r<=(RANGE+IOTA); r+=STEP) {
        Dipole.p[X_] = Vd[X_];
        Dipole.p[Y_] = Vd[Y_];
        Dipole.p[Z_] = r + Vd[Z_];

        // compute source z^2 & save to file
        DIPSolve2D(&Dipole, Bp, Cinv, &Header, Channel, ORDER, &span);
        SAMsolve(C, Cinv, Bp, W, Noise2, &Src, &Nse);
        fprintf(fp, "%f\t%f\t%f\n", 1000.*r, 1.0e+09*sqrt(Src), Src/Nse);
        fflush(fp);

        // heartbeat
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }
    fclose(fp);

	// we're done now!
	if(vflg == TRUE) {
	    printf(" - done\n");
		printf("'SAMsimulate' done\n");
		fflush(stdout);
	}
	exit(0);
}
