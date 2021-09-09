// CTFsim -- simulate solution of cortical source activity &
//  independent source resolution, given an array of sensors
//  (CTF-275) & an array of test dipole locations. Thus simulation
//  employs the forward solution for a dipole in a homogeneous
//  spherical conductor (Sarvas).
//
//  Note: The average skull thickness for adult males is 6.5 mm
//      & for adult females is 7.1 mm.
//
//	Author:	Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>
#include <siglib.h>
#include <samlib.h>
#include <DataFiles.h>
#include <samutil.h>
#include <filters.h>

#define HEAD_RADIUS     0.09        // 9 cm head radius
#define DURATION        300.0       // 5 minute duration
#define SAMPLERATE      600.0       // 600 Hz sample rate
#define BANDWIDTH       300.0       // 300 Hz bandwidth
#define NOISE           5.0e-15     // 5 fT/√Hz noise density
#define DIPMOMENT       3.5e-08     // default 35 nA-m rms dipole moment
#define OUTERRAD        0.08        // 8 cm outer radius
#define INNERRAD        0.06        // 6 cm inner radius
#define NUMNOISE        2000       // 10,000 brain noise sources
#define BRAINMAX        1.0e-09     // maximum brain noise moment
#define BRAINMIN        1.0e-12     // minimum brain noise
#define BRAINPOW        6.0         // exponent for approximating brain noise
#define ORDER           3           // 3rd-order integration over coil area
#define MODULUS         8           // Modulus value for correlated sources
#define IOTA		1.0e-09	    // small distance for loop completion

int
main(
    int             argc,
    char            **argv
)

{
    // variables
    HeaderInfo  Header;         // header information data
    ChannelInfo *Channel;       // Channel[M] -- channel information data
    ChannelInfo *NewChannel;    // NewChannel[M] -- channel information data
	EpochInfo   *Epoch;         // Epoch[E] -- epoch information data
    FIELD       Q6;             // 6-parameter test dipole
    FIELD       Q5;             // 5-parameter test dipole
    VOXELINFO   *Dipole;        // Dipole[D] -- array of test dipoles
    XFORM       Where;          // translation rotation matrix
    gsl_matrix  *C;             // MxM covariance matrix
    gsl_matrix  *Cinv;          // MxM inverse covariance matrix
    gsl_vector  *Bp;            // Mx1 forward solution vector
    gsl_vector  *W;             // Mx1 beamformer coeffcient vector
    double 	CorticalDepth=0.;  // depth of test dipole in noise shell
    double      rand_rad3;        // / Random radius between INNERRAD^3 and OUTERRAD^3
    double      **data;         // data[M][T] -- MEG data time-series
    double      **dip;          // dip[D][T] -- dipole source time series
    double      **src;          // src[D][T] -- SAM time-series
    double      **forward;      // forward[D][M] -- dipole forward solution matrix
    double      *Brain;         // Brain[T] -- brain noise time-series
    double      *BrainScaled;         // BrainScaled[T] -- scaled brain noise time-series
    double      *Gain;          // Gain[D] -- recovered waverform gain -- relative to simulated waveform
    double      **distance;     // distance[D][D] -- distances between dipoles
    double      **correlation;  // correlation[D][D] -- correlation between dipole & source time-series
    double      Vc[3];          // cartesian vector
    double      Rotation[3];    // rotation (Euler angles)
    double      Vs[3];          // spherical vector
    double      Vb[3];          // baseline vector
    double      Vx[3];          // x' tangential vector
    double      Vy[3];          // y' tangential vector
    double      Dx, Dy;         // test dipole coordinates
    double      MaxSensor;      // maximum z-coordinate of sensor array
    double      MaxDipole;      // maximum z-coordinate of test dipole array
    double      interval;       // dipole spacing
    double      CTFNoise;       // rms CTF sensor noise
    double      MaxBrainNoise;  // rms maximum brain noise
    double      MinBrainNoise;  // rms minimum brain noise
    double      Noise2;         // sensor noise power
    double      BrainNoise;     // just more noise, as far as I'm concerned...
    double      Src;            // source power
    double      Nse;            // noise power
    double      span;           // span of eigensystem solution
    double      sx;             // sum of x
    double      sx2;            // sum of x^2
    double      sy;             // sum of y
    double      sxy;            // sum of x*y
    double      isr;            // independent source resolution
    double      r, r2, rvar;
    double      rms;            // rms
    double      angle;          // angle (degrees)
    double      MomVec;         // tangential moment vector (radians)
    double      Moment = DIPMOMENT; // default dipole moment
    double      percent;        // percent done
    double      edge;           // source array edge
    double      range;          // image range (edge + 10mm)
    double      tmp;            // accumulator of sums
    int         D;              // number of dipoles
    long        seed;           // random number seed
    int         M;              // number of primary sensors
    int         R;              // number of reference sensors
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
    int         pct;            // integer percentage done
    int         LongIndex;      // position of arguments
    int         dflg = FALSE;   // dataset flag
    int         eflg = FALSE;   // edge flag
    int         fflg = FALSE;   // fixed source grid flag
    int         nflg = FALSE;   // source number flag
    int         rflg = FALSE;   // random moment vector flag
    int         vflg = FALSE;   // verbose flag
    int         sflg = FALSE;   // random seed flag
    int         zflg = FALSE;   // output filename stem
    int         errflg = FALSE; // error flag
    char        DSName[64];     // dataset name
    char        prefix[32];     // prefix for output files
	char        output[32];     // output file name
	FILE        *fp;            // generic file pointer
	FILE        *flog;          // log file pointer
	extern char *optarg;
	extern int  opterr;
    char        ShortOpts[] = "d:e:n:m:s:c:z:frv";
	static struct option	LongOpts[] = {                  // command line options
		{"dataset_name", required_argument, NULL, 'd'},     // name of file with sensor coordinates
		{"array_edge", required_argument, NULL, 'e'},       // source array edge (mm)
		{"number_dipoles", required_argument, NULL, 'n'},   // number of dipole sources
		{"dipole_moment", required_argument, NULL, 'm'},    // dipole moment (nA-m)
		{"fixed_dipoles", no_argument, NULL, 'f'},          // fixed dipole grid
		{"random", no_argument, NULL, 'r'},                 // random source moment vector
		{"verbose", no_argument, NULL, 'v'},
                {"seed", required_argument, NULL, 's'},
                {"cortical_depth",required_argument, NULL, 'c'},    // depth of test dipole in noise shell
                {"prefix", required_argument, NULL, 'z'},
		{NULL, no_argument, NULL, 0}
	};

	// parse command line parameters
	opterr = 1;				// enable error reporting
	while((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
		switch (c) {
			case 'd':		// sensor file name
				sprintf(DSName, "%s", optarg);
				dflg = TRUE;
				break;
			case 'e':       // area edge
			    edge = 0.001 * atof(optarg);
			    eflg = TRUE;
			    break;
			case 'n':		// number of sources
			    D = atoi(optarg);
				nflg = TRUE;
				break;
			case 'm':       // dipole moment (optional)
			    Moment = 1.0e-09 * atof(optarg);
			    break;
			case 'c':       // cortical depth  (optional)
			    CorticalDepth = .001 * atof(optarg);
			    break;
            		case 'f':       // fixed dipole grid -- otherwise random positions
                	    fflg = TRUE;
                	    break;
 			case 's': 	// use a specified seed instead of system time
 		 	    sflg = TRUE;
			    seed = atoi(optarg);
		 	    break;
 			case 'z': 	// use a specified seed instead of system time
 		 	    zflg = TRUE;
			    sprintf(prefix, "%s", optarg);
		 	    break;
			case 'r':       // generate random source moment vectors -- else parallel moments
			    rflg = TRUE;
			    break;
			case 'v':       // verbose output
			    vflg = TRUE;
			    break;
			default:
				errflg = TRUE;
				break;
		}
	}
	if(errflg == TRUE || dflg == FALSE || nflg == FALSE || eflg == FALSE) {
		fprintf(stderr, "CTFsim [arguments]\tv1.0\n");
		fprintf(stderr, "\t-d <dataset name>\n");
		fprintf(stderr, "\t-n <number of dipole sources\n");
		fprintf(stderr, "\t-m <dipole moment (nA-m)\n");
		fprintf(stderr, "\t-e <edge (mm)>\n");
		fprintf(stderr, "\t-f -- fixed dipole grid (not random)\n");
		fprintf(stderr, "\t-r -- generate random dipole moment vectors\n");
		fprintf(stderr, "\t-v -- verbose\n");
		exit(-1);
	}

	// announce program
	if (vflg) {
	    printf("CTFMsim version 1.0\n");
	    fflush(stdout);
	}

	// get data with sensor structures
	printf("reading sensor information from CTF 275-channel dataset");
	fflush(stdout);
	GetDsInfo(DSName, &Header, &Channel, &Epoch, NULL, TRUE);
	GetDsInfo(DSName, &Header, &NewChannel, &Epoch, NULL, TRUE);    // cloned ChannelInfo
    M = Header.NumPri;
    R = Header.NumRef;

	// initialize header parameters
	if (vflg) {
	    printf(" - done\n");
	    printf("initializing HeaderInfo with simulation constants");
	    fflush(stdout);
	}
	Header.NumEpochs = 1;
	Header.MaxSamples = (int)(SAMPLERATE * DURATION);
	Header.SampleRate = SAMPLERATE;
	Header.PowerFreq = 60.;
	Header.DataFormat = 0;

    // initialize sensor integration points & set local-sphere origins to zero
    if (vflg) {
	    printf(" - done\n");
	    printf("read %d primary & %d reference sensors\n", M, R);
	    printf("translating sensors & setting sensor coil integration points");
	    fflush(stdout);
    }
    Where.Translate[X_] = 0.0;    // shift sensors 50 mm posterior
    Where.Translate[Y_] = 0.;
    Where.Translate[Z_] = -0.065;    // shift sensors 120 mm superior
    Rotation[PSI_] = 0.;
    Rotation[THETA_] = 0.;
    Rotation[PHI_] = 0.;
    EtoR(Rotation, Where.Rotate);
    MoveMEG(NewChannel, Channel, Header.SqIndex, Header.NumSquid, &Where);
    for (m=0, MaxSensor=0.; m<Header.NumChannels; m++)
        if (Channel[m].ChannelType == TYPE_REF_MAG || Channel[m].ChannelType == TYPE_REF_GRAD || Channel[m].ChannelType == TYPE_MEG) {

            // correct sensor polarities & set order to raw
            Channel[m].BalOrder = 0;
            Channel[m].NumBal = 0;
            Channel[m].Geom.MEGSensor.Coil[0].SenseTurns = 1.;
            if (Channel[m].Geom.MEGSensor.NumCoils == 2)
                Channel[m].Geom.MEGSensor.Coil[1].SenseTurns = -1.;

            // zero sphere origins
            for (v=X_; v<=Z_; v++)
                Channel[m].Geom.MEGSensor.LocalSphere[v] = 0.;

            // set coiil integration points
            FwdSensIntPnt(&Channel[m].Geom.MEGSensor, ORDER);

            // get maximum z-coordinate of sensor array
            if (Channel[m].ChannelType == TYPE_MEG)
                if (Channel[m].Geom.MEGSensor.origin.p[Z_] > MaxSensor)
                    MaxSensor = Channel[m].Geom.MEGSensor.origin.p[Z_];
        }

    //initialize some constants
    srand48((long)time((time_t *)0));   // seed random number generator
    T = Header.MaxSamples;
    CTFNoise = NOISE * sqrt(BANDWIDTH);
    MaxBrainNoise = BRAINMAX * sqrt(BANDWIDTH);
    MinBrainNoise = BRAINMIN * sqrt(BANDWIDTH);
    Noise2 = CTFNoise * CTFNoise;       // noise power

	// allocate space for data (excluding the reference channels)
	if (vflg) {
	    printf(" - done\n");
	    printf("allocating data array");
	    fflush(stdout);
	}
    data = new_matrixE(M, T, "data[M][T]");

	// clear data time series for primary sensors, only
	for (m=0; m<M; m++)
	    for (t=0; t<T; t++)
	        data[m][t] = 0.;

	// simulate brain noise in shell
	// ***** NOTE: THE PARAMETERS USED TO SIMULATE BRAIN NOISE WERE DERIVED EXPERIMENTALLY  *****
	// ***** TO RESEMBLE MEASURED BRAIN NOISE. CHANGE PARAMETERS AT YOUR OWN RISK           *****
	if (vflg) {
	    printf(" - done\n");
	    printf("simulating brain noise time-series\n");
	    fflush(stdout);
	}
	Bp = gsl_vector_alloc(M);

	Brain = new_arrayE(double, T, "Brain[T]");
	BrainScaled = new_arrayE(double, T, "BrainScaled[T]");
	for (i=0, pct=(-1); i<=NUMNOISE; i++) {

		// generate a random unit dipole within the shell, in spherical coordinates
                Q5.p[THETA_] = 2. * M_PI * drand48();       // 0 <= theta <= 2*pi
                Q5.p[PHI_] = M_PI;                          // initially set phi in the Southern hemisphere
                while(Q5.p[PHI_] > M_PI / 2. )              // keep picking phi values until you get one in the Northern hemisphere
                  Q5.p[PHI_] = acos(2. * drand48() -1.);    // northern hemisphere: 0 <= phi <= pi/2
                rand_rad3 = (pow(OUTERRAD,3) - pow(INNERRAD,3)) * drand48() + pow(INNERRAD,3);
                Q5.p[RHO_] = cbrt(rand_rad3);               // area of a sphere scales like r^3 - need equal dipoles per unit area!

		Q5.v[PSI_] = 2. * M_PI * drand48();				// 0 <= psi <= 2*pi
		Q5.v[Q_] = 1.;									// |Q| = 1.0

		// convert this to 6-parameter cartesian dipole
		P5toP6(&Q5, &Q6);

		// for each time sample compute brain noise time-series -- this is the tricky part!
		BrainNoise = (pow((double)i / (double)NUMNOISE, BRAINPOW) * (MaxBrainNoise-MinBrainNoise)) + MinBrainNoise;

                if (i % MODULUS == 0)
		    for (t=0; t<T; t++)
			Brain[t] = RandomGauss(0., 1);

                for (t=0; t<T; t++)
                        BrainScaled[t] = BrainNoise * Brain[t];

		// simulate source field using Sarvas equation for unit (1.0 A-m dipole)
		for(m=0; m<M; m++) {
		    mm = Header.PsIndex[m];                     // make sure we are selecting primary sensors
		    gsl_vector_set(Bp, m, FwdSensor(&Q6, &Channel[mm].Geom.MEGSensor, ORDER));
		}

		// sum this dipole MEG signal
		for (m=0; m<M; m++)
			for(t=0; t<T; t++)
			    data[m][t] += BrainScaled[t] * gsl_vector_get(Bp, m);

		// show progress (so I can be bored to tears)
		if (vflg) {
		    percent = 100. * (double)i / (double)NUMNOISE;
		    if(rint(percent) > (double)pct) {
			    pct = (int)rint(percent);
			    printf("\r\t%d percent done", pct);
			    fflush(stdout);
		    }
		}
	}
	free(Brain);

	// open log file
	if (zflg)
          sprintf(output, "%s_CTFsim_Log.txt", prefix);
        else
	  sprintf(output, "CTFsim_Log,a%d.txt", M);
	if((flog = fopen(output, "w")) == NULL)
		Cleanup("can't open CTFsim log for write");

    // compute rms sensor noise plus brain noise
    for (m=0, tmp=0.; m<M; m++)
        for (t=0; t<T; t++)
            tmp += data[m][t] * data[m][t];
    rms = sqrt(tmp / (double)(M * T));
    fprintf(flog, "brain noise (%d sources): %6.1f fT rms\n", NUMNOISE, 1.0e+15*rms);
    fprintf(flog, "brain noise density: %5.1f fT/rt-Hz\n", 1.0e+15*rms/sqrt(BANDWIDTH));
    fflush(flog);

    // simulate sensor noise
    if (vflg) {
        printf("\nadding simulated sensor noise time-series");
        fflush(stdout);
    }
    for (m=0; m<M; m++)         // for each sensor...
        for (t=0; t<T; t++)     // generate Gaussian random time-series with zero mean & standard deviation of CTFNoise
            data[m][t] += RandomGauss(0., CTFNoise);

    // read dipoles
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
    Dipole = new_arrayE(VOXELINFO, D, "Dipole[D]");
    if (fflg) {
        interval = edge / (sqrt((double)D) - 1.);
        Dx = -edge / 2.;
        Dy = -edge / 2.;
    }

    if (sflg) {
        printf("\nRe-setting  user selected seed to generate random test dipoles\n");
        srand48(seed);   // seed random number generator
    }

    for (d=0, MaxDipole=0.; d<D; d++) {
        if(fflg) {
            Vc[X_] = Dx;
            Vc[Y_] = Dy;
            Vc[Z_] = OUTERRAD-CorticalDepth;
        } else {
            Vc[X_] = drand48() * edge - edge / 2.;
            Vc[Y_] = drand48() * edge - edge / 2.;
            Vc[Z_] = OUTERRAD-CorticalDepth;
        }
        CtoS(Vc, Vs);                   // convert dipole coordinates to spherical coordinates
        Vs[RHO_] = OUTERRAD-CorticalDepth;
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

        // get maximum z-coordinate of test dipole
        if (Dipole[d].p[Z_] > MaxDipole)
            MaxDipole = Dipole[d].p[Z_];
    }

    // output sensor & dipole positions to file
    if (vflg) {
        printf(" - done\n");
        printf("gap between sensor & dipole arrays is %4.1f mm\n", 1000.*(MaxSensor-MaxDipole));
	    printf("writing sensor & dipole geometry");
	    fflush(stdout);
    }
    if (zflg)
      sprintf(output, "%s_Geometry.txt", prefix);
    else
      sprintf(output, "CTF_Geometry,a%d,d%d,e%4.1f.txt", M, D, 1000.*edge);
    if((fp = fopen(output, "w")) == NULL)
        Cleanup("can't open geometry file for write");
    for (m=0; m<M; m++) {
        mm = Header.PsIndex[m];
        fprintf(fp, "%7.3f\t%7.3f\t%7.3f\n",
          1000.*Channel[mm].Geom.MEGSensor.origin.p[X_], 1000.*Channel[mm].Geom.MEGSensor.origin.p[Y_], 1000.*Channel[mm].Geom.MEGSensor.origin.p[Z_]);
    }
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

	// simulate dipole amplitudes
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
		    mm = Header.PsIndex[m];
		    forward[d][m] = FwdSensor(&Q6, &Channel[mm].Geom.MEGSensor, ORDER);
		}

		// sum this dipole MEG signal after scaling by moment
		for (m=0; m<M; m++)
			for (t=0; t<T; t++)
			    data[m][t] += dip[d][t] * forward[d][m];

		// heartbeat
		if (vflg) {
		    printf(".");
		    fflush(stdout);
		}
	}

    // compute spatial angle matrix
    if (vflg) {
        printf(" - done\n");
        printf("computing spatial angle matrix");
        fflush(stdout);
    }
    if (zflg)
      sprintf(output, "%s_SpatialAngles.txt", prefix);
    else
      sprintf(output, "CTF_SpatialAngles,a%d,d%d,e%4.1f.txt", M, D, 1000.*edge);
    if((fp = fopen(output, "w")) == NULL)
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
	for (i=0; i<M; i++)
		for (j=0; j<M; j++) {
			for (t=0, tmp=0.; t<T; t++)
			    tmp += data[i][t] * data[j][t];
            gsl_matrix_set(C, i, j, tmp / (double)T);
		}
    pinv(C, Cinv);

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
        DIPSolve2D(&Dipole[d], Bp, Cinv, &Header, Channel, 0, &span);
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
           printf("writing SAM image file\n");
           fflush(stdout);
        }
        range = edge / 2. + 0.005;
        tmp = 1000. * edge + 11.;
        tmp *= tmp;
        if (zflg)
          sprintf(output, "%s_Image.txt", prefix);
        else
          sprintf(output, "CTF_Image,a%d,d%d,e%4.1f.txt", M, D, 1000.*edge);
        if ((fp = fopen(output, "w")) == NULL)
          Cleanup("can't open 'CTF_Image.txt' file for write");
        Vc[Z_] = OUTERRAD-CorticalDepth;
        for (Vc[X_]=(-range), i=0, pct=(-1); Vc[X_]<=(range+IOTA); Vc[X_]+=0.001) {
          for (Vc[Y_]=(-range); Vc[Y_]<=(range+IOTA); Vc[Y_]+=0.001, i++) {
            CtoS(Vc, Vs);           // convert to spherical coordinates
            Vs[RHO_] = OUTERRAD-CorticalDepth;    // set image radius
            StoC(Vs, Vb);           // convert to cartesian coordinates
            for (v=X_; v<=Z_; v++) {
                Dipole[0].p[v] = Vb[v];
                Dipole[0].v[v] = 0.;
            }
            Dipole[0].Solve = TRUE;
            Dipole[0].ROI = TRUE;

            // compute source z^2 & save to file
            DIPSolve2D(&Dipole[0], Bp, Cinv, &Header, Channel, 0, &span);
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

    gsl_matrix_free(C);
    gsl_matrix_free(Cinv);
    gsl_vector_free(Bp);
    gsl_vector_free(W);
    free(Dipole);

    // open & write gain file
    if((fp = fopen("CTF_Gain.txt", "w")) == NULL)
        Cleanup("can't open 'CTF_Gain.txt' file for write");
    for (d=0; d<D; d++)
        fprintf(fp, "%8.6f\n", Gain[d]);
    fclose(fp);
    free(Gain);

    // open 'Solutions' file
    if (zflg)
      sprintf(output, "%s_Matrix.txt", prefix);
    else
      sprintf(output, "CTF_Matrix,a%d,d%d,e%4.1f.txt", M, D, 1000.*edge);
    if((fp = fopen(output, "w")) == NULL)
        Cleanup("can't open 'CTF_Matrix.txt' file for write");

    // compute correlation coefficient of SAM estimate with original dipole time-series
    if (vflg) {
        printf(" - done\n");
        printf("computing correlation of beamformer dipole estimates with original dipoles\n");
        fflush(stdout);
    }
    correlation = new_matrixE(D, D, "correlation[D][D]");
    for (d=0, pct=(-1); d<D; d++)   // dipole index
        for (dd=0; dd<D; dd++) {    // source index
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
    fclose(fp);
    free_matrix(src);
    free_matrix(dip);

    // open & write isr file
    if (zflg)
      sprintf(output, "%s_ISR.txt", prefix);
    else
      sprintf(output, "CTF_ISR,a%d,d%d,e%4.1f.txt", M, D, 1000.*edge);
    if((fp = fopen(output, "w")) == NULL)
        Cleanup("can't open 'CTF_ISR.txt' file for write");
    for (d=0; d<D; d++)
        for (dd=0; dd<D; dd++)
            fprintf(fp, "%9.6f\t%9.6f\n", 1000.*distance[d][dd], fabs(correlation[d][dd]));
    fclose(fp);

    // compute independent source resolution based on mean diagonal minus mean off-diagonal correlations
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
    free_matrix(correlation);

    // that's all, folks!
    fprintf(flog, "ISR = %f\n", isr);
    fprintf(flog, "residual variance = %f percent\n", rvar);
    fprintf(flog, "gap = %4.1f mm\n", 1000.*(MaxSensor-MaxDipole));
    fclose(flog);
    if (vflg) {
        printf(" - done\n\n");
        printf("CTF simulation done\n");
        fflush(stdout);
    }
    exit(0);
}
