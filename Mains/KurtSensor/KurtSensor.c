// KurtSensor -- computes the excess kurtosis of all sensors
//
//	flags:
//		-r <dataset name>			-- MEG run name
//		-m <parameter file name>	-- parameters
//		-v							-- verbose mode
//
//		Author: Stephen E. Robinson
//				MEG Core Facility
//				NIMH
//


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <string.h>
#include <fcntl.h>
#include <geoms.h>
#include <samlib.h>
#include <siglib.h>
#include <filters.h>
#include <rvelib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <samutil.h>
#include <filters.h>
#include <version.h>

#define MINOR_REV		0
#define TRUE			1
#define FALSE			0
#define SIG_ORDER		4


int	main(
	int		argc,
	char	**argv
)

{
	FFTSPEC			FFTData;				// FFT data filter structure
	FFTSPEC			FFTNotch;				// FFT notch filter structure (sum of notches)
	HeaderInfo		Header;					// MEG data header
	ChannelInfo		*Channel;				// MEG channel info
	EpochInfo		*Epoch;					// MEG epoch info
	double			*x;						// x[TT] -- primary sensor data
	double          *In;                    // In[TT] -- input time-series
	double			*Out;					// Out[TT] -- output time-series
	double          *Seg;                   // Seg[TW]
	double			p2;						// sample^2
	double			p4;						// sample^4
	double			sx2;					// sum of sample^2
	double			sx4;					// sum of sample^4
	double			var;					// variance
	double			g2;						// excess kurtosis
	double          max;                    // maximum kurtosis
	double          HPfreq;                 // highpass edge
	double          LPfreq;                 // lowpass edge
	double			BWData;					// data filter bandwidth
	double          SegTime;                // segment duration
	double			tmp;					// temporary
	int				E;						// total number of epochs
	int             M;                      // number of MEG or reference channels used
	int             T;                      // number of samples
	int				TT;						// E*T
	int             c;                      // input character
	int             e;                      // epoch index
	int             m;                      // channel index
	int             n;                      // trigger marker index
	int             t;                      // sample index
	int				tt;						// alternate sample index
	int             NumSegs;                // number of segments
	int             TW;                     // samples per segment
	int             Step;                   // half segment step
	int             w;                      // segment index
	int             LongIndex;              // position of arguments
	static int      eflg = FALSE;           // command-line error flag
	static int		rflg = FALSE;			// dataset name flag
	static int      mflg = FALSE;           // parameter file name flag
	static int      vflg = FALSE;           // verbose mode flag
	extern char     *optarg;
	extern int      opterr;
	char            fpath[256];             // general path name
	char            DSName[256];            // MEG dataset name
	char            DSpath[256];            // MEG dataset path
	char            SAMpath[256];           // SAM subdirectory path
	char			ParmName[256];			// parameter file name
	char			Name[8];				// channel name
    char            Line[256];              // line buffer
    char            Label[32];              // parameter label
	double          scale;                  // bits-->Tesla conversion
	unsigned char   **Bad;                  // Bad[E][T] -- flags for bad samples
	char            ShortOpts[] = "r:m:v";
	static struct option    LongOpts[] = {  // command line options
		{"dataset_name", required_argument, NULL, 'r'},
		{"analysis_name", required_argument, NULL, 'm'},
		{"verbose", no_argument, NULL, 'v'},
		{NULL, no_argument, NULL, 0}
	};
	FILE            *fp;

	// parse command line parameters
	opterr = 1;             // enable error reporting
	while ((c = getopt_long(argc, argv, ShortOpts, LongOpts, &LongIndex)) != EOF) {
		switch(c) {
			case 'r':       // run name
				sprintf(DSName, "%s", optarg);
				rflg = TRUE;
				break;
			case 'm':       // analysis parameter file name
				sprintf(ParmName, "%s", optarg);
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
	if (eflg || !rflg || !mflg) {
		fprintf(stderr, "KurtSensor\t-r <run name>\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
		fprintf(stderr, "\t-m <parameter file name>\n");
		fprintf(stderr, "\t-v -- verbose mode\n");
		exit(-1);
	}

	// get dataset information structures
	if (vflg) {
		printf("S E N S O R - S P A C E   E X C E S S   K U R T O S I S\t%s\n", PRG_REV);
		printf("opening data file");
		fflush(stdout);
	}

	// get data with sensor structures
	GetDsInfo(DSName, &Header, &Channel, &Epoch, &Bad, TRUE);
	sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
	sprintf(SAMpath, "%s/SAM", DSpath);

	// extract constants
	E = Header.NumEpochs;
	M = Header.NumPri;
	T = Header.MaxSamples;
	TT = E * T;					// total samples

	// parse analysis parameter specification file in current working directory
	if (vflg) {
		printf(" - done\n");
		printf("parsing '%s' parameter file", ParmName);
		fflush(stdout);
	}
    sprintf(fpath, "%s", ParmName);
    if ((fp = fopen(fpath, "r")) == NULL)
        Cleanup("can't open simulation parameter file");
    while (fgets(Line, 128, fp) != NULL) {              // scan parameter values
        if (strlen(Line) == 1)                          // skip empty lines
            continue;
        if (Line[0] == '#')                             // skip comment lines	
            continue;
        if (sscanf(Line, "%s", Label) != 1)             // read the label field
            Cleanup("can't read label field");
        if (!strncmp("BandPass", Label, 8)) {           // read 'BandPass'
            if (sscanf(Line, "%s%lf%lf", Label, &HPfreq, &LPfreq) != 3)
                Cleanup("can't read 'BandPass' frequencies");
         //   FilterType = BANDPASS;
        } else if (!strncmp("SegTime", Label, 7)) {     // read 'SegTime'
            if (sscanf(Line, "%s%lf", Label, &SegTime) != 2)
                Cleanup("can't read 'SegTime'");
        } else {
            fprintf(stderr, "\nBad label: '%s'\n", Label);
            fflush(stdout);
        }
    }
    fclose(fp);

    // determine segmentation
    TW = (int)floor(SegTime * Header.SampleRate);   // floor!
    Step = TW / 2;
    NumSegs = T / Step - 1;         // for 50% overlap

	// allocate memory for arrays
	if (vflg) {
		printf(" - done\n");
		printf("allocating memory");
		fflush(stdout);
	}
	x = new_arrayE(double, TT, "x[TT]");
	In = new_arrayE(double, TT, "In[T]");
	Out = new_arrayE(double, TT, "Out[T]");
	Seg = new_arrayE(double, TW, "Seg[TW]");

	// compute data filter
	if (vflg) {
		printf(" - done\n");
		printf("computing %6.3f to %6.3f Hz data & notch filter", HPfreq, LPfreq);
		fflush(stdout);
	}
	Butterworth(&FFTData, T, HPfreq, LPfreq, Header.SampleRate, BANDPASS, MAXORDER, &BWData);
    mknotch(&FFTNotch, HPfreq, LPfreq, Header.SampleRate, T, 60.);  

	// open file for saving kurtosis
	sprintf(fpath, "%s,Kurt.txt", DSName);
	if((fp = fopen(fpath, "w")) == NULL)
		Cleanup("can't open kurtosis file for write");

	// for each sensor...
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("computing excess kurtosis for all primary sensors:\n");
		fflush(stdout);
	}
	for(m=0; m<M; m++) {

		// get channel number & name
		n = Header.PsIndex[m];
		strncpy(Name, Channel[n].ChannelName, 5);
		if(vflg == TRUE) {
			printf("\t%s", Name);
			fflush(stdout);
		}

		// read & filter channel data
		scale = Channel[n].Scale;
		for(e=0; e<E; e++) {;
			GetDsData(&Header, e, n, scale, In);
			FFTfilter(In, Out, &FFTData, T);
			FFTfilter(Out, In, &FFTNotch, T);
			demean(In,T);
			for(t=0, tt=e*T; t<T; t++, tt++)
				x[tt] = In[t];
		}

        // parse data segments
        for (w=1, max=0.; w<NumSegs-1; w++) {       // ignore starting & ending segments!
            for (t=0, tt=w*Step; t<TW; t++, tt++)
                Seg[t] = In[tt];
         
            // compute excess kurtosis
            demean(Seg, TW);
            for (t=0,sx2=sx4=0.; t<TW; t++) {
                tmp = Seg[t];
			    p2 = tmp * tmp;
			    p4 = p2 * p2;
			    sx2 += p2;
			    sx4 += p4;
		    }
		    var = sx2 / (double)TW;
		    g2 = sx4 / ((double)TW * var * var) - 3.;
		    if (g2 > max)
		        max = g2;
		}

		// log excess kurtosis
		fprintf(fp, "%d\t%5s\t%f\n", m, Name, max);
		fflush(fp);

		// done with this sensor
		if(vflg == TRUE) {
			printf(" - done\n");
			fflush(stdout);
		}
	}			// for(m...
	fclose(fp);

	// we're done now!
	if(vflg == TRUE) {
		printf("'KurtSensor' done\n");
		fflush(stdout);
	}
	exit(0);
}
