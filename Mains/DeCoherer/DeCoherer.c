// DeCoherer -- adaptive spectral noise reduction (related to the optimal Wiener filter)
//
//  flags:
//      -r <dataset name>	-- MEG run name
//      -l <list name>		-- file containing list of reference channels
//      -s <segment duration>	-- time duration of data segments (optional for multi-epoch data)
//      -d			-- delete backup file (optional)
//      -v			-- verbose mode
//
//  Procedures:
//      Principle is analogous to using multiple regression to find the
//      gains required for subtraction of correlated signals. Here one
//      performs these steps in the frequency domain -- one spectral line
//      at a time (after J. Vrba). The trickery part (shhhhhh...!) is that
//      we do this _without_ applying a window function to the time-series.
//
//  Algorithm:
//      Let us assume primary sensor time-series x[t] & multiple reference sensor time-series y[i][t].
//      After FFT of these, we obtain complex spectra X[s] & Y[i][s]. For each spectral line, we
//      compute the regression of X on Y[i]. The regression slope is the complex "gain" for subtraction
//      of the reference. This is given by:
//
//      Gain[s] = Cxy[i][s] / Ayy[i][s]     (note that gain is type fftw_complex)
//
//      where Cxy[i][s] is the crosspower of X[s] times the complex conjugate of Y[i][s], &
//      Ayy[i][s] is the autopower of Y[i][s]. Note that the crosspower & autopower are statistics.
//      This means they must be accumulated across multiple data segments or trials. The reference
//      spectral line is subtracted from the primary sensor spectrum -- on a segment-by-segment basis:
//
//      X[s] = X[s] - Gain[s] * Y[i][s]
//
//      The order in which references are subtracted is determined by the primary-reference coherence:
//
//      Rxy[i][s] = (Cxy[i][s] * Cxy[i][s]*) / (Axx[s] * Ayy[i][s])
//
//      where Axx[s] is the autopower of the primary sensor. The reference sensor having the highest
//      coherence is called the "pivot" reference. Coherence is analogous to the correlation
//      coefficient, but is always positive.
//
//      After each pivot reference is subtracted from the primary sensor, the remaining references are
//      "deflated" by subtracting the pivot reference spectrum from each remaining reference. This is
//      achieved using an identical procedure as with the primary sensor. After reference deflation.
//      there is one fewer reference, and the next pivot reference is found and subtracted from the
//      primary sensor. Pivot reference channels are marked "used" so that they do not enter into
//      subsequent calculations for the remaining references. This process is, of course, row reduction!
//
//      The "cleaned" primary sensor spectrum is then converted back to the time domain, & written
//      back to the data file.
//
//      In revision 1, the statistics and gains are computed from windowed FFT segments. The
//      complex gains for each spectral line are applied to the unwindowed FFTs of the primary
//      & reference sensors, by simple liinear interpolatiion.
//
//      Author: Stephen E. Robinson
//              MEG Core Facility
//              NIMH
//
//      (after a suggestion by J. Vrba, & based on code from CTFI
//

#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <fcntl.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <byteswap.h>
#include <samlib.h>
#include <DataFiles.h>
#include <complex.h>
#include <fftw3.h>
#include <samutil.h>
#include <version.h>

#define MINOR_REV   1       // initial revision
#define CX0     (0.+0.*I)   // complex zero


// Let us begin...
main(
    int     argc,
    char    **argv
)

{
    HeaderInfo      Header;         // header information data
    ChannelInfo     *Channel;       // Channel[M] -- channel information data
    EpochInfo       *Epoch;         // Epoch[E] -- epoch information data
    fftw_plan       XUFwd;          // plan for primary unwindowed forward FFT
    fftw_plan       XUInv;          // plan for primary inverse FFT
    fftw_plan       YUFwd;          // plan for reference unwindowed forward FFT
    fftw_plan       XSegFwd;        // plan for primary forward FFT
    fftw_plan       YSegFwd;        // plan for refence forward FFT
    fftw_complex    *X;             // X[TW] -- primary sensor spectrum
    fftw_complex    *Y;             // Y[TW/2+1] -- reference sensor spectrum
    fftw_complex    **Xo;           // Xo[G][S] -- primary sensor windowed spectrum
    fftw_complex    *Xu;            // Xu[SS] -- primary sensor un-windowed spectrum
    fftw_complex    *Yfwd;          // Yfwd[S] -- reference sensor un-windowed spectrum
    fftw_complex    **Yu;           // Yu[R][SS] -- reference sensor un-windowed spectrum
    fftw_complex    *Xs;            // Xs[G] -- primary sensor single line
    fftw_complex    **XYgain;       // XYgain[R][S] -- gain for subtraction of unsegmented data
    fftw_complex    ***Yo;          // Yo[R][G][S] -- reference sensor windowed spectra
    fftw_complex    **Ys;           // Ys[R][G] -- reference sensor single lines
    fftw_complex    *Cxy;           // Cxy[R] -- primary-reference crosspowers
    fftw_complex    Cyy;            // pivot_reference-reference crosspower
    fftw_complex    Gain;           // complex gain
    double          *Seg;           // Seg[TW] -- data segment time series
    double          *Ayy;           // Ayy[R] -- reference autopowers
    double          Axx;            // primary autopower
    double          Rxy;            // coherence
    double          *In;            // In[T] -- single channel time-series
    double          *Out;           // Out[T] -- single channel time-series
    double          *Hann;          // Hann[TW] -- Hanning window
    double          max;            // maximum coherence (or maximum signal)
    double          SegTime;        // segmentation window time
    double          binHz;          // frequency resolution
    double          scale;          // bits-->Tesla
    double          fract;          // fraction for interpolating spectrum
    int32_t         *xint;          // xint[T] -- one channel, epoch time-series data
    int32_t         tmp32;          // for byteswap
    int             *RList;         // RList[R] -- reference channel index list
    int             *Reference;     // Reference[R] -- reference pivot usage list
    int             c;              // command line flag variable
    int             M;              // number of channels
    int             G;              // segments
    int             T;              // file samples/channel/epoch
    int             TW;             // samples/time-window
    int             Step;           // step for overlapping windows
    int             R;              // number of reference sensors
    int             P;              // number of primary sensors
    int             S;              // number of spectral lines in segment
    int             SS;             // number of spectral lines for all data samples
    int             g;              // segment index
    int             i;              // general purpose index
    int             j;              // another crucial index
    int             p;              // primary index
    int             r;              // reference index
    int             s;              // spectral bin index
    int             ss;             // alternate spectral bin index
    int             t;              // sample index
    int             tt;             // alternate sample index
    int             m;              // sensor index
    int             n;              // pivot loop index
    int             pivot;          // pivot index
    int             fdin;           // input file descriptor
    int             fdout;          // output file descriptor
    int             fdtmp;          // temporary file descriptor
    int             NumBytes;       // data bytes
    int             found;          // flag for finding named reference index
    int             bflg = FALSE;   // backup flag
    int             dflg = FALSE;   // dataset flag
    int             eflg = FALSE;   // command-line error flag
    int             lflg = FALSE;   // reference list flag
    int             sflg = FALSE;   // segmentation flag
    int             vflg = FALSE;   // verbose mode flag
    extern char     *optarg;
    extern int      optind;
    extern int      opterr;
    char            fpath[256];     // general path name
    char            SetName[256];   // MEG/OPM dataset name
    char            TmpName[256];   // temporary file name
    char            ListName[256];  // reference list file name
    char            RefName[32];    // reference channel name
    char            OldName[256];   // old '.meg4' file name
    char            NewName[256];   // new '.meg4' file name
    char            MEG4Hdr[] =     // MEG 4.1 header string
        "MEG41CP";
    FILE            *fp;            // general file pointer
    off_t           offset;         // number of bytes offset

    // parse command line parameters
    opterr = 1;     // enable error reporting
    while ((c = getopt(argc, argv, "r:l:s:dv")) != EOF) {
        switch(c) {
            case 'r':   // dataset name (mandatory)
                sprintf(SetName, "%s", optarg);
                dflg = TRUE;
                break;
            case 'l':   // reference list (mandatory)
                sprintf(ListName, "%s", optarg);
                lflg = TRUE;
                break;
            case 's':   // duration of data segments (seconds)
                SegTime = atof(optarg);
                sflg = TRUE;
                break;
            case 'd':   // delete backup
                bflg = TRUE;
                break;
            case 'v':   // verbose
                vflg = TRUE;
                break;
            default:
                eflg = TRUE;
                break;
        }
    }
    if (eflg || !dflg  || !lflg || !sflg) {
        fprintf(stderr, "DeCoherer [arguments]\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
        fprintf(stderr, "\t-r <dataset name>\n");
        fprintf(stderr, "\t-l <reference list file>\n");
        fprintf(stderr, "\t-s <segment duration (s)>\n");
        fprintf(stderr, "\t-d -- delete backup file\n");
        fprintf(stderr, "\t-v -- verbose mode\n");
        exit(-1);
    }

    // announce & get data with sensor structures
    if (vflg) {
        printf("S P E C T R A L   A D A P T I V E   N O I S E   R E D U C T I O N\t%s%0d\n", PRG_REV, MINOR_REV);
        printf("opening dataset");
        fflush(stdout);
    }
    GetDsInfo(SetName, &Header, &Channel, &Epoch, NULL, TRUE);
    if (Header.NumEpochs > 1)
        Cleanup("single continuous dataset req'd");
    M = Header.NumChannels;     // total number of channels in dataset
    P = Header.NumPri;
    R = Header.NumRef;          // to be updated after reading reference list
    T = Header.MaxSamples;      // number of samples per channel per epoch
    TW = (int)floor(SegTime * Header.SampleRate);   // floor!
    Step = TW / 2;              // overlap for windows
    G = T / Step - 1;           // number of overlapping segments for 50% overlap
    if (G < 10)
        Cleanup("insufficent segments -- decrease segment time");
    S = TW / 2 + 1;             // number of frequency bins in segment spectrum
    SS = T / 2 + 1;             // number of spectral lines for all samples
    binHz = Header.SampleRate / (double)TW;

    // read reference list file
    if (vflg) {
        printf(" - done\n");
        printf("reading reference list file");
        fflush(stdout);
    }
    if ((fp = fopen(ListName, "r")) == NULL)
        Cleanup("can't open reference list file for read");
    if (fscanf(fp, "%d", &i) != 1)
        Cleanup("can't read number of reference channels");
    if (i > R)
        Cleanup("reference list file contained more than 32 channels");
    R = i;          // update reference count
    RList = new_arrayE(int, R, "RList[R]");
    for (r=0; r<R; r++) {
        if (fscanf(fp, "%s", RefName) != 1)
            Cleanup("can't read reference channel name");
        for (j=0, found=FALSE; j<M; j++)
            if (!strcmp(Channel[j].ChannelName, RefName)) {
                RList[r] = j;
                found = TRUE;
                break;
            }
        if (!found)
            Cleanup("could not find listed reference channel");
    }
    fclose(fp);

    // summarize analysis
    if (vflg) {
        printf(" - done\n");
        printf("cleaning dataset '%s', using parameters:\n", SetName);
        printf("\tM: %d total channels\n", M);
        printf("\tP: %d primary channels\n", P);
        printf("\tR: %d reference channels\n", R);
        printf("\tT: %d samples/epoch/channel\n", T);
        printf("\tG: %d total data segments\n", G);
        printf("\tU: %d samples/segment\n", TW);
        printf("\tS: %d segment step\n", Step);
        printf("\tF: %6.1f Hz sample rate\n", Header.SampleRate);
        printf("\tB: %9.6f Hz frequency resolution\n", binHz);
        fflush(stdout);
    }
        
    // determine array sizes & allocate memory for arrays
    if (vflg) {
        printf(" - done\n");
        printf("allocating memory");
        fflush(stdout);
    }
    In = new_arrayE(double, T, "In[T]");
    Out = new_arrayE(double, T, "Out[T]");
    Seg = new_arrayE(double, TW, "Seg[TW]");
    Hann = new_arrayE(double, TW, "Hann[TW]");
    X = (fftw_complex *)fftw_malloc(S * sizeof(fftw_complex));
    Y = (fftw_complex *)fftw_malloc(S * sizeof(fftw_complex));
    Xo = (fftw_complex **)fftw_malloc(G * sizeof(fftw_complex *));
    for (g=0; g<G; g++)
        Xo[g] = (fftw_complex *)fftw_malloc(S * sizeof(fftw_complex));
    Xu = (fftw_complex *)fftw_malloc(SS * sizeof(fftw_complex));
    Yfwd = (fftw_complex *)fftw_malloc(SS * sizeof(fftw_complex));
    Yu = (fftw_complex **)fftw_malloc(R * sizeof(fftw_complex *));
    for (r=0; r<R; r++)
        Yu[r] = (fftw_complex *)fftw_malloc(SS * sizeof(fftw_complex));
    Xs = (fftw_complex *)fftw_malloc(G * sizeof(fftw_complex));
    Yo = (fftw_complex ***)fftw_malloc(R * sizeof(fftw_complex **));
    for (r=0; r<R; r++) {
        Yo[r] = (fftw_complex **)fftw_malloc(G * sizeof(fftw_complex *));
        for (g=0; g<G; g++)
            Yo[r][g] = (fftw_complex *)fftw_malloc(S * sizeof(fftw_complex));
    }
    Ys = (fftw_complex **)fftw_malloc(R * sizeof(fftw_complex *));
    for (r=0; r<R; r++)
        Ys[r] = (fftw_complex *)fftw_malloc(G * sizeof(fftw_complex));
    XYgain = (fftw_complex **)fftw_malloc(R * sizeof(fftw_complex *));
    for (r=0; r<R; r++)
        XYgain[r] = (fftw_complex *)fftw_malloc(S * sizeof(fftw_complex));
    Cxy = (fftw_complex *)fftw_malloc(R * sizeof(fftw_complex));
    Ayy = new_arrayE(double, R, "Ayy[R]");
    xint = new_arrayE(int32_t, T, "xint[T]");
    Reference = new_arrayE(int, R, "Reference[R]");

    // create FFT plans
    if (vflg) {
        printf(" - done\n");
        printf("initializing forward & inverse FFT transforms");
        fflush(stdout);
    }
    XSegFwd = fftw_plan_dft_r2c_1d(TW, Seg, X, FFTW_ESTIMATE);
    YSegFwd = fftw_plan_dft_r2c_1d(TW, Seg, Y, FFTW_ESTIMATE);
    XUFwd = fftw_plan_dft_r2c_1d(T, Out, Xu, FFTW_ESTIMATE);
    XUInv = fftw_plan_dft_c2r_1d(T, Xu, Out, FFTW_ESTIMATE);
    YUFwd = fftw_plan_dft_r2c_1d(T, Out, Yfwd, FFTW_ESTIMATE);

    // generate Hanning window
    Hanning(Hann, TW);

    // compute FFT of reference sensors
    if (vflg) {
        printf(" - done\n");
        printf("processing reference sensors:\n");
        fflush(stdout);
    }
    for (r=0; r<R; r++) {
        m = RList[r];
        scale = Channel[m].Scale;
        GetDsData(&Header, 0, m, scale, In);
        demean(In, T);

        // compute unwindowed FFT for entire time-series
        for (t=0; t<T; t++)     // save In, so it is not over-written
            Out[t] = In[t];
        demean(Out, T);
        fftw_execute(YUFwd);
        for (ss=0; ss<SS; ss++)
            Yu[r][ss] = Yfwd[ss];

        // for each overlapping window
        for (g=0; g<G; g++) {

            // compute windowed FFT & copy to spectral array
            for (t=0, tt=g*Step; t<TW; t++, tt++)
                Seg[t] = Hann[t] * In[tt];
            fftw_execute(YSegFwd);     // compute reference channel complex spectrum
            for (s=0; s<S; s++)
                Yo[r][g][s] = Y[s];
        }
        if (vflg) {
            printf(" - done\n");
            fflush(stdout);
        }
    }

    // do coherence balancing
    if (vflg) {
        printf("processing primary sensors:\n");
        fflush(stdout);
    }
    for (p=0; p<P; p++) {
        m = Header.PsIndex[p];
        if (vflg) {
            printf("\tprimary sensor '%s' - reading/FFT", Channel[m].ChannelName);
            fflush(stdout);
        }
        scale = Channel[m].Scale;
        GetDsData(&Header, 0, m, scale, In);
        demean(In, T);

        // compute unwindowed FFT for entire time-series
        for (t=0; t<T; t++)     // save In, so it is not over-written
            Out[t] = In[t];
        demean(Out, T);
        fftw_execute(XUFwd);

        // for each segment
        for (g=0; g<G; g++) {

            // compute windowed FFT for accumulating statistics
            for (t=0, tt=g*Step; t<TW; t++, tt++)
                Seg[t] = Hann[t] * In[tt];
            fftw_execute(XSegFwd);                      // compute primary channel complex spectrum
            for (s=0; s<S; s++)
                Xo[g][s] = X[s];
        }

        // for each spectral line...
        if (vflg) {
            printf(" - cleaning");
            fflush(stdout);
        }
        for (s=0; s<S; s++) {

            // make copies of primary & reference's single spectral line
            for (g=0; g<G; g++) {
                Xs[g] = Xo[g][s];
                for (r=0; r<R; r++)
                    Ys[r][g] = Yo[r][g][s];
            }

            // compute primary autopower
            for (g=0, Axx=0.; g<G; g++)
                Axx += Xs[g] * conj(Xs[g]);
            Axx /= (double)G;                               // normalize primary autopower by segment count

            // initialize reference flags
            for (r=0; r<R; r++)
                Reference[r] = TRUE;

            for (n=0; n<R; n++) {	// deflate R times

                // compute primary-reference crosspower, reference autopowers, & determine pivot channel
                for (r=0, max=0.; r<R; r++)
                    if (Reference[r]) {

                        // compute primary-reference crosspower & reference autopower
                        for (g=0, Ayy[r]=0., Cxy[r]=CX0; g<G; g++) {
                            Cxy[r] += Xs[g] * conj(Ys[r][g]);
                            Ayy[r] += Ys[r][g] * conj(Ys[r][g]);
                        }
                        Cxy[r] /= (double)G;                // normalize primary-reference crosspower by segment count
                        Ayy[r] /= (double)G;                // normalize reference autopower by segment count

                        // compute primary-reference coherence & determine pivot from largest coherence
                        Rxy = (Cxy[r] * conj(Cxy[r])) / (Axx * Ayy[r]);
#ifdef DEBUG
                        printf("\ns=%d, r=%d, Rxy=%f\n", s, r, Rxy);
                        fflush(stdout);
#endif
                        if (Rxy > max) {
                            max = Rxy;
                            pivot = r;
                        }
                    }
                Reference[pivot] = FALSE;		// important: mark the pivot channel as used...

                // compute primary-pivot_reference complex gain
                XYgain[n][s] = Cxy[pivot] / Ayy[pivot];     // use 'n' -- deflation count
#ifdef DEBUG
                printf("XYgain[%d][%d]=%e %e\n", n, s, creal(XYgain[n][s]), cimag(XYgain[n][s]));
                fflush(stdout);
#endif
                // compute reference-pivot_reference crosspower & deflate all non-pivot references
                for (r=0; r<R; r++)
                    if (Reference[r]) {                     // if this reference has not yet been used
                        for (g=0, Cyy=CX0; g<G; g++)
                            Cyy += Ys[pivot][g] * conj(Ys[r][g]);
                        Cyy /= (double)G;
                        Gain = Cyy / Ayy[r];
                        for (g=0; g<G; g++)                 // deflate reference
                            Ys[r][g] -= Ys[pivot][g] * conj(Gain);
                    }
            }   // end pivot	

        }   // end spectral lines

        // interpolate spectral bins before subtraction
        for (ss=0; ss<SS; ss++) {
            i = ss / S;
            j = ss % S;
            fract = (double)(S - j) / (double)S;
            for (r=0; r<R; r++) {
                Gain = fract * XYgain[r][i] + (1. - fract) * XYgain[r][i+1];
                    Xu[ss] -= Yu[r][ss] * conj(Gain);
            }
        }

        // compute inverse FFT of cleaned primary sensors & normalize
        if (vflg) {
            printf(" - iFFT/writing");
            fflush(stdout);
        }
        fftw_execute(XUInv);
        for (t=0; t<T; t++)         // normalize by number of samples? -- according to FFTW docs, yes...
            Out[t] /= (double)T;
        demean(Out, T);

        // open temporary data file
        sprintf(fpath, "%s.ds/%s.tmp", SetName, Channel[m].ChannelName);
        if ((fdtmp = open(fpath, O_WRONLY | O_CREAT | O_TRUNC, 0666)) == -1)
            Cleanup("can't open '.tmp' file for write");
        NumBytes = T * sizeof(int32_t);     // this is number of bytes per channel per epoch

        // convert segment to 32-bit integer -- note that CTF format uses byte-reversed integers
        for (t=0; t<T; t++) {
            tmp32 = (int32_t)((double)Out[t] / Channel[m].Scale);
            xint[t] = bswap_32(tmp32);
        }

        // write normalized 32-bit integer output to disk
        if (write(fdtmp, (void *)xint, NumBytes) != NumBytes)
            Cleanup("can't write data to '.tmp' file");

        close(fdtmp);
        if (vflg) {
            printf(" - done\n");
            fflush(stdout);
        }
    }       // end processing primary sensors

    // open output '.meg4' file
    if (vflg) {
        printf("writing output file");
        fflush(stdout);
    }
    sprintf(fpath, "%s.ds/%s.meg4", SetName, SetName);
    if ((fdin = open(fpath, O_RDONLY)) == -1)
        Cleanup("can't open '.meg4' file for reaf");
    sprintf(fpath, "%s.ds/%s.meg4.NEW", SetName, SetName);
    if ((fdout = open(fpath, O_WRONLY | O_CREAT | O_TRUNC, 0666)) == -1)
        Cleanup("can't open '.meg4' file for write");
    if (write(fdout, MEG4Hdr, 8) != 8)
        Cleanup("can't write header to '.meg4' file");
    NumBytes = T * sizeof(int32_t);		// this is number of bytes per channel per epoch

    // write all channels in the same order as the original dataset
    for (m=0; m<M; m++) {

        // if this is a primary (cleaned) sensor
        if (Channel[m].ChannelType == TYPE_MEG) {

            // open & read temporary file
            sprintf(fpath, "%s.ds/%s.tmp", SetName, Channel[m].ChannelName);
            if ((fdtmp = open(fpath, O_RDONLY)) == -1)
                Cleanup("can't open '.tmp' file for write");
            if (read(fdtmp, (char *)xint, NumBytes) != NumBytes)
                Cleanup("can't read data from '.tmp' file");
            close(fdtmp);

            // write cleaned primary sensor to .meg4 file
            if (write(fdout, (char *)xint, NumBytes) != NumBytes)
                Cleanup("can't write data to '.meg4' file");  

        // simply copy channel to output
        } else {

            // seek & read input
            offset = (off_t)(NumBytes * m) + (off_t)8;
            if (lseek(fdin, offset, SEEK_SET) != offset)
                Cleanup("'lseek()' to channel failed");
            if (read(fdin, (char *)xint, NumBytes) != NumBytes)
                Cleanup("can't read data from 'meg4' file");

            // copy to output
            if (write(fdout, (char *)xint, NumBytes) != NumBytes)
                Cleanup("can't write data to '.meg4' file");
        }	

        // heartbeat
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }
    close(fdin);
    close(fdout);

    // delete temporary files
    if (vflg) {
        printf(" - done\n");
        printf("deleting temporary data files");
        fflush(stdout);
    }
    for (p=0; p<P; p++) {
        m = Header.PsIndex[p];
        sprintf(TmpName, "%s.ds/%s.tmp", SetName, Channel[m].ChannelName);
        if (unlink(TmpName) == -1)
            Cleanup("can't delete temporary data file");
    }

    // rename original '.meg4' file
    if (vflg) {
        printf(" - done\n");
        printf("moving & renaming '.meg4' file");
        fflush(stdout);
    }
    sprintf(OldName, "%s.ds/%s.meg4", SetName, SetName);
    sprintf(NewName, "%s.ds/%s.meg4.bak", SetName, SetName);
    if (rename(OldName, NewName) == -1)
        Cleanup("backup of original '.meg4' file failed");
    sprintf(OldName, "%s.ds/%s.meg4.NEW", SetName, SetName);
    sprintf(NewName, "%s.ds/%s.meg4", SetName, SetName);
    if (rename(OldName, NewName) == -1)
        Cleanup("rename of filtered '.meg4' file failed");
    if (bflg) {
        if (vflg) {
            printf(" - done\n");
            printf("deleting backup (original) '.meg4' file");
            fflush(stdout);
        }
        sprintf(OldName, "%s.ds/%s.meg4.bak", SetName, SetName);
        if (unlink(OldName) == -1)
            Cleanup("can't delete backup data file");
    }
    if (vflg) {
        printf(" - done\n");
        printf("DeCoherer: done\n");
        fflush(stdout);
    }
    exit(0);
}
