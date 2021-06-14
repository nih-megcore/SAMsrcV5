// ArtRedux -- minimize artifacts using independent reference channels.
//  This program ambitiously uses multiple-regression to allow
//  multiple simultaneous artifact references.
//
//      y[i] = b[1] * x[i][1] + b[2] * x[i][2] + ... +  b[n] * x[i][n]
//
//       n   p                              n
//      Sum Sum x[i][j] * x[i][k] * b[k] = Sum x[i][j] * y[i]       j = 1...p
//      i=1 k=1                            i=1
//
//  in matrix form: (X'X)b = X'Y
//
//  solution:       b = (X'X)^-1 X'Y
//
//  The artifact file has the format:
//      number of artifact channels
//      <name of artifact channel 1> <highpass frequency> <lowpass frequency>
//          etc.
//
//  flags:
//      -r <dataset name>       -- MEG run name
//      -a <artifact list file> -- list of artifact reference channels & parameters
//      -d                      -- delete backup
//      -v                      -- verbose mode
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//
//          (based on code written while at CTF Systems, Inc.)
//

#include <stdlib.h>
#include <unistd.h>
#include <byteswap.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <fcntl.h>
#include <ctype.h>
#include <string.h>
#include <siglib.h>
#include <samlib.h>
#include <DataFiles.h>
#include <samutil.h>
#include <version.h>

#define MAINS       60.     // power mains frequency
#define FFTORDER    4       // 4th-order FFT filters
#define MINOR_REV   0

typedef struct	{           // artifact specification structure
    char    Name[16];       // artifact channel name
    int     Index;          // artifact channel index
    double  HPFreq;         // artifact high-pass frequency
    double  LPFreq;         // artifact low-pass frequency
    FFTSPEC	Filter;         // artifact FFT filter
} ARTSPECS;


main(
    int     argc,
    char    **argv
)

{
    ARTSPECS    *ArtSpec;       // ArtSpec[NA] -- artifact channel parameter list
    HeaderInfo  Header;         // header information data
    ChannelInfo *Channel;       // Channel[M] -- channel information data
    EpochInfo   *Epoch;         // Epoch[E] -- epoch information data
    FFTSPEC     NotchFFT;       // power mains notch filters
    gsl_matrix  *A;             // A -- X'X
    gsl_vector  *B;             // B -- X'Y
    gsl_vector  *C;             // C -- gain coefficient
    gsl_permutation *perm;      // LU permutation
    double      **Data;         // Data[M][T] -- complete data time-series
    double      *SosP;          // SosP[NA] -- sums of squares of primaries
    double      *SosA;          // SosA[NA] -- sums of squares of artifacts
    double      *In;            // In[T] --  filter input
    double      *Out;           // Out[T] -- filter output
    double      scale;          // bits-->Tesla scale factor
    double      BandWidth;      // filter bandwidth
    double      tmp;            // a universally temporary variable
    int         c;              // command line flag variable
    int         t;              // time index
    int         m;              // sensor index
    int         mp;             // primary channel index
    int         mr;             // reference channel index
    int         mi, mj;
    int         n;              // alternate sensor index
    int         p;              // primary channel index
    int         r;              // reference channel index
    int         i;              // matrix index
    int         j;              // another matrix index
    int         signum;         // sign of LU permutation
    int         P;              // number of primary channels used
    int         R;              // number of artifact channels used
    int         M;              // total number of channels used
    int         T;              // maximum number of samples
    int         FilterType;     // filter type
    int         flag;           // name-found flag
    int         NumBytes;       // data bytes
    int         fdin;           // input file descriptor
    int         fdout;          // output file descriptor
    int         fdtmp;          // temporary file descriptor
    int         aflg = FALSE;   // artifact file flag
    int         dflg = FALSE;   // delete backup flag
    int         eflg = FALSE;   // command line error flag
    int         rflg = FALSE;   // dataset flag
    int         vflg = FALSE;   // verbose mode flag
    int32_t     *xint;          // xint[T] -- channel data
    int32_t     tmp32;          // temporary 32-bit integer for byte-swap
    extern char *optarg;
    extern int  optind;
    extern int  opterr;
    char        fpath[256];     // general path name
    char        SetName[256];   // MEG dataset name
    char        ArtName[256];   // artifact list file name
    char        OldName[256];   // artifact list file name
    char        NewName[256];   // artifact list file name
    char        MEG4Hdr[] =     // MEG 4.1 header string
        "MEG41CP";
    FILE        *fp;            // file pointer for artifact list file
    off_t       offset;         // number of bytes offset

    // parse command line parameters
    opterr = 1;         // enable error reporting
    while ((c = getopt(argc, argv, "r:a:d:v")) != EOF) {
        switch (c) {
            case 'r':   // data set name (mandatory)
                sprintf(SetName, "%s", optarg);
                rflg = TRUE;
                break;
            case 'a':   // artifact list file name
                sprintf(ArtName, "%s", optarg);
                aflg = TRUE;
                break;
            case 'd':   // delete .meg4 backup
                dflg = TRUE;
                break;
            case 'v':   // verbose
                vflg = TRUE;
                break;
            default:
                eflg = TRUE;
                break;
        }
    }
    if (!rflg || !aflg || eflg) {
        fprintf(stderr, "usage: artRedux\t-r <dataset name>\t\t%s rev-%0d, %s\n", PRG_REV, MINOR_REV, __DATE__);
        fprintf(stderr, "\t\t-a <artifact channel list file>\n");
        fprintf(stderr, "\t\t-v -- verbose mode\n");
        exit(-1);
    }

    // announce & get data with sensor structures
    if (vflg) {
        printf("A D A P T I V E    A R T I F A C T    C A N C E L L A T I O N\t%s%0d\n", PRG_REV, MINOR_REV);
        printf("opening data file");
        fflush(stdout);
    }
    GetDsInfo(SetName, &Header, &Channel, &Epoch, NULL, TRUE);
    if (Header.NumEpochs > 1)
        Cleanup("single continuous dataset req'd");
    M = Header.NumChannels;     // total number of channels in dataset
    P = Header.NumPri;          // total number of primary channels (to be cleaned)
    R = Header.NumRef;          // to be updated after reading reference list
    T = Header.MaxSamples;      // number of samples per channel per epoch
    NumBytes = T * sizeof(int32_t); // number of bytes per channel

    // read artifact list file
    if (vflg) {
        printf(" - done\n");
        printf("reading artifact list file");
        fflush(stdout);
    }
    sprintf(fpath, "%s", ArtName);
    if ((fp = fopen(fpath, "r")) == NULL)
        Cleanup("can't open artifact channel list for read");
    if (fscanf(fp, "%d", &i) != 1)
        Cleanup("can't read number of artifact channels");
    if (i > R)
        Cleanup("too many channels in reference list file");
    R = i;          // update reference count
    ArtSpec = new_arrayE(ARTSPECS, R, "ArtSpec[R]");
    for (r=0; r<R; r++)
        if (fscanf(fp, "%s%lf%lf", &ArtSpec[r].Name, &ArtSpec[r].HPFreq, &ArtSpec[r].LPFreq) != 3)
            Cleanup("can't read artifact parameters");
    fclose(fp);

    // find reference channel indices
    if (vflg) {
        printf(" - done\n");
        printf("searching for reference channels");
        fflush(stdout);
    }
    for (r=0; r<R; r++) {
        n = strlen(ArtSpec[r].Name);
        for (m=0, flag=FALSE; m<Header.NumChannels; m++)
            if (!strncmp(Channel[m].ChannelName, ArtSpec[r].Name, n)) {
                ArtSpec[r].Index = m;
                flag = TRUE;
                break;
            }
        if (!flag)
           Cleanup("reference channel name not found");
    }

    // determine array sizes & allocate memory for arrays
    if (vflg) {
        printf(" - done\n");
        printf("determining array sizes");
        fflush(stdout);
    }
    if (vflg) {
        printf(" - done\n");
        printf("number of primary channels: %d, number of artifact channels: %d\n", P, R);
        printf("number of samples: %d\n", T);
        printf("allocating memory");
        fflush(stdout);
    }

    // allocate matrices & arrays
    if (vflg) {
        printf(" - done\n");
        printf("allocating arrays");
        fflush(stdout);
    }
    A = gsl_matrix_alloc(R, R);
    B = gsl_vector_alloc(R);
    C = gsl_vector_alloc(R);
    perm = gsl_permutation_alloc(R);
    xint = new_arrayE(int32_t, T, "xint[T]");
    Data = new_matrixE(M, T, "Data[M][T]");
    SosP = new_arrayE(double, P, "SosP[P]");
    SosA = new_arrayE(double, R, "SosA[R]");
    In = new_arrayE(double, T, "In[T]");
    Out = new_arrayE(double, T, "Out[T]");

    // generate filters
    if (vflg) {
        printf(" - done\n");
        printf("computing artifact filters:\n");
        fflush(stdout);
    }
    mknotch(&NotchFFT, 0., 175., Header.SampleRate, T, MAINS);
    for (r=0; r<R; r++) {
        if (vflg) {
            printf("artifact channel '%s'", ArtSpec[r].Name);
            fflush(stdout);
        }
        if (ArtSpec[r].HPFreq == 0. && ArtSpec[r].LPFreq == 0.) {
            if (vflg) {
                printf(" - AllPass\n");
                fflush(stdout);
            }
            FilterType = ALLPASS;
        } else if (ArtSpec[r].HPFreq > 0. && ArtSpec[r].LPFreq > 0.) {
            if (vflg) {
                printf(" - BandPass (%g to %g Hz)\n", ArtSpec[r].HPFreq, ArtSpec[r].LPFreq);
                fflush(stdout);
            }
            FilterType = BANDPASS;
        } else if (ArtSpec[r].HPFreq == 0. && ArtSpec[r].LPFreq > 0.) {
            if (vflg) {
                printf(" - LowPass (%g to %g Hz)\n", ArtSpec[r].HPFreq, ArtSpec[r].LPFreq);
                fflush(stdout);
            }
            FilterType = LOWPASS;
        } else if (ArtSpec[r].HPFreq > 0. && ArtSpec[r].LPFreq < Header.SampleRate/2.) {
            if (vflg) {
                printf(" - HighPass (%g to %g Hz)\n", ArtSpec[r].HPFreq, ArtSpec[r].LPFreq);
                fflush(stdout);
            }
            FilterType = HIGHPASS;
        } else {
            Cleanup("error in reference filter frequencies");
        }
        Butterworth(&ArtSpec[r].Filter, T, ArtSpec[r].HPFreq, ArtSpec[r].LPFreq, Header.SampleRate, FilterType, FFTORDER, &BandWidth);
	}

    // read all data
    if (vflg) {
        printf("reading data & applying %g Hz notch filter", MAINS);
        fflush(stdout);
    }
    for (m=0; m<M; m++) {
        scale = Channel[m].Scale;
        GetDsData(&Header, 0, m, scale, In);
        demean(In, T);
        for (r=0; r<R; r++)
        FFTfilter(In, Out, &NotchFFT, T);
        for (t=0; t<T; t++)
            Data[m][t] = Out[t];
        if (vflg) {
            printf(".");
            fflush(stdout);
        }
    }

    // apply additional filters to artifact channels
    if (vflg) {
        printf(" - done\n");
        printf("filtering reference channels;\n");
        fflush(stdout);
    }
    for (r=0; r<R; r++) {
        mr = ArtSpec[r].Index;
        if (vflg) {
            printf("'%s'\n", ArtSpec[r].Name);
            fflush(stdout);
        }
        for (t=0; t<T; t++)
            In[t] = Data[mr][t];
        FFTfilter(In, Out, &ArtSpec[r].Filter, T);
        for (t=0; t<T; t++)
            Data[mr][t] = Out[t];
    }

    // compute sum of squares of each reference channel
    for (r=0; r<R; r++) {
        mr = ArtSpec[r].Index;
        for (t=0, tmp=0.; t<T; t++)
            tmp += Data[mr][t] * Data[mr][t];
        SosA[r] = tmp;      // don't normalize sums of squares
    }

    // compute artifact channel regression matrix & perform LU decomposition
    for (i=0; i<R; i++) {
        mi = ArtSpec[i].Index;
        for (j=0; j<R; j++) {
            mj = ArtSpec[j].Index;
            for (t=0, tmp=0.; t<T; t++)
                tmp += Data[mi][t] * Data[mj][t];
            gsl_matrix_set(A, i, j, tmp / SosA[j]);
        }
    }
    gsl_linalg_LU_decomp(A, perm, &signum);

    // clean each primary channel & write temporary file
    if (vflg) {
        printf("decorrelating artifacts:\n");
        fflush(stdout);
    }
    for (p=0; p<P; p++) {
        mp = Header.PsIndex[p];
        if (vflg) {
            printf("cleaning '%s'\n", Channel[mp].ChannelName);
            fflush(stdout);
        }

        // compute sum of squares
        for (t=0, tmp=0.; t<T; t++)
            tmp += Data[mp][t] * Data[mp][t];
        SosP[p] = tmp;          // don't normalize sums of squares

        // compute primary - artifact channel regression matrix
        for (r=0; r<R; r++) {
            mr = ArtSpec[r].Index;
            for (t=0, tmp=0.; t<T; t++)
                tmp += Data[mp][t] * Data[mr][t];
            gsl_vector_set(B, r, tmp / SosA[r]);
        }

        // solve gain coefficient & regress out references from primary sensor
        gsl_linalg_LU_solve(A, perm, B, C);
#ifdef DEBUG
        for (r=0; r<R; r++) {
            printf("\tgain[%d]=%f\n", r, gain);
            fflush(stdout);
        }
#endif
        for (r=0; r<R; r++) {
            mr = ArtSpec[r].Index;
            for (t=0; t<T; t++)
                Data[mp][t] -= gsl_vector_get(C, r) * Data[mr][t];
        }

       // convert segment to 32-bit byte-reversed integers
        for (t=0; t<T; t++) {
            tmp32 = (int32_t)((double)Data[mp][t] / Channel[mp].Scale);
            xint[t] = bswap_32(tmp32);
        }

        // write file
        sprintf(fpath, "%s.ds/%s.tmp", SetName, Channel[mp].ChannelName);
        if ((fdtmp = open(fpath, O_WRONLY | O_CREAT | O_TRUNC, 0666)) == -1)
            Cleanup("can't open '.tmp' file for write");

        // write unscaled 32-bit integer output to disk
        if (write(fdtmp, (void *)xint, NumBytes) != NumBytes)
            Cleanup("can't write data to '.tmp' file");

        close(fdtmp);
    }

    // open '.meg4.NEW' file
    if (vflg) {
        printf("copying temporary files to .meg4.NEW\n");
        fflush(stdout);
    }
    sprintf(fpath, "%s.ds/%s.meg4", SetName, SetName);
    if ((fdin = open(fpath, O_RDONLY)) == -1)
        Cleanup("can't open '.meg4' file for read");
    sprintf(fpath, "%s.ds/%s.meg4.NEW", SetName, SetName);
    if ((fdout = open(fpath, O_WRONLY | O_CREAT | O_TRUNC, 0666)) == -1)
        Cleanup("can't open '.meg4' file for write");
    if (write(fdout, MEG4Hdr, 8) != 8)
        Cleanup("can't write header to '.meg4' file");

    // write all channels in the same order as the original dataset
    for (m=0; m<M; m++) {

        // if this is a primary (cleaned) sensor
        if (Channel[m].ChannelType == TYPE_MEG) {

            // open & read temporary file
            sprintf(fpath, "%s.ds/%s.tmp", SetName, Channel[m].ChannelName);
            if ((fdtmp = open(fpath, O_RDONLY)) == -1)
                Cleanup("can't open '.tmp' file for write");
            if (read(fdtmp, (void *)xint, NumBytes) != NumBytes)
                Cleanup("can't read data from '.tmp' file");
            close(fdtmp);

            // write cleaned primary sensor to .meg4 file
            if (write(fdout, (char *)xint, NumBytes) != NumBytes)
                Cleanup("can't write data to '.meg4' file");  

        } else {

            // seek & read input
            offset = (off_t)(NumBytes * m) + (off_t)8;
            if (lseek(fdin, offset, SEEK_SET) != offset)
                Cleanup("'lseek()' to channel failed");
            if (read(fdin, (void *)xint, NumBytes) != NumBytes)
                Cleanup("can't read data from 'meg4' file");

            // copy to output
            if (write(fdout, (void *)xint, NumBytes) != NumBytes)
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
    if (dflg) {     // delete backup
        if (vflg) {
            printf(" - done\n");
            printf("deleting backup (original) '.meg4' file");
            fflush(stdout);
        }
        sprintf(OldName, "%s.ds/%s.meg4.bak", SetName, SetName);
        if (unlink(OldName) == -1)
            Cleanup("can't delete backup data file");
    }
	if(vflg == TRUE) {
		printf(" - done\n");
		printf("ArtRedux: done\n");
		fflush(stdout);
	}
	exit(0);
}
