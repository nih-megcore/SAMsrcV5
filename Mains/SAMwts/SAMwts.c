// SAMwts -- linear synthetic aperture magnetometry coefficient solution.
//  Linear solution for the Cartesian moment vector, using multiple
//  local spheres. Orientation vector is determined from a separate
//  covariance matrix.
//
// Features of this version include: reading a cortical atlas file, reading
//  a coordinate transform, & generating target weights using atlas nearest
//  neighbor coordinates. Note that the presence of an atlas over-rides the
//  Talairach transformed weights.
//
// Additional changes require that the MRI directory pathname is provided
//  so that hull.shape & atlas mesh files can be found. This enforces the
//  directory hierarchy where all MRIs are under a single directory and
//  identified by a prefix such as the 8-character hashcode.
//
//  flags:
//      -r <run name>       -- MEG dataset name
//      -d <data file name> -- MEG data file name (4D, only)
//      -m <parameters>     -- parameter file name
//      -B                  -- output primary sensor forward solutions for target list (along with weights)
//      -Z                  -- output SAM weights for normalized (SNR) source projection
//      -v                  -- Verbose mode
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <float.h>
#include <string.h>
#include <strings.h>
#include <fcntl.h>
#include <model.h>
#include <voxel.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <version.h>
#include "samutil.h"
#include "sam_parse.h"

//#define MINOR_REV       1           // changed Wd from double ** to gsl_matrix * -- plus rouutines to handle it
//#define MINOR_REV       2           // when no ROI and no targets are designated in parameters, targets are in <dataset>/SAM/Max
#define MINOR_REV       3           // new command line parser & multithreading

#define ETA             1.0e-12     // a small fraction of a step

typedef struct {
    gsl_matrix  *Cov;               // Cov - MxM covariance matrix
    gsl_matrix  *Cinv;              // Cinv - inverse covariance matrix
    gsl_matrix  *Wd;                // Wd - VxM primary sensor weight vector
    double      *CN;                // CN[V] - condition number for each vertex
    double      Noise;              // noise power/variance -- including optional regularization
    char        Name[64];           // short name
    char        CovFileName[256];   // covariance pathname
    char        WgtFileName[256];   // weight pathname
    char        CNFileName[256];    // where to store condition numbers
    char        NiiName[64];        // nifti weight file name
    COV_HDR     WgtCovHdr;          // weight covariance header
    SAM_HDR     SAMHeader;          // weight file header
    int         Used;               // flag for computing weights
} COV_WGTS;

static int vflg = FALSE;            // verbose mode flag (default: silent)
                                    // global for fixBounds()
static void fixBounds(PARMINFO *p);

extern int CoeffFlag;               // initialization for Nolte solution
extern int ORDER;                   // spherical harmonic expansion order

static void saveCN(char *name, double *cn, int n)
{
    int i;
    FILE *f;

    f = fileopen(name, "w");
    for (i = 0; i < n; i++) {
        fprintf(f, "%g\n", cn[i]);
    }
    fclose(f);
}

int main(
    int             argc,
    char            **argv)
{
    COV_WGTS        *Stats;             // analysis parameters for covariance files
    PARMINFO        Params;             // analysis parameters
    PARM            *p;                 // pointer into active parameter list
    HULL            Hull;               // cortical hull model -- vertices & normals
    HeaderInfo      Header;             // MEG data header
    ChannelInfo     *Channel;           // MEG channel info
    ChannelInfo     *OrigChannel;       // original MEG channel info (before transform)
    EpochInfo       *Epoch;             // MEG epoch info
    VOXELINFO       *Voxel;             // Voxel[V] -- structure for voxel info (including targets)
    VOXELINFO       *Atlas;             // Atlas[A] -- structure for atlas voxel info
    VOXELINFO       *Target;            // Target[V] -- structure for reading targets
    VOXELINFO       VoxDelta;           // interpolated voxel
    GIIATLAS        *atlas;             // gifti atlas info
    COEFFS          *coeff;             // structure to store detection points and coefficients
    struct dirent   *dp;                // directory entry structure
    DIR             *dirp;              // directory pointer
    gsl_vector      *Wgt;               // Wgt - Mx1 beamformer coefficients
    gsl_vector      *Bd;                // Bd - Mx1 primary sensor lead field
    gsl_matrix      *Bdarray;           // VxM primary sensor lead field for each voxel
    gsl_vector      *IntW;              // IntW - Mx1 integrated beamformer coefficients
    gsl_vector_view diag;               // view of the diagonal of a matrix
    double          HPFreq = 0.;        // highpass frequency
    double          LPFreq = 0.;        // lowpass frequency
    double          Src;                // voxel source power
    double          Nse;                // voxel noise power
    double          Mu;                 // regularization
    double          percent;            // floating percentage
    double          span;               // singular value span for solution
    double          min;                // minimum distance
    double          r2;                 // distance^2
    int             ndim = 0;           // # of dimensions for pinv()
    register double tmp;                // what else?
    int             *ChanIndex;         // ChanIndex[M] -- index of channels used for covariance
    int             A;                  // number of atlas voxels
    int             M;                  // number of MEG or reference channels used
    int             R;                  // number of physical reference SQUID channels
    int             S;                  // number of physical SQUID channels
    int             T;                  // number of target voxels
    int             V;                  // number of voxels
    int             c;                  // virtual channel
    int             m;                  // channel index
    int             n;                  // covariance type index
    int             v;                  // vector index
    int             i;                  // matrix index
    int             j;                  // matrix index
    int             ix;                 // x-interpolation index
    int             iy;                 // y-interpolation index
    int             iz;                 // z-interpolation index
    int             index = 0;          // index for nearest neighbor search
    int             pct;                // integer percentage
    int             N;                  // pct counter
    int             a;                  // atlas index
    int             t;                  // target index
    int             NumCov;             // how many covariance matrices we make
    int             *vidx;              // array for shuffling voxel indices
    static int      Order = 3;          // coil integration order -- default 3rd-order
    static int      aflg = FALSE;       // atlas flag
    static int      gflg = FALSE;       // use gifti
    static int      bflg = FALSE;       // forward solution flag
    static int      eflg = FALSE;       // command-line error flag
    static int      iflg = FALSE;       // interpolation flag
    static int      mflg = FALSE;       // parameter file flag
    static int      pflg = FALSE;       // use pseudo-inverse
    static int      rflg = FALSE;       // dataset name flag
    static int      tflg = FALSE;       // target file flag
    static int      wflg = FALSE;       // make Z-weights
    time_t          secs;               // random number seed
    char            *s, *s0;            // string pointers
    char            fpath[256];         // general path name
    char            *DSName = NULL;     // MEG dataset name
    char            DSpath[256];        // MEG dataset path
    char            SAMpath[256];       // SAM subdirectory path
    char            HullPath[256];      // hull.shape path
    char            AtlasPath[256];     // atlas path
    char            TargetPath[256];    // target path to where targets are found
    char            TargetName[256];    // target file name
    char            TlrcPath[256];      // ortho+tlrc path
    char            CovDirName[256];    // covariance subdirectory name
    char            WtsDirName[256];    // weights subdirectory name
    char            *CovName = NULL;    // name used in input file names
    char            *WtsName = NULL;    // name used in weights file names
    char            CurrentDir[256];    // current directory
    char            ByYourCommand[256]; // string for system call
    char            FwdID[] =           // SAM fwd solution ID
        "SAMFWDSL";
    FILE            *fp;
    FILE            *np;
    FILE            *bp = NULL;         // forward solution pointer

    // system-dependent variables
#if BTI
    static int      dflg;               // pdf file name flag
    char            *PDFName;           // name of data file
#else
    HEAD_MODEL      Model;              // model structure for .hdm file
    int             found;              // flag for finding sensor name matches
    int             k;                  // k-index
#endif

    // Register the parameters used by this program.

    reg_usage(argv[0], MINOR_REV, __DATE__);
    reg_std_parm();

    reg_parm("Marker");
    reg_parm("SegFile");

    reg_parm("CovBand");
    reg_parm("NoiseBand");
    reg_parm("Mu");

    reg_parm("XBounds");
    reg_parm("YBounds");
    reg_parm("ZBounds");
    reg_parm("ImageStep");

    reg_parm("MRIDirectory");
    reg_parm("MRIPattern");
    reg_parm("PrefixLength");
    reg_parm("HullName");

    reg_parm("Normalize");
    set_parm_help("Normalize", "normalize the SAM weights by sqrt(noise)");
    reg_parm("Field");

    reg_parm("Pinv");

    reg_parm("Model");
    reg_parm("Order");
    reg_parm("doMultiSphere");  // since Model isn't allowed on the command line
    reg_parm("doNolte");        // allow --domulti or --donolte

    reg_parm("AtlasName");
    reg_parm("TargetName");
    reg_parm("ImageFormat");
    reg_parm("Extent");

    reg_parm("Transform");

    reg_parm("CovName");
    reg_parm("WtsName");
    reg_parm("OutName");
    set_parm_help("OutName", "use this name for --CovName and --WtsName");

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

    p = get_parm("Normalize");      // -n
    wflg = p->set;
    p = get_parm("Field");          // -B
    bflg = p->set;
    p = get_parm("Pinv");           // -p
    pflg = p->set;
    if (pflg) {
        ndim = *(int *)p->ptr;
    }

#if BTI
    if (eflg || !rflg || !dflg || !mflg) {
#else
    if (eflg || !rflg || !mflg) {
#endif
        msg("dataset (-r) and parameter file (-m) are required\n");
        do_help();  // doesn't return
    }

    CovName = Params.ParmName;
    p = get_parm("OutName");
    if (p->set) {
        CovName = (char *)p->ptr;
    }
    p = get_parm("CovName");
    if (p->set) {
        CovName = (char *)p->ptr;
    }
    WtsName = CovName;
    p = get_parm("WtsName");
    if (p->set) {
        WtsName = (char *)p->ptr;
    }

    // check for model override
    j = 0;
    p = get_parm("doMultiSphere");
    if (p->set) {
        Params.Model = MSPHERE;
        j++;
    }
    p = get_parm("doNolte");
    if (p->set) {
        Params.Model = NOLTE;
        j++;
    }
    if (j == 2) {
        Cleanup("can't specify both MultiSphere and Nolte");
    }

    // set the random number seed
    time(&secs);
    srandom((unsigned int)secs);

    // announce
    if (vflg) {
        printf("S A M   C O E F F I C I E N T   G E N E R A T O R\t%s\n", PRG_REV);
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
    GetDsInfo(DSName, &Header, &Channel, &Epoch, NULL, TRUE);
    sprintf(DSpath, "%s/%s.ds", Header.DsPath, Header.SetName);
#endif
    GetFilePath(DSDIR, SAMpath, sizeof(SAMpath), &Params, "SAM", 0);

    // set channel & epoch count constants
    M = Header.NumPri;
    R = Header.NumRef;
    S = Header.NumPri + Header.NumRef;

    ORDER = Params.Order;

    // test for MRI directory pathname for use with hull name, cortical atlas, & Talairach transform
    p = get_parm("MRIDirectory");
    if (!p->set) {
        if (Params.Model != SSPHERE)                // MRI not req'd for single sphere model
            fatalerr("MRIDirectory root req'd for use with Nolte or multi-sphere model");
    } else {
        p = get_parm("AtlasName");
        if (p->set) {
            GetMRIPath(AtlasPath, sizeof(AtlasPath), &Params, (char *)p->ptr, TRUE);
            Params.ImageFormat = ORIG;              // @@@ force ORIG if atlas is present
            aflg = TRUE;

            /* Check for GIFTI atlas. */

            i = strlen(AtlasPath);
            if (i > 4 && strcmp(&AtlasPath[i - 4], ".gii") == 0) {
                gflg = TRUE;
            }
        }
        p = get_parm("TargetName");
        if (p->set) {
            GetMRIPath(TargetName, sizeof(TargetName), &Params, (char *)p->ptr, TRUE);
            tflg = TRUE;
        } else {                                    // no target file
            if (!aflg && Params.SAMStep == -999.) { // if no ROI, then target is in <datasetname>/SAM/Max

                // read directory entries for target
                if (vflg) {
                    printf(" - done\n");
                    printf("searching for target file");
                    fflush(stdout);
                }
                GetFilePath("%d/SAM/Max", TargetPath, sizeof(TargetPath), &Params, "", 0);
                dirp = opendir(TargetPath);
                if (dirp == NULL) {
                    fatalerr("No ROI and target directory %s not found", TargetPath);
                }
                while ((dp = readdir(dirp)) != NULL) {
                    s0 = dp->d_name;
                    i = strlen(s0);
                    if (!strcmp(s0 + i - 4, ".max")) {          // Yes, this is a .max file
                        s = copy_string(s0);
                        s[i - 4] = '\0';                        // strip off .max
                        Params.TargetName = s;
                        GetFilePath("%d/SAM/Max/%s", TargetName, sizeof(TargetName), &Params, s0, 0);
                        if (vflg) {
                            printf(" - found");
                            fflush(stdout);
                        }
                        tflg = TRUE;
                        break;
                    }
                }
                closedir(dirp);
                if (vflg && !tflg) {
                    printf(" - done\n");
                    fatalerr("target file not found");
                }
            }
        }
        if (Params.Model == NOLTE || Params.Model == MSPHERE) {
            s = "hull.shape";
            p = get_parm("HullName");
            if (p->set) {
                s = (char *)p->ptr;
            }
            GetMRIPath(HullPath, sizeof(HullPath), &Params, s, TRUE);
            // read cortical hull model
            if (vflg) {
                printf(" - done\n");
                printf("reading cortical hull model");
                fflush(stdout);
            }
            GetHull(HullPath, &Hull);               // read hull model
        }
        if (Params.ImageFormat == TLRC) {           // make sure that we can access brain+tlrc for @auto_tlrc step
            GetMRIPath(TlrcPath, sizeof(TlrcPath), &Params, "brain+tlrc.HEAD", TRUE);
            GetMRIPath(TlrcPath, sizeof(TlrcPath), &Params, "brain+tlrc", FALSE);
            if (Params.ImageRes < (1000. * Params.SAMStep) || Params.ImageRes > 20.)
                Cleanup("Talairach resolution must be >= ortho image step-size & < 20 mm");
        }
    }

    // if transform is present, move all MEG primary and reference sensors
    p = get_parm("Transform");
    if (p->set) {
        if (vflg) {
            printf(" - done\n");
            printf("reading Transform");
            fflush(stdout);
        }
        GetFilePath("%d/SAM/%s.xfm", fpath, sizeof(fpath), &Params, (char *)p->ptr, 0);
        OrigChannel = new_array(ChannelInfo, Header.NumChannels);
        memcpy(OrigChannel, Channel, sizeof(ChannelInfo) * Header.NumChannels);
        MoveFrame(&Header, OrigChannel, Channel, fpath);
    }

    // determine weight interpolation step
    if (Params.Extent > 0.) {
        if (Params.Extent <= 0.5 * Params.SAMStep) {
            if (Params.Model != NOLTE)
                iflg = TRUE;
                iflg = FALSE;
        } else
            Cleanup("'Extent' step size must be <= 'SAMstep'");
    }

    if (vflg) {
        printf(" - done\n");
    }

    // Make sure that the ROI is an integral number of voxels.

    if (get_parm("XBounds")->set && get_parm("YBounds")->set &&
        get_parm("ZBounds")->set && get_parm("ImageStep")->set) {
        fixBounds(&Params);
    }

    // define forward solution model
    if (Params.Model == -1)
        Cleanup("'Model' specification missing from parameter file");
    if (vflg) {
        printf("modeling forward solution using ");
    }
    switch (Params.Model) {

        // single sphere origin -- note: dither the sphere origins by a small random offset so that the eigensystem solution for moment vector succeeds!
        case SSPHERE:
            if (vflg) {
                printf("single local sphere origin: { %6.3f, %6.3f, %6.3f } (cm)\n", 100. * Params.Sp[X_], 100. * Params.Sp[Y_], 100. * Params.Sp[Z_]);
            }
            for (i=0; i<M; i++) {
                j = Header.PsIndex[i];
                for (v=X_; v<=Z_; v++)
                    Channel[j].Geom.MEGSensor.LocalSphere[v] = Params.Sp[v] + RandomGauss(0., 1.0e-6);
            }
            for (i=0; i<R; i++) {
                j = Header.RsIndex[i];
                for (v=X_; v<=Z_; v++)
                    Channel[j].Geom.MEGSensor.LocalSphere[v] = Params.Sp[v] + RandomGauss(0., 1.0e-6);
            }
            Hull.nv = 1;
            if ((Hull.vertex = (FIELD *)malloc(2 * sizeof(FIELD))) == NULL)
                Cleanup("'Hull.vertex[1]' allocation failed");
            for (v=X_; v<=Z_; v++)
                Hull.Vo[v] = Params.Sp[v];
            break;

        // multiple local spheres
        case MSPHERE:
            if (vflg) {
                printf("multiple local spheres\n");
            }
#if BTI
            // read hs_file into 'Model'
            HStoMS(DSpath, &Header, Channel);
#else
            // read .hdm file into 'Model'
            sprintf(fpath, "%s/default.hdm", DSpath);
            GetHDM(fpath, &Model);

            // test if the model matches the number of sensors
            if (Model.NumSensors != S)
                Cleanup("number of sensors in '.hdm' file does not match that of the dataset");

            // set multiple sphere origins
            for (i=0; i<S; i++) {
                k = (int)strlen(Model.MultiSphere[i].ChanName) - 1;
                for (j=0, found=FALSE; j<Header.NumChannels; j++)
                    if (!strncmp(Model.MultiSphere[i].ChanName, Channel[j].ChannelName, k)) {
                        for (v=X_; v<=Z_; v++)
                            Channel[j].Geom.MEGSensor.LocalSphere[v] = Model.MultiSphere[i].Origin[v];
                        found = TRUE;
                    }
                if (!found)
                    Cleanup("could not find matching channel name in '.HDM' file");
            }
#endif
            break;

        // Nolte realistic head model
        case NOLTE:
            if (vflg) {
                printf("%dth-order realistic head model based on\n\t'%s'\n", ORDER, HullPath);
            }
            CoeffFlag = TRUE;
            break;

        default:
            break;
    }

    // check for atlas
    if (aflg) {
        if (vflg) {
            printf("reading '%s' atlas file instead of using ROI", Params.Atlas);
            fflush(stdout);
        }
        if (gflg) {
            atlas = GetGiiAtlas(AtlasPath);
            Atlas = atlas->Voxel;
            A = atlas->nrh + atlas->nlh;
            if (vflg) {
                printf("\n%d atlas points\n", A);
            }
        } else {
            GetVoxelsFS(AtlasPath, &Atlas, &A, Params.SAMStart, Params.SAMEnd, &Params.SAMStep);
            if (vflg) {
                printf("\n%d voxels\n", A);
                printf("XBounds: %g to %g cm\n", 100.*Params.SAMStart[X_], 100.*Params.SAMEnd[X_]);
                printf("YBounds: %g to %g cm\n", 100.*Params.SAMStart[Y_], 100.*Params.SAMEnd[Y_]);
                printf("ZBounds: %g to %g cm\n", 100.*Params.SAMStart[Z_], 100.*Params.SAMEnd[Z_]);
                printf("ImageStep %g cm\n", 100.*Params.SAMStep);
            }
        }
        if (!tflg) {
            Voxel = Atlas;
            V = A;
        }
    }

    // allocate voxels with boundary flags
    if (!tflg) {    // no target files specified in parameter file
        if (vflg) {
            printf("generating ROI voxel array");
            fflush(stdout);
        }
        if (aflg) {
            printf(" using '%s' atlas", Params.Atlas);
            fflush(stdout);
        } else {                                    // if there NO Target file and NO Atlas found, 'Voxel' is allocated!
            if (Params.SAMStart[X_] == -999. || Params.SAMStart[Y_] == -999. || Params.SAMStart[Z_] == -999.)        // no need to test SAMend...
                Cleanup("ROI 'XBounds/YBounds/ZBounds' specifications required");
            if (Params.SAMStep == -999.)
                Cleanup("ROI 'ImageStep' specification required");
#if BTI
            GetVoxelsHS(&Voxel, &V, DSpath, Params.SAMStart, Params.SAMEnd, Params.SAMStep);
#else
            GetVoxels(&Voxel, &V, &Hull, Params.SAMStart, Params.SAMEnd, Params.SAMStep);
#endif
            if (vflg) {
                printf(" of %d voxels", V);
                fflush(stdout);
            }
        }
    } else {                                        // a Target file is present!
        if (vflg) {
            printf("reading target list file");
            fflush(stdout);
        }

        // atlas present -- find nearest neighbors, & set orientations
        //  NOTE:   'Target' is an array of target coordinates obtained by 'GetTargets()' where T is the number of target coordinates
        //          'Atlas' is an array of atlas coordinates obtained by 'GetVoxelsFS()' where A is the number of atlas coordinates
        if (aflg) {
            GetTargets(TargetName, &Target, &T);
            if (vflg) {
                printf(" - read %d targets", T);
                fflush(stdout);
            }

            // for each target voxel...
            for (t=0; t<T; t++) {

                // scan atlas for nearest neighbor
                for (a=0, min=HUGE; a<A; a++) {
                    for (v=X_, r2=0.; v<=Z_; v++) {
                        tmp = Atlas[a].p[v] - Target[t].p[v];
                        r2 += tmp * tmp;
                    }
                    if (r2 < min) {
                        min = r2;
                        index = a;
                    }
                }
                if (sqrt(min) > 0.02)
                    msg("\nWarning: target %d was more than 20 mm from nearest atlas voxel\n", t);
                for (v=X_; v<=Z_; v++) {         // set the target position & moment vectors to nearest neighbor in atlas
                    Target[t].p[v] = Atlas[index].p[v];
                    Target[t].v[v] = Atlas[index].v[v];
                }
            }
            Voxel = Target;     // set Voxel pointer to Target
            V = T;
        } else {
            GetTargets(TargetName, &Voxel, &V);
            if (vflg) {
                printf(" - read %d targets", V);
                fflush(stdout);
            }
        }
    }

    // always use CovBand frequencies
    if (Params.CovHP >= 0. && Params.CovLP >= 0.) {
        HPFreq = Params.CovHP;
        LPFreq = Params.CovLP;
    } else {
        Cleanup("CovBand not specified!");
    }

    // initialize flags
    NumCov = Params.NumMark + MARKER1_;
    Stats = new_array(COV_WGTS, NumCov);
    for (n=0; n<NumCov; n++)
        Stats[n].Used = FALSE;

    // generate covariance subdirectory full pathname string by concatenating SAMpath, CovName, & bandpass
    sprintf(fpath, "%s,%-d-%-dHz", CovName, (int)HPFreq, (int)LPFreq);
    GetFilePath("%d/SAM/%s", CovDirName, sizeof(CovDirName), &Params, fpath, 0);
    sprintf(fpath, "%s,%-d-%-dHz", WtsName, (int)HPFreq, (int)LPFreq);
    GetFilePath("%d/SAM/%s", WtsDirName, sizeof(WtsDirName), &Params, fpath, 0);

    // create Wts subdirectory if different from CovDir
    if (strcmp(WtsName, CovName) != 0) {
        if (mkdir(WtsDirName, 0755) == -1)
            if (errno != EEXIST)
                Cleanup("can't create weights subdirectory");
    }

    // count actual number of covariance matrices
    if (vflg) {
        printf(" - done\n");
        printf("determining covariance & weight names");
        fflush(stdout);
    }

    // populate full pathnames for Global
    sprintf(Stats[GLOBAL_].Name, "Global");
    sprintf(Stats[GLOBAL_].CovFileName, "%s/Global.cov", CovDirName);
    sprintf(Stats[GLOBAL_].WgtFileName, "%s/Global.nii", WtsDirName);
    sprintf(Stats[GLOBAL_].CNFileName, "%s/GlobalCN.dat", WtsDirName);
    sprintf(Stats[GLOBAL_].NiiName, "Global.nii");
    Stats[GLOBAL_].Used = TRUE;

    // populate full pathname for Orient
    sprintf(Stats[ORIENT_].Name, "Orient");
    sprintf(Stats[ORIENT_].CovFileName, "%s/Orient.cov", CovDirName);
    sprintf(Stats[ORIENT_].WgtFileName, "NULL");
    sprintf(Stats[ORIENT_].CNFileName, "%s/OrientCN.dat", WtsDirName);
    sprintf(Stats[ORIENT_].NiiName, "NULL");
    Stats[ORIENT_].Used = TRUE;

    // check for Noise covariance
    if (Params.NoiseHP != -999. && Params.NoiseLP != -999.) {
        sprintf(Stats[NOISE_].Name, "Noise");
        sprintf(Stats[NOISE_].CovFileName, "%s/Noise.cov", CovDirName);
        sprintf(Stats[NOISE_].WgtFileName, "%s/Noise.nii", WtsDirName);
        sprintf(Stats[NOISE_].CNFileName, "%s/NoiseCN.dat", WtsDirName);
        sprintf(Stats[NOISE_].NiiName, "Noise.nii");
        Stats[NOISE_].Used = TRUE;
    }

    // now check markers
    if (Params.NumMark > 0) {

        // include Sum covariance
        sprintf(Stats[SUM_].Name, "Sum");
        sprintf(Stats[SUM_].CovFileName, "%s/Sum.cov", CovDirName);
        sprintf(Stats[SUM_].WgtFileName, "%s/Sum.nii", WtsDirName);
        sprintf(Stats[SUM_].CNFileName, "%s/SumCN.dat", WtsDirName);
        sprintf(Stats[SUM_].NiiName, "Sum.nii");

        for (i=0; i<Params.NumMark; i++) {
            s = Params.Marker[i].MarkName2;
            if (s == NULL) {
                s = Params.Marker[i].MarkName;
            }
            j = MARKER1_ + i;                   // begin at marker1
            sprintf(Stats[j].Name, "%s", s);
            sprintf(Stats[j].CovFileName, "%s/%s.cov", CovDirName, s);
            sprintf(Stats[j].WgtFileName, "%s/%s.nii", WtsDirName, s);
            sprintf(Stats[j].CNFileName, "%s/%sCN.dat", WtsDirName, s);
            sprintf(Stats[j].NiiName, "%s.nii", s);
            Stats[j].Used = TRUE;
            if (Params.Marker[i].Sum) {
                Stats[SUM_].Used = TRUE;
            }
        }
    }

    // allocate memory for arrays
    if (vflg) {
        printf(" - done\n");
        printf("allocating memory");
        fflush(stdout);
    }
    for (n=0; n<NumCov; n++) {
        if (Stats[n].Used) {
            Stats[n].Cov = gsl_matrix_alloc(M, M);
            Stats[n].Cinv = gsl_matrix_alloc(M, M);
            Stats[n].Wd = gsl_matrix_alloc(V, M);
            Stats[n].CN = new_array(double, V);
        }
    }

    // read covariance matrices
    if (vflg) {
        printf(" - done\n");
        printf("reading covariance files:\n");
    }

    for (n=0; n<NumCov; n++)
        if (Stats[n].Used) {

            // read covariance matrix
            if (vflg) {
                printf("\t'%s'", Stats[n].Name);
                fflush(stdout);
            }
            GetCov(Stats[n].CovFileName, &Stats[n].WgtCovHdr, &ChanIndex, Stats[n].Cov);

            // read noise files
            sprintf(fpath, "%s/%s_Noise", CovDirName, Stats[n].Name);
            if ((np = fopen(fpath, "r")) == NULL)
                Cleanup("can't open noise file");
            if (fscanf(np, "%le", &Stats[n].Noise) != 1)
                Cleanup("can't read noise file");
            fclose(np);

            diag = gsl_matrix_diagonal(Stats[n].Cov);

            // if req'd, regularize covariance matrix
            switch (Params.Operator) {
                case ADD_MU:    // add constant noise power to diagonal of covariance matrix
                    tmp = 1.0e-30 * Params.Mu * Params.Mu;      // convert fT/root-Hz to T^2/Hz
                    Mu = Stats[n].WgtCovHdr.BWFreq * tmp;       // multiply by covariance bandwidth to get T^2
                    gsl_vector_add_constant(&diag.vector, Mu);
                    Stats[n].Noise += Mu;
                    break;
                case MPY_MU:    // add proportionally scaled noise power to diagonal of covariance matrix
                    Mu = Stats[n].Noise * Params.Mu;
                    gsl_vector_add_constant(&diag.vector, Mu);
                    Stats[n].Noise += Mu;
                    break;
                default:        // no regularizarion -- do nothing
                    break;
            }

            if (vflg) {
                printf(" - done\n");
            }

            // write new noise files after regularization
            if (Params.Operator != -1) {
                sprintf(fpath, "%s/%s_MuNoise", WtsDirName, Stats[n].Name);
                if ((np = fopen(fpath, "w")) == NULL)
                    Cleanup("can't write noise file");
                fprintf(np, "%e\n", Stats[n].Noise);
                fclose(np);
            }
        }

    if (wflg && vflg) {     // normalize
        printf("normalizing weights by sqrt(regularized noise power)\n");
        for (n=0; n<NumCov; n++) {
            if (Stats[n].Used) {
                printf("normalizing %s weights by %g\n", Stats[n].Name, sqrt(Stats[n].Noise));
            }
        }
    }

    // announce regularization used
    if (vflg) {
        if (Params.Operator != -1) {
            printf("regularizing covariance matrices by ");
            switch (Params.Operator) {
                case ADD_MU:
                    printf("adding %6.3f fT/root-Hz to diagonal\n", Params.Mu);
                    break;
                case MPY_MU:
                    printf("multiplying diagonal by %6.3f\n", Params.Mu);
                    break;
                default:
                    break;
            }
        }
    }

    // Invert covariance matrices.

    for (n=0; n<NumCov; n++) {
        if (Stats[n].Used) {
            if (pflg) {
                pinvN(Stats[n].Cov, Stats[n].Cinv, ndim);
            } else {
                pinv(Stats[n].Cov, Stats[n].Cinv);
            }
        }
    }

    // create shape file (used for checking that subject's head is fully inside helmet)
    if (vflg) {
        printf("saving '%s.shape'", Header.SetName);
        fflush(stdout);
    }
    GetFilePath("%d/%s.shape", fpath, sizeof(fpath), &Params, Header.SetName, 0);
    fp = fileopen(fpath, "w");
    fprintf(fp, "%d\n", M);
    for (m=0; m<M; m++) {
        j = ChanIndex[m];
        fprintf(fp, "%7.3f\t%7.3f\t%7.3f\n",
            100. * Channel[j].Geom.MEGSensor.origin.p[X_], 100. * Channel[j].Geom.MEGSensor.origin.p[Y_], 100. * Channel[j].Geom.MEGSensor.origin.p[Z_]);
    }
    fclose(fp);

    // if req'd, open SAM forward solution files & write header
    if (bflg) {
        if (vflg) {
            printf(" - done\n");
            printf("opening forward solution file for write");
            fflush(stdout);
        }
        n = GLOBAL_;
        Stats[n].SAMHeader.Version = SAM_REV;                   // file version number
        Stats[n].SAMHeader.NumChans = M;                        // number of channels used by SAM
        Stats[n].SAMHeader.NumWeights = V;                      // number of SAM virtual sensors
        Stats[n].SAMHeader.XStart = Params.SAMStart[X_];        // x-start coordinate (m)
        Stats[n].SAMHeader.XEnd = Params.SAMEnd[X_];            // x-end coordinate (m)
        Stats[n].SAMHeader.YStart = Params.SAMStart[Y_];        // y-start coordinate (m)
        Stats[n].SAMHeader.YEnd = Params.SAMEnd[Y_];            // y-end coordinate (m)
        Stats[n].SAMHeader.ZStart = Params.SAMStart[Z_];        // z-start coordinate (m)
        Stats[n].SAMHeader.ZEnd = Params.SAMEnd[Z_];            // z-end coordinate (m)
        if (tflg || gflg) {
                Stats[n].SAMHeader.StepSize = 0.;
        } else {
                Stats[n].SAMHeader.StepSize = Params.SAMStep;   // voxel step size (m)
        }
        sprintf(fpath, "%s/%s.fwd", SAMpath, WtsName);
        bp = fileopen(fpath, "w");
        if (fwrite((void *)FwdID, 8, 1, bp) != 1)
            Cleanup("can't write ID to forward solution file");
        if (fwrite((void *)&Stats[GLOBAL_].SAMHeader, sizeof(SAM_HDR), 1, bp) != 1)
            Cleanup("can't write header to forward solution file");
        if (fwrite((void *)ChanIndex, sizeof(int), M, bp) != M)
            Cleanup("can't write channel index to forward solution file");

        // array to store forward solutions
        Bdarray = gsl_matrix_alloc(V, M);
    }

    // compute integration points and coefficients
    if (vflg) {
        printf(" - done\n");
        printf("computing sensor integration points and coefficients");
        fflush(stdout);
    }
    if (Params.Model == NOLTE) {
        coeff = new(COEFFS);
        ECDIntPnt(&Header, Channel, &Hull, coeff);
    } else {
        for (i = 0; i < Header.NumRef; i++) {
            j = Header.RsIndex[i];
            FwdSensIntPnt(&Channel[j].Geom.MEGSensor, Order);
        }
        for (i = 0; i < Header.NumPri; i++) {
            j = Header.PsIndex[i];
            FwdSensIntPnt(&Channel[j].Geom.MEGSensor, Order);
        }
    }

    // run the voxel loop in random order, so the voxels with .ROI
    // set to false are spread out evenly across threads
    vidx = new_array(int, V);
    for (c = 0; c < V; c++) {
        vidx[c] = c;
    }
    shuffle(vidx, V);

    // compute beamformer coefficients for each marker
    if (vflg) {
        printf(" - done\n");
        printf("computing beamformer coefficients for:\n");
    }
    for (n=0; n<NumCov; n++) {      // loop on each valid marker
        if (Stats[n].Used && n != ORIENT_) {

            if (vflg) {
                printf("\t%s:\n", Stats[n].Name);
                fflush(stdout);
            }

            // loop on voxels
            pct = -1;
            N = 0;
#pragma omp parallel private(Wgt, Bd, IntW)
{
            Wgt = gsl_vector_alloc(M);
            Bd = gsl_vector_alloc(M);
            if (iflg) {
                IntW = gsl_vector_alloc(M);
            }

#pragma omp for private(c, span, Src, Nse, tmp, v, VoxDelta, ix, iy, iz)
            for (j=0; j<V; j++) {
                c = vidx[j];

                // estimate SAM voxel only if ROI flag is true
                if (Voxel[c].ROI) {

                    // set voxel coordinate & solve for optimal beamformer orientation using orientation covariance
                    if (!iflg) {

                        // no interpolation
                        switch (Params.Model) {
                            case SSPHERE:
                            case MSPHERE:
                                DIPSolve(&Voxel[c], Bd, Stats[ORIENT_].Cinv, &Header, Channel, Order, &span);
                                break;
                            case NOLTE:
                                ECDSolve0(&Voxel[c], Bd, Stats[ORIENT_].Cinv, &Header, Channel, &Hull, &span, coeff);
                                break;
                            default:
                                break;
                        }
                        SAMsolve(Stats[n].Cov, Stats[n].Cinv, Bd, Wgt, Stats[n].Noise, &Src, &Nse);
                        if (wflg) {     // noise normalized weights
                            tmp = 1. / sqrt(Stats[n].Noise);
                            gsl_vector_scale(Wgt, tmp);
                        }
                        gsl_matrix_set_row(Stats[n].Wd, c, Wgt);
                        Stats[n].CN[c] = span;
                    } else {    // else interpolate

                        // copy orientation if using Atlas & set flags
                        if (aflg) {
                            for (v=X_; v<=Z_; v++)
                                VoxDelta.v[v] = Voxel[c].v[v];
                            VoxDelta.Solve = FALSE;
                        } else {
                            VoxDelta.Solve = TRUE;
                        }
                        VoxDelta.ROI = TRUE;

                        // zero integrated weight vector
                        gsl_vector_set_zero(IntW);

                        // loop on interpolations
                        for (ix=(-1); ix<=1; ix++) {
                            VoxDelta.p[X_] = Voxel[c].p[X_] + (double)ix * Params.Extent;
                            for (iy=(-1); iy<=1; iy++) {
                                VoxDelta.p[Y_] = Voxel[c].p[Y_] + (double)iy * Params.Extent;
                                for (iz=(-1); iz<=1; iz++) {
                                    VoxDelta.p[Z_] = Voxel[c].p[Z_] + (double)iz * Params.Extent;

                                    switch (Params.Model) {
                                        case SSPHERE:
                                        case MSPHERE:
                                            DIPSolve(&VoxDelta, Bd, Stats[ORIENT_].Cinv, &Header, Channel, Order, &span);
                                            break;
                                        case NOLTE:
                                            ECDSolve0(&VoxDelta, Bd, Stats[ORIENT_].Cinv, &Header, Channel, &Hull, &span, coeff);
                                            break;
                                        default:
                                            break;
                                    }

                                    // solve for weights & sum
                                    SAMsolve(Stats[n].Cov, Stats[n].Cinv, Bd, Wgt, Stats[n].Noise,  &Src, &Nse);
                                    if (wflg) {     // noise normalized weights
                                        tmp = 1. / sqrt(Stats[n].Noise);
                                        gsl_vector_scale(Wgt, tmp);
                                    }

                                    // sum weights
                                    gsl_vector_add(IntW, Wgt);

                                }   // end iz loop
                            }       // end iy loop
                        }           // end ix loop

                        // average weights & set Wd row
                        gsl_vector_scale(IntW, 1./27.);
                        gsl_matrix_set_row(Stats[n].Wd, c, IntW);

                    }   // end iflg

                // if outside bounds
                } else {                // zero forward solution & coefficients for this voxel
                    gsl_vector_set_zero(Bd);
                    gsl_vector_set_zero(Wgt);
                    gsl_matrix_set_row(Stats[n].Wd, c, Wgt);
                    Stats[n].CN[c] = 0.;
                }

                // if req'd, save forward solution vector
                if (bflg) {
                    gsl_matrix_set_row(Bdarray, c, Bd);
                }

                // percent done
#pragma omp critical
                if (vflg) {
                    percent = 100. * (double)N / (double)V;
                    N++;
                    if (rint(percent) > (double)pct) {
                        pct = (int)rint(percent);
                        printf("\r\t\t%d percent done", pct);
                        fflush(stdout);
                    }
                }
            }       // end c-loop

            gsl_vector_free(Wgt);
            gsl_vector_free(Bd);
            if (iflg) {
                gsl_vector_free(IntW);
            }
} // end parallel

            if (vflg) {
                printf("\n");
            }
        }           // end if stats TRUE

        if (bflg) {
            Bd = gsl_vector_alloc(M);
            for (v = 0; v < V; v++) {
                gsl_matrix_get_row(Bd, Bdarray, v);
                gsl_vector_fwrite(bp, Bd);
            }
            gsl_vector_free(Bd);
            gsl_matrix_free(Bdarray);
            fclose(bp);
            bflg = FALSE;
        }
    }

    // if req'd. write legacy weights
    if (tflg) {

        // write legacy weight file to SAM subdirectory
        if (vflg) {
            printf(" - done\n");
            printf("writing legacy weights");
            fflush(stdout);
        }
        for (n=0; n<NumCov; n++)
            if (Stats[n].Used && n != ORIENT_) {
                sprintf(fpath, "%s/%s,%s,%s.wts", SAMpath, WtsName, Params.TargetName, Stats[n].Name);
                PutWts(fpath, Stats[n].Wd, Header, Stats[n].WgtCovHdr);
            }

        // write NIFTI target weights
        if (vflg) {
            printf(" - done\n");
            printf("writing NIFTI weights for targets");
            fflush(stdout);
        }
        for (n=0; n<NumCov; n++)
            if (Stats[n].Used && n != ORIENT_)
                PutNIFTITargets(Stats[n].WgtFileName, Stats[n].Wd, Stats[n].Noise);

    } else {    // no targets -- write weights for ROI

        // write NIFTI weights
        if (vflg) {
            printf(" - done\n");
            printf("writing %cIFTI weights for ROI", gflg ? 'G' : 'N');
            fflush(stdout);
        }
        for (n=0; n<NumCov; n++) {
            if (Stats[n].Used && n != ORIENT_) {
                if (gflg) {
                    PutGIFTIWts(Stats[n].WgtFileName, Stats[n].Wd, Stats[n].Noise);     // write the beamformer coefficients in GIFTI format
                } else {
                    PutNIFTIWts(Stats[n].WgtFileName, &Params, Stats[n].Wd, Stats[n].Noise);    // write the beamformer coefficients to the NIFTI file

                    // save the condition number separately.

                    saveCN(Stats[n].CNFileName, Stats[n].CN, V);
                }
            }
        }
    }

    // if TLRC, execute system call for @auto_tlrc
    if (!tflg && Params.ImageFormat == TLRC) {
        if (vflg) {
            printf(" - done\n");
            printf("computing Talairach transform of weights at %4.1f mm resolution:\n", Params.ImageRes);
        }
        if (getcwd(CurrentDir, sizeof(CurrentDir)) == NULL)
            Cleanup("can't get current working directory");

        /* Create an absolute pathname for -apar if TlrcPath isn't already absolute. */

        s0 = TlrcPath;
        if (TlrcPath[0] != '/') {
            s0 = new_string(strlen(CurrentDir) + strlen(TlrcPath) + 1);
            s = strecpy(s0, CurrentDir);
            *s++ = '/';
            strcpy(s, TlrcPath);
        }

        if (chdir(WtsDirName) == -1)
            Cleanup("can't change to directory %s", WtsDirName);
        for (n=0; n<NumCov; n++) {
            if (Stats[n].Used && n != ORIENT_) {    // note that the verbose output of @auto_tlrc is now eliminated
                if (vflg) {
                    printf("\t%s", Stats[n].Name);
                    fflush(stdout);
                }

                /* @auto_tlrc -overwrite doesn't always work; explicity remove the output name. */

                s = new_string(strlen(Stats[n].Name) + 7); // "_at.nii"
                sprintf(s, "%s_at.nii", Stats[n].Name);
                unlink(s);
                free(s);

                sprintf(ByYourCommand, "@auto_tlrc -rmode NN -dxyz %4.1f -apar %s -input ./%s -overwrite > /dev/null 2>&1", Params.ImageRes, TlrcPath, Stats[n].NiiName);
                if (system(ByYourCommand) != 0)
                    fatalerr("command failed: %s\n", ByYourCommand);
                if (vflg) {
                    printf(" - done\n");
                }
            }
        }
        if (chdir(CurrentDir) == -1)
            Cleanup("can't change back to current directory");
    } else {
        if (vflg) {
            printf(" - done\n");
        }
    }

    log_params(SAMpath);

    // all done!
    if (vflg) {
        printf("'%s' done\n", Progname);
    }
    exit(0);
}

/* Convert ROI bounds so that there are always an integral number of steps
and the computed ROI encloses the given one. This way the voxel boundaries
line up exactly with the coordinate axes. (The origin is at the corner of
a voxel.) lowhigh is -1 or 1 for low or high bounds, resp. */

static double fixBound(double b, double step, int lowhigh)
{
    double t, r;

    if (lowhigh < 0) {
        /* lower bound */
        t = floor(b / step);
        r = (t + .5) * step;
        if (r > b) {
            r = (t - .5) * step;
        }
    } else {
        /* higher bound */
        t = ceil(b / step);
        r = (t - .5) * step;
        if (r < b) {
            r = (t + .5) * step;
        }
    }
    return r;
}

static void fixBounds(PARMINFO *p)
{
    int i, j;
    double tmp;

    j = 0;
    for (i = 0; i < 3; i++) {
        tmp = fixBound(p->SAMStart[i], p->SAMStep, -1);
        if (!j && fabs(tmp - p->SAMStart[i]) > ETA) {
            j = 1;
        }
        p->SAMStart[i] = tmp;
        tmp = fixBound(p->SAMEnd[i], p->SAMStep, 1);
        if (!j && fabs(tmp - p->SAMEnd[i]) > ETA) {
            j = 1;
        }
        p->SAMEnd[i] = tmp;
    }

    if (vflg) {
        if (j) {
            printf("\nROI bounds adjusted to accommodate step size of %g\n\n", p->SAMStep * 100.);
            printf("XBounds %g %g\n", p->SAMStart[0] * 100., p->SAMEnd[0] * 100.);
            printf("YBounds %g %g\n", p->SAMStart[1] * 100., p->SAMEnd[1] * 100.);
            printf("ZBounds %g %g\n\n", p->SAMStart[2] * 100., p->SAMEnd[2] * 100.);
        }
    }
}
