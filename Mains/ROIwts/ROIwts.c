// ROIwts -- compute beamformer weights for extended source defined by ROI list
//  It is assumed that the subject's MRI has been segmented using FreeSurfer.
//  Next, FSnormal.py is run to generate a cortical mesh for the pial, white,
//  & smooth wight matter surfaces. The mesh includes the unit normal vector
//  at each vertex. This version uses the realistic head model (Nolte) forward
//  solution! The cortical normal vector is used for the dipole source vector
//  instead of solving the linear eigensystem equations for the vector giving
//  the largest S/N. It's important to run SAMcoreg prior to this!
//
// In this version, if radial 'Extent' is given in the parameter file, the
//  solution is restricted to vertices within the radial extent from the
//  ROI centroid. This is the fastest solution. If 'Extent' is not given,
//  then the beamformer coefficients for the entire ROI are estimated. This
//  is very slow!
//
//  flags:
//  -r <run name>       -- MEG Run name
//  -m <parameters>     -- parameter file
//  -s                  -- solve unconstrained moment vector
//  -v                  -- Verbose mode
//
//  Author: Stephen E. Robinson
//      MEG Core Facility
//      NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <gsl/gsl_matrix.h>
#include <model.h>
#include <samlib.h>
#include <siglib.h>
#include <lfulib.h>
#include <voxel.h>
#include <DataFiles.h>
#include <SAMfiles.h>
#include <sam_parse.h>
#include <sam_param.h>
#include <nifti1.h>
#include <version.h>

#define MINOR_REV   1                   // new command line parser
#define MIN_VERTEX  10                  // minimum number of vertices for solution

// global variable
extern int  CoeffFlag;                  // flag for reevaluating Nolte coefficients
extern int  ORDER;                      // spherical harmonic expansion order

int main(
    int         argc,
    char        **argv)
{
    HeaderInfo      Header;             // MEG data header
    HULL            Hull;               // cortical hull model -- vertices & normals
    PARMINFO        Params;             // analysis parameters
    PARM            *p;                 // pointer into active parameter list
    VOXELINFO       *Vertex;            // Vertex[N] -- ROI atlas vertex structure list
    ChannelInfo     *OrigChannel;       // untransformed MEG channel info
    ChannelInfo     *Channel;           // cloned MEG channel info for transformation
    EpochInfo       *Epoch;             // MEG epoch info
    COV_HDR         CovHeader;          // covariance file header
    gsl_matrix      *C;                 // C - MxM covariance matrix
    gsl_matrix      *Co;                // Co - MxM orientation covariance matrix
    gsl_matrix      *W;                 // W - VxM SAM weights
    gsl_vector      *Wv;                // Wv - Mx1 SAM weight vector
    gsl_vector      *A;                 // A - Mx1 temporary vector
    double          Vcentroid[3];       // ROI centroid
    double          s2;                 // source power
    double          n2;                 // noise power
    int             *ChanIndex;         // ChanIndex[M] -- index of channels used for covariance
    int             rank;               // number of bases
    int             ROInum;             // ROI index number
    int             ROIhemi;            // ROI hemisphere -- 1 == L, -1 == R
    int             Nused;              // number of ROI vertices used for solution
    int             M;                  // number of primary channels
    int             N;                  // number of vertices in ROI
    int             V;                  // number of targets
    int             c;                  // input character
//  int             i;                  // i-index
//  int             j;                  // j-index
    int             n;                  // vertex index
    int             f;                  // index of 1st primary sensor (for recorded bandpass)
    int             v;                  // vector index
    int             eflg = FALSE;       // command-line error flag
    int             rflg = FALSE;       // dataset flag
    int             mflg = FALSE;       // parameter flag
    int             vflg = FALSE;       // verbose mode flag (default: silent)
    int             uflg;               // unconstrained solution flag
    int             xflg;               // extent flag
    int             tflg;               // transform flag
    char            *s;                 // string pointer
    char            *DSName = NULL;     // MEG dataset name
    char            SAMpath[256];       // data set path
    char            NormPath[256];      // path to .norm atlas
    char            fpath[256];         // general path name
    char            *ParmName = NULL;   // parameter file name
    char            Line[128];          // space for one line of ROI names
    unsigned char   **Bad;              // Bad[E][T] -- flags for bad samples
    FILE            *fp;                // file pointer
    FILE            *fl;                // log file pointer

    // Register the parameters used by this program.

    reg_usage(argv[0], MINOR_REV, __DATE__);
    reg_std_parm();

    reg_parm("CovBand");
    reg_parm("OrientBand");
    reg_parm("MRIDirectory");
    reg_parm("MRIPattern");
    reg_parm("PrefixLength");
    reg_parm("HullName");
    reg_parm("OutName");
    reg_parm("ROIList");
    reg_parm("Surface");
    reg_parm("Transform");
    reg_parm("Model");
    set_parm_arg_help("Model", "Nolte [order]");
    set_parm_help("Model", "forward solution type (must be Nolte)");
    reg_parm("Unconstrained");
    reg_parm("Extent");

    // Parse the arguments, environment variables, and parameter files.

    do_parse_args(argc, argv);

    // Extract specified parameters.

    eflg = get_params(&Params);

    p = get_parm("DataSet");
    if (p->set) {
        rflg = TRUE;
        DSName = copy_string(Params.DataSetName);
    }

    p = get_parm("verbose");
    vflg = p->set;
    p = get_parm("param");
    mflg = p->set;

    if (eflg || !rflg || !mflg) {
        msg("dataset (-r) and parameter file (-m) are required\n");
        do_help();  // doesn't return
    }

    ParmName = Params.ParmName;
    p = get_parm("OutName");
    if (p->set) {
        ParmName = (char *)p->ptr;
    }
    p = get_parm("Unconstrained");
    uflg = p->set;
    p = get_parm("Extent");
    xflg = p->set;
    p = get_parm("Transform");
    tflg = p->set;

    // get data with sensor structures
    if (vflg) {
        printf("S A M   C O E F F I C I E N T S   F O R   R O I\t%s\n", PRG_REV);
        printf("opening data file '%s'", Header.SetName);
        fflush(stdout);
    }
    GetDsInfo(DSName, &Header, &OrigChannel, &Epoch, &Bad, TRUE);
    GetDsInfo(DSName, &Header, &Channel, &Epoch, NULL, TRUE);   // clone the channel structure
    GetFilePath(DSDIR, SAMpath, sizeof(SAMpath), &Params, "SAM", 0);
    if (access(SAMpath, F_OK) != 0)
        Cleanup("no SAM subdirectory; did you first execute sam_cov?");

    // derive constants
    M = Header.NumPri;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index

    ORDER = Params.Order;       // default should be 20th-order spherical harmonic expansion
    if (ORDER < 0) {
        ORDER = 20;
    }

    p = get_parm("MRIDirectory");
    if (!p->set)
        Cleanup("'MRIDirectory' path req'd");
    p = get_parm("ROIList");
    if (!p->set)
        Cleanup("'ROIList' req'd");
    p = get_parm("Surface");
    if (!p->set)
        Cleanup("'Surface' req'd");

    // sanity check for parameters

    if (Params.CovHP == -999. || Params.CovLP == -999.)
        Cleanup("CovBand bandpass parameters required");
    if (Params.CovHP < Channel[f].AnalogHpFreq)
        Cleanup("CovBand highpass cutoff frequency out of recorded range");
    if (Params.CovLP > Channel[f].AnalogLpFreq)
        Cleanup("CovBand lowpass cutoff frequency out of recorded range");
    if (uflg) {
        if (Params.OrientHP == -999. || Params.OrientLP == -999.)
            Cleanup("OrientBand bandpass parameters required");
        if (Params.OrientHP < Channel[f].AnalogHpFreq)
            Cleanup("OrientBand highpass cutoff frequency out of recorded range");
        if (Params.OrientLP > Channel[f].AnalogLpFreq)
            Cleanup("OrientBand lowpass cutoff frequency out of recorded range");
    }

    // announce surface & read cortical hull model
    if (vflg) {
        printf(" - done\n");
        if (!xflg && vflg) {
            msg("warning! -- computing beamformer coefficients for full ROI will be very slow!\n");
        }
        if (!uflg)
            printf("using '%s.norm' surface norms\n", Params.Surface);
        else
            printf("using linear solution for source orientation\n");
        printf("reading cortical hull model");
        fflush(stdout);
    }
    s = "hull.shape";
    p = get_parm("HullName");
    if (p->set) {
        s = (char *)p->ptr;
    }
    GetMRIPath(fpath, sizeof(fpath), &Params, s, TRUE);
    GetHull(fpath, &Hull);      // read hull model from MRI directory

    // read covariance file
    if (vflg) {
        printf(" - done\n");
        printf("reading covariance file(s)");
        fflush(stdout);
    }
    if (uflg) {
        sprintf(fpath, "%s/%s,%-d-%-dHz/Orient.cov", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
        Co = gsl_matrix_alloc(M, M);
        GetCov(fpath, &CovHeader, &ChanIndex, Co);
    } else {
        Co = NULL;
    }
    sprintf(fpath, "%s/%s,%-d-%-dHz/Global.cov", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
    C = gsl_matrix_alloc(M, M);
    GetCov(fpath, &CovHeader, &ChanIndex, C);

    // if transform is present, move all MEG primary and reference sensors
    if (tflg) {
        sprintf(fpath, "%s/%s.xfm", SAMpath, Params.Surface);
        if (vflg) {
            printf(" - done\n");
            printf("reading '%s.xfm' & transforming MEG sensor coordinates", Params.Surface);
            fflush(stdout);
        }
        MoveFrame(&Header, OrigChannel, Channel, fpath);
    }

    // read list of ROI targets
    if (vflg) {
        printf(" - done\n");
        printf("reading ROI list");
        fflush(stdout);
    }
    GetMRIPath(fpath, sizeof(fpath), &Params, Params.ROIList, 0);
    if ((fp = fopen(fpath, "r")) == NULL)
        Cleanup("can't open list of ROIs '%s'", fpath);
    V = 0;
    while (fgets(Line, sizeof(Line), fp) != NULL)   // count lines
        V++;
    if (V == 0)
        Cleanup("empty ROI file");
    rewind(fp);

    // allocate memory for coefficients
    if (vflg) {
        printf(" - done\n");
        printf("allocating memory");
        fflush(stdout);
    }
    W = gsl_matrix_alloc(V, M);
    Wv = gsl_vector_alloc(M);
    A = gsl_vector_alloc(M);

    // open log file for write
    sprintf(fpath, "%s,%s.log", ParmName, Params.ROIName);
    if ((fl = fopen(fpath, "w")) == NULL)
        Cleanup("can't open log file '%s' for write", fpath);
    fprintf(fl, "ROI\tRank\tNused\ts2\t\tn2\t\tz2\n");

    // compute ROI beamformer coefficients
    if (vflg) {
        printf(" - done\n");
        printf("computing beamformer coefficients for %d ROIs:\n", V);
    }

    s = strecpy(fpath, Params.Surface);
    strcpy(s, ".norm");
    GetMRIPath(NormPath, sizeof(NormPath), &Params, fpath, TRUE);
    for (c=0; c<V; c++) {

        // evaluate coefficients once
        CoeffFlag = FALSE;
        if (c == 0)
            CoeffFlag = TRUE;

        // read line from ROIfile
        if ((fscanf(fp, "%d%d", &ROInum, &ROIhemi)) != 2)
            cleanup("can't read ROI specs");
        if (vflg) {
            printf("\t%d ", ROInum);
            if (ROIhemi == 1)
                printf("L: ");
            else
                printf("R: ");
            fflush(stdout);
        }

        // read atlas & parse ROI vertices
        GetROIVertices(&Vertex, NormPath, ROInum, ROIhemi, &N);

        // solve weights based on entire ROI or radial extent
        if (vflg) {
            printf("solving weights");
            fflush(stdout);
        }
        if (xflg) {     // extent defined -- find ROI centroid & get vertices within extent

            // find ROI centroid
            for (v=X_; v<=Z_; v++)
                Vcentroid[v] = 0.;
            for (n=0; n<N; n++)
                for (v=X_; v<=Z_; v++)
                    Vcentroid[v] += Vertex[n].p[v];
            for (v=X_; v<=Z_; v++)
                Vcentroid[v] /= (double)N;
            GetPatchVertices(Vertex, N, Vcentroid, Params.Extent, &Nused);
            N = Nused;
        }

        // solve ROI weights
        if (N >= MIN_VERTEX) {
            ROISolveWts(Vertex, N, Wv, C, Co, &Header, Channel, &Hull, &rank, uflg);
            gsl_matrix_set_row(W, c, Wv);
            if (vflg) {
                printf(" - done\n");
            }
        } else {
            printf("\nempty ROI\n");
        }

        // solve s2 & n2
        gsl_blas_dgemv(CblasNoTrans, 1., C, Wv, 0., A);     // A = C * Wv
        gsl_blas_ddot(Wv, A, &s2);
        gsl_blas_ddot(Wv, Wv, &n2);
        n2 *= CovHeader.Noise;

        if (vflg) {
            printf("\t\trank=%d (%d): s2=%e, n2=%e, z2=%f\n", rank, N, s2, n2, s2/n2);
        }
        if (ROIhemi == -1)
            fprintf(fl, "%d -R\t%d\t%d\t%e\t%e\t%f\n", ROInum, rank, N, s2, n2, s2/n2);
        else
            fprintf(fl, "%d -L\t%d\t%d\t%e\t%e\t%f\n", ROInum, rank, N, s2, n2, s2/n2);
    }
    fclose(fp);     // close target list

    // save legacy beamformer coefficients (for DataEditor virtual sensors)
    sprintf(fpath, "%s/%s.wts", SAMpath, ParmName);
    for (c=0; c<V; c++)
        PutWts(fpath, W, Header, CovHeader);

    // save NIFTI coefficients
    sprintf(fpath, "%s/%s,%-d-%-dHz/Global.nii", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
    PutNIFTITargets(fpath, W, CovHeader.Noise);

    // finit
    log_params(SAMpath);

    gsl_matrix_free(C);
    if (Co)
        gsl_matrix_free(Co);
    fclose(fl);

    if (vflg) {
        printf("'%s' done\n", Progname);
    }
    exit(0);
}
