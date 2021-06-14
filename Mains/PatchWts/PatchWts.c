// PatchWts -- compute LCMV beamformer coefficients for extended patches
//  of vertices & their normal vectors (from Freesurfer). Unlike 'ROIwts',
//  this version solves for a target coordinate & extent. If no target
//  file is given, 'PatchWts' solves for a rectangular array of voxels
//  (also using the vertices). Otherwise, the vertices are
//  found by the nearest neighbor to the targets. This version
//  uses the realistic head model (Nolte) forward solution!
//
//  flags:
//  -r <run name>           -- MEG Run name
//  -m <parameters>         -- parameter file
// @@@
//  -n <parameter name>     -- name for finding targets in ../Max subdirectory (optional)
//  -s                      -- solve moment vectors instead of using cortical normals
// @@@
//  -v                      -- Verbose mode
//
//  Author: Stephen E. Robinson
//      MEG Core Facility
//      NIMH
//

//#define HAVE_INLINE 1
//#define GSL_RANGE_CHECK_OFF 1

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
#include <nifti1.h>
#include <version.h>
#include "samutil.h"
#include "sam_param.h"
#include "sam_parse.h"

//#define MINOR_REV   1                   // optional name for finding targets in ../Max subdirectory
#define MINOR_REV   2                   // new command line parser
#define MAX_NOISE   1.0e-19             // maximum noise power

extern int CoeffFlag;                   // initialization for Nolte solution
extern int ORDER;                       // spherical harmonic expansion order

int main(
    int             argc,
    char            **argv
) {
    HeaderInfo      Header;             // MEG data header
    HULL            Hull;               // cortical hull model -- vertices & normals
    PARMINFO        Params;             // analysis parameters
    PARM            *p;                 // pointer into active parameter list
    VOXELINFO       *Vertex;            // Vertex[N] -- atlas vertex structure list
    VOXELINFO       *Target;            // Target[V] -- target list
    ChannelInfo     *OrigChannel;       // untransformed MEG channel info
    ChannelInfo     *Channel;           // cloned MEG channel info for transformation
    EpochInfo       *Epoch;             // MEG epoch info
    COV_HDR         CovHeader;          // covariance file header
    gsl_matrix      *C;                 // C - MxM source covariance matrix
    gsl_matrix      *Co;                // Co - MxM orientation covariance matrix
    gsl_matrix      *W;                 // W - VxM beamformer coefficients
    gsl_vector      *Wv;                // Wv - Mx1 beamformer coefficient vector
    double          Vcentroid[3];       // nearest-neighbor centroid
    double          min;                // minimum distance^2
    double          s2;                 // source power
    double          n2;                 // noise power
    double          snr;                // snr
    double          r2;                 // r^2
    double          vmin[3], vmax[3];   // bounding box
    double          percent;            // floating percentage
    double          tmp;                // what else?
    float           *power;             // power[V] -- image voxels for S/N
    int             *ChanIndex;         // ChanIndex[M] -- index of channels used for covariance
    int             rank;               // number of bases
    int             pct;                // integer percentage
    int             M;                  // number of primary channels
    int             N;                  // number of total vertices
    int             V;                  // number of targets
    int             Nuse;               // total number of vertices used
    int             i;                  // i-index
    int             j;                  // j-index
    int             n;                  // vertex index
    int             f;                  // index of 1st primary sensor (for recorded bandpass)
    int             v;                  // vector index
    int             eflg = FALSE;       // command-line error flag
    int             mflg = FALSE;       // parameter flag
    int             nflg = FALSE;       // prior parameter name flag
    int             rflg = FALSE;       // dataset flag
    int             tflg = FALSE;       // target file flag
    int             uflg = FALSE;       // unconstrained source vector flag
    int             vflg = FALSE;       // verbose mode flag (default: silent)
    char            *name;              // random filename
    char            *DSName = NULL;     // MEG dataset name
    char            CovPath[256];       // covariance file path name
    char            HullPath[256];      // hull path name
    char            SAMpath[256];       // data set path
    char            TargetPath[256];    // target (ROI) path name
    char            VertexPath[256];    // vertex file full path name
    char            fpath[256];         // general path name
    char            *ParmName = NULL;   // parameter file name
    char            PriorName[64];      // prior parameter name
    char            Prefix[64];         // dataset prefix code
    unsigned char   **Bad;              // Bad[E][T] -- flags for bad samples

    // Register the parameters used by this program.

    reg_usage(argv[0], MINOR_REV, __DATE__);
    reg_std_parm();

    reg_parm("Marker");

    reg_parm("CovBand");
    reg_parm("OrientBand");

    reg_parm("MRIDirectory");
    reg_parm("MRIPattern");
    reg_parm("PrefixLength");
    reg_parm("HullName");

    reg_parm("Model");
    reg_parm("doMultiSphere");  // since Model isn't allowed on the command line
    reg_parm("doNolte");        // allow --domulti or --donolte

    reg_parm("AtlasName");
    reg_parm("TargetName");
    reg_parm("ImageFormat");
    reg_parm("Extent");

    reg_parm("Transform");
    reg_parm("Surface");        // @@@ use this instead of GIFTI atlas??

    reg_parm("OutName");

    // Parse the arguments, environment variables, and parameter files.

    do_parse_args(argc, argv);

    // Extract specified parameters.

    eflg = get_params(&Params);

    p = get_parm("DataSet");
    if (p->set) {
        rflg = TRUE;
        DSName = copy_string(Params.DataSetName);   // modified!
    }
    p = get_parm("verbose");
    vflg = p->set;
    p = get_parm("param");
    mflg = p->set;

/*
    p = get_parm("Normalize");      // -n
    wflg = p->set;
    p = get_parm("Field");          // -B
    bflg = p->set;
    p = get_parm("Pinv");           // -p
    pflg = p->set;
    if (pflg) {
        ndim = *(int *)p->ptr;
    }
*/

    if (eflg || !rflg || !mflg) {
        msg("dataset (-r) and parameter file (-m) are required\n");
        do_help();  // doesn't return
    }

    ParmName = Params.ParmName;
    p = get_parm("OutName");
    if (p->set) {
       ParmName = (char *)p->ptr;
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

    // announce
    if (vflg) {
        printf("B E A M F O R M E R   C O E F F I C I E N T S   F O R   P A T C H E S\t%s\n", PRG_REV);
        printf("opening data file '%s'", DSName);
        fflush(stdout);
    }

    // get data with sensor structures
    name = copy_string(DSName);
    GetDsInfo(DSName, &Header, &OrigChannel, &Epoch, &Bad, TRUE);
    GetDsInfo(name, &Header, &Channel, &Epoch, NULL, TRUE);  // clone the channel structure
    GetFilePath(DSDIR, SAMpath, sizeof(SAMpath), &Params, "SAM", 0);
    if (access(SAMpath, F_OK) != 0)
        Cleanup("no SAM subdirectory; did you first execute SAMcov?");
    free(name);

    // set ORDER
    ORDER = Params.Order;       // default should be 20th-order spherical harmonic expansion
    if (ORDER < 0) {
        ORDER = 20;
    }

    // derive constants
    M = Header.NumPri;
    f = Header.PsIndex[0];      // 1st MEG primary sensor channel index

    p = get_parm("MRIDirectory");
    if (!p->set)
        fatalerr("MRIDirectory req'd");

    // set prefix to NumPrefix characters (default is 8-character dataset hashcode)
    memset(Prefix, 0, sizeof(Prefix));
    memcpy(Prefix, Header.SetName, Params.NumPrefix);

    // sanity check parameters
    if (Params.Extent == -999.)
        Cleanup("Extent specification required");
    if (Params.CovHP == -999. || Params.CovLP == -999.)
        Cleanup("CovBand bandpass parameters required");
    if (Params.CovHP < Channel[f].AnalogHpFreq)
        Cleanup("CovBand highpass cutoff frequency out of recorded range");
    if (Params.CovLP > Channel[f].AnalogLpFreq)
        Cleanup("CovBand lowpass cutoff frequency out of recorded range");
    if (uflg) {     // @@@  there is no uflg
        if (Params.OrientHP == -999. || Params.OrientLP == -999.)
            Cleanup("OrientBand bandpass parameters required");
        if (Params.OrientHP < Channel[f].AnalogHpFreq)
            Cleanup("OrientBand highpass cutoff frequency out of recorded range");
        if (Params.OrientLP > Channel[f].AnalogLpFreq)
            Cleanup("OrientBand lowpass cutoff frequency out of recorded range");
    }

    // read cortical hull model
    name = "hull.shape";
    p = get_parm("HullName");
    if (p->set) {
        name = (char *)p->ptr;
    }
    GetMRIPath(HullPath, sizeof(HullPath), &Params, name, TRUE);
    if (vflg) {
        printf(" - done\n");
        printf("reading cortical hull model");
        fflush(stdout);
    }
    GetHull(HullPath, &Hull);               // read hull model

    // read covariance file
    if (vflg) {
        printf(" - done\n");
        printf("reading covariance file '%s'", Params.Marker[0].MarkName);
        fflush(stdout);
    }
    if (uflg) {
        sprintf(CovPath, "%s/%s,%-d-%-dHz/Orient.cov", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
        Co = gsl_matrix_alloc(M, M);
        GetCov(CovPath, &CovHeader, &ChanIndex, Co);
    }
    if (Params.NumMark == 1)
        sprintf(CovPath, "%s/%s,%-d-%-dHz/%s.cov", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP, Params.Marker[0].MarkName);
    else
        sprintf(CovPath, "%s/%s,%-d-%-dHz/Global.cov", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
    C = gsl_matrix_alloc(M, M);
    GetCov(CovPath, &CovHeader, &ChanIndex, C);

    // if transform is present, move all MEG primary and reference sensors
    if (Params.Transform) {
        GetFilePath("%d/SAM/%s.xfm", fpath, sizeof(fpath), &Params, Params.Surface, TRUE);
        if (vflg) {
            printf(" - done\n");
            printf("reading '%s.xfm' & transforming MEG sensor coordinates", Params.Surface);
            fflush(stdout);
        }
        MoveFrame(&Header, OrigChannel, Channel, fpath);
    }

    @@@

    // read vertices from 'Surface'
    if (vflg) {
        printf(" - done\n");
        printf("reading vertices from '%s.norm'", Params.Surface);
        fflush(stdout);
    }
    sprintf(VertexPath, "%s/%s/%s.norm", Params.MRIdirectory, Prefix, Params.Surface);
    GetSurfer(VertexPath, &Vertex, &N);

    // check surface bounds
    for(v=X_; v<=Z_; v++) {
        vmin[v] = HUGE;
        vmax[v] = 0.;
    }
    for(n=0; n<N; n++) {
        for(v=X_; v<=Z_; v++) {
            if(Vertex[n].p[v] > vmax[v])
                vmax[v] = Vertex[n].p[v];
            if(Vertex[n].p[v] < vmin[v])
                vmin[v] = Vertex[n].p[v];
        }
    }
    if(vflg == TRUE) {
        printf("\nread %d vertices ranging from:\n", N);
        printf("\tx=%7.3f to %7.3f (cm)\n", 100.*vmin[X_], 100.*vmax[X_]);
        printf("\ty=%7.3f to %7.3f (cm)\n", 100.*vmin[Y_], 100.*vmax[Y_]);
        printf("\tz=%7.3f to %7.3f (cm)\n", 100.*vmin[Z_], 100.*vmax[Z_]);
        fflush(stdout);
    }

    // check for targets -- 1st, see if rectangular ROI bounds are specified either in Max subdirectory, named by TargetName, or rectangular ROI with x,y,z bounds
    if(Params.SAMStart[X_] != -999. && Params.SAMStart[Y_] != -999. && Params.SAMStart[Z_] != -999. && Params.SAMStep != -999.) {  // no need to test SAMend...
        tflg = FALSE;
    } else if(strncmp("NULL", Params.Target, 4)) {      // no.. then test for named target file & verify path
        sprintf(TargetPath, "%s/%s/%s", Params.MRIdirectory, Prefix, Params.Target);
        if(access(TargetPath, F_OK) != 0)
            Cleanup("can't access target file");
        tflg = TRUE;
    } else if(nflg == TRUE) {                           // finally, look for automatic target in Max subdirectory
        sprintf(TargetPath, "%s/Max/%s,%s.nii.max", SAMpath, Header.SetName, PriorName);
        if(access(TargetPath, F_OK) == 0)
            tflg = TRUE;
        else
            Cleanup("invalid target path name");
    } else {                                            // nothing found!
            Cleanup("no ROI or targets found");
    }

    // read target list
    if(tflg == TRUE) {
        if(vflg == TRUE) {
            printf("reading target list");
            fflush(stdout);
        }
        GetTargets(TargetPath, &Target, &V);
    } else {
        if(vflg == TRUE) {
            printf("generating voxel list");
            fflush(stdout);
        }
        GetVoxels(&Target, &V, &Hull, Params.SAMStart, Params.SAMEnd, Params.SAMStep);
    }

    // allocate memory for coefficients
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("allocating memory");
        fflush(stdout);
    }
    W = gsl_matrix_alloc(V, M);     // V is the number of targets!
    Wv = gsl_vector_alloc(M);       // M is the number of sensors
    if(tflg == FALSE)
        if((power = (float *)malloc(V * sizeof(float))) == NULL)
            allocfailed("power[]");

    // compute ROI beamformer coefficients
    if(vflg == TRUE) {
        printf(" - done\n");
        printf("computing beamformer coefficients for %d targets:\n", V);
        fflush(stdout);
    }
    for(c=0, pct=(-1); c<V; c++) {

        // evaluate spherical harmonic coefficients once
        if(c == 0)
            CoeffFlag = TRUE;
        else
            CoeffFlag = FALSE;

        if(tflg == TRUE) {

            for(n=0, min=HUGE; n<N; n++) {      // find nearest vertex centroid for target
                for(v=X_, r2=0.; v<=Z_; v++) {
                    tmp = Target[c].p[v] - Vertex[n].p[v];
                    r2 += tmp * tmp;
                }
                if(r2 < min) {
                    min = r2;
                    i = n;
                }
            }
            for(v=X_; v<=Z_; v++)
                Vcentroid[v] = Vertex[i].p[v];
        } else {
            for(v=X_; v<=Z_; v++)
                Vcentroid[v] = Target[c].p[v];
        }
        GetPatchVertices(Vertex, N, Vcentroid, Params.Extent, &Nuse, uflg);
        if(Nuse > 0) {

            // solve for optimal ROI beamformer coefficients
            if(vflg == TRUE && tflg == TRUE) {
                printf("target %d beamformer coefficients", c+1);
                fflush(stdout);
            }
            ROISolveWts(Vertex, Nuse, Wv, C, Co, &Header, Channel, &Hull, &rank, uflg);
            if(vflg == TRUE && tflg == TRUE) {
                printf(" - done\n");
                fflush(stdout);
            }

            // solve source power and projected noise power
            if(gsl_vector_isnull(Wv) == 0) {
                for(i=0, s2=0.; i<M; i++)
                    for(j=0; j<M; j++)
                        s2 += gsl_vector_get(Wv, i) * gsl_matrix_get(C, i, j) * gsl_vector_get(Wv, j);
                for(i=0, n2=0.; i<M; i++)
                    n2 += gsl_vector_get(Wv, i) * CovHeader.Noise * gsl_vector_get(Wv, i);
                snr = (n2 > 0.)? (s2 - n2) / n2: 0.;
                if(snr <= 0.)       // set coefficients to zero if (s-n)/n <= 0
                    gsl_vector_set_zero(Wv);
            } else {
                s2 = n2 = snr = 0.;
            }
            gsl_matrix_set_row(W, c, Wv);
        } else {
            if(tflg == TRUE)
                fprintf(stderr, "too few vertices for target %d", c);
        }
        if(tflg == TRUE) {
            if(vflg == TRUE) {
                printf("\trank=%d (%d): s2=%e, n2=%e, z2=%f\n", rank, Nuse, s2, n2, snr);
                fflush(stdout);
            }
        } else {
            if(Nuse > 0 && snr > 0.) {      // (s-n)/n > 0
                if(vflg == TRUE) {
                    printf("c=%d: rank=%d (%d): s2=%e, n2=%e, z2=%f\n", c, rank, Nuse, s2, n2, snr);
                    fflush(stdout);
                }
            }
            tmp = (n2 > 0.)? 1.0e+18 * (s2 - n2): 0.;
            power[c] = (Nuse > 0 && n2 < MAX_NOISE)? (float)tmp: 0.;        // image (nA-m)^2
        }

        // percent done
        if(vflg == TRUE && tflg == FALSE) {
            percent = 100. * (double)c / (double)V;
            if(rint(percent) > (double)pct) {
                pct = (int)rint(percent);
                printf("\r\t%d percent done", pct);
                fflush(stdout);
            }
        }

    }

    // save legacy beamformer coefficients
    if(nflg == TRUE)
        sprintf(fpath, "%s/%s,%s.wts", SAMpath, ParmName, PriorName);
    else
        sprintf(fpath, "%s/%s,%s.wts", SAMpath, ParmName, Params.Target);
    PutWts(fpath, W, Header, CovHeader);

    // save NIFTI coefficients
    sprintf(fpath, "%s/%s,%-d-%-dHz/Global.nii", SAMpath, ParmName, (int)Params.CovHP, (int)Params.CovLP);
    if(tflg == TRUE) {
        PutNIFTITargets(fpath, W, CovHeader.Noise);
    } else {
        PutNIFTIWts(fpath, &Params, W, CovHeader.Noise);
        sprintf(fpath, "%s/%s/%s,Image.nii", Params.MRIdirectory, Prefix, ParmName);
        PutNIFTIImage(power, fpath, &Params, "power");
    }

    // finit
    gsl_matrix_free(W);
    gsl_vector_free(Wv);
    if(tflg == FALSE)
        free((void *)power);
    if(vflg == TRUE) {
        printf("\n'PatchWts' done\n");
        fflush(stdout);
    }
    exit(0);
}
