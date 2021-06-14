#include "samutil.h"
#include "sam_param.h"
#include "sam_parse.h"
#include "coreg.h"

/* This table describes all the parameters used by all the programs. */

PARM_TABLE Parm_table[] = {
    { 1, "help", "h", helpfn, NULL, "show this help", NULL },
    { 1, "verbose", "v", boolfn, NULL, "verbose output", NULL },
    { 1, "DataSet", "r", stringfn, "DSNAME", "MEG dataset name", "ds" },
    { 1, "PDFName", "d", stringfn, "PDFNAME", "MEG PDF file (4D only)", "pdf" },
    { 1, "param", "m", stringfn, "PFILE", "parameter file name (optionally ending\nin \".param\")", "param" },
    { 0, "%include", NULL, includefn, "", "", NULL },

    { 1, "NumMarkers", NULL, intfn, "N", "number of markers (deprecated)", NULL },
    { 0, "Marker", NULL, markerfn, "MARKNAME T0 T1 SUMFLAG [COVNAME]", "marker name, time window relative to marker,\nwhether to include this marker in the\nSUM covariance (TRUE|FALSE), and an optional name\nto use instead of the marker name", NULL },

    //{ 1, "MarkFile", NULL, markfilefn, "MARKNAME FILE T0 T1 SUMFLAG", "file containing 'trial time' lines for a mark", NULL },
    { 1, "SegFile", NULL, segfilefn, "MARKNAME FILE SUMFLAG", "file containing 'trial T0 T1' lines for a segment", NULL },

    { 1, "Baseline", NULL, double2fn, "T0 T1", "baseline time window relative to Marker", NULL },
    { 1, "DataSegment", NULL, double2fn, "T0 T1", "allows expanding time windows relative\nto all markers", NULL },
    { 1, "SignSegment", NULL, double2fn, "T0 T1", "polarity time window relative to Marker", NULL },
    { 1, "OutRange", NULL, trangefn, "T0 T1", "output time range", NULL },

    { 1, "FilterType", NULL, stringfn, "FFT|IIR", "Fast Fourier (default) or Infinite Impulse Response", NULL },
    { 1, "CovBand", NULL, double2fn, "LO HI", "covariance bandwidth limits in Hz", NULL },
    { 1, "ImageBand", NULL, double2fn, "LO HI", "imaging covariance bandwidth limits in Hz", NULL },
    { 1, "OrientBand", NULL, double2fn, "LO HI", "band to use for orientation instead of Global (Hz)", NULL },
    { 1, "NoiseBand", NULL, double2fn, "LO HI", "bandwidth limits for estimating noise (Hz)", NULL },
    { 1, "SmoothBand", NULL, double2fn, "LO HI", "used to smooth Hilbert envelope, kurtosis, or RVE (Hz)", NULL },

    { 1, "Notch", NULL, boolfn, NULL, "apply a powerline notch\n(including harmonics)", NULL },
    { 1, "Hz", NULL, doublefn, "HZ", "frequency for powerline (Hz, default 60)", NULL },

    { 1, "XBounds", NULL, double2fn, "START END", "ROI X bounds in cm (positive is anterior)", NULL },
    { 1, "YBounds", NULL, double2fn, "START END", "ROI Y bounds in cm (positive is left)", NULL },
    { 1, "ZBounds", NULL, double2fn, "START END", "ROI Z bounds in cm (positive is superior)", NULL },
    { 1, "ImageStep", NULL, doublefn, "STEP", "image voxel size (cm)", NULL },
    //{ 1, "Interpolate", NULL, doublefn, "DIST", "interpolation distance (cm) for weights", NULL },

    { 1, "TargetName", "t", stringfn, "TARGETFILE", "file containing target coordinates", NULL },
    { 1, "ROIList", NULL, stringfn, "ROILISTFILE", "file containing regions of interest", NULL },
    { 1, "AtlasName", "A", stringfn, "ATLASFILE", "file containing target coordinates and orientations", NULL },
    { 1, "XformName", NULL, stringfn, "XFORMNAME", "Optional name root for transform file", NULL },
    { 1, "Transform", NULL, stringfn, "XFORM", "apply the transform in $ds/SAM/XFORM.xfm", NULL },

    { 1, "MRIDirectory", "i", stringfn, "MRIDIR", "MRI directory root", "mridir" },
    { 1, "ImageDirectory", "o", stringfn, "IMAGEDIR", "image output directory", "imagedir" },
    { 1, "HullName", "H", stringfn, "HULL", "hull name, default 'hull.shape'", NULL },
    { 1, "MRIPattern", NULL, stringfn, "PATTERN", "MRI pathname pattern, default '%M/%P/%s'", NULL },
    { 1, "PrefixLength", NULL, prefixfn, "N", "DataSet prefix may be specified as a number\nor as a prefix delimiter (default \"_\")", NULL },
    { 1, "CovName", "C", stringfn, "NAME", "use this name for covariance file names", NULL },
    { 1, "WtsName", "W", stringfn, "NAME", "use this name for weights file names", NULL },
    { 1, "OutName", "N", stringfn, "NAME", "use this name instead of the parameter file name\nin output file names", NULL },

    { 0, "ImageFormat", NULL, imgformatfn, "ORIG|TLRC RES", "original (default) or TLRC transform with\nRES (mm) sized voxels", NULL },
    { 1, "MaskName", NULL, stringfn, "MASKFILE", "mask filename for TLRC transform", NULL },

    { 1, "Normalize", "n", boolfn, NULL, "normalize the output", NULL },
    { 1, "ZNormalize", "Z", boolfn, NULL, "convert the output to z-scores", NULL },
    { 1, "Field", "B", boolfn, NULL, "output primary sensor forward solutions\n(along with weights)", NULL },
    { 1, "Read_binary", "b", boolfn, NULL, "read binary instead of recomputing", NULL },
    { 1, "Absolute", "a", boolfn, NULL, "take absolute value of voxels", NULL },
    { 1, "AUC", "u", boolfn, NULL, "area under the curve mode", NULL },

    { 1, "Unconstrained", "s", boolfn, NULL, "compute unconstrained moment vector", NULL },
    { 1, "Extent", NULL, doublefn, "DIST", "radial extent from ROI voxels or centroid (mm)", NULL },

    { 1, "CovType", NULL, stringfn, "GLOBAL|SUM|ALL", "which covariance matrix to use for the analysis", NULL },

    { 1, "Mu", "u", mufn, "[+|*]MU", "regularization adds a constant (+) or scaled (*)\nnoise to the diagonal of the covariance matrix", NULL },

    { 0, "Model", NULL, modelfn, "SingleSphere x y z|MultiSphere|Nolte", "forward solution type", NULL },
    { 0, "ImageMetric", NULL, imgmetricfn, "METRICSPEC", "imaging metric", NULL },
    { 1, "doMultiSphere", NULL, boolfn, NULL, "set Model to MultiSphere", NULL },
    { 1, "doNolte", NULL, boolfn, NULL, "set Model to Nolte", NULL },
    { 1, "Order", NULL, intfn, "ORDER", "spherical harmonic order for Nolte solution", NULL },

    { 1, "Pinv", NULL, intfn, "NDIM", "use a pseudo-inverse that removes the smallest\nNDIM dimensions", NULL },

    { 1, "TimeStep", NULL, doublefn, "STEP", "time step", NULL },
    { 1, "TimeInt", NULL, doublefn, "TIME", "integration time", NULL },
    { 1, "BoxWidth", NULL, doublefn, "TIME", "time window for boxcar smoothing", NULL },

    { 1, "RemoveBaseline", NULL, stringfn, "NONE|BYVOXEL|GLOBAL", "baseline removal method", NULL },
    { 1, "BlockName", NULL, stringfn, "MARKNAME", "epoch marker name", NULL },

    // These are for NIFTIPeak
    { 1, "Threshold", "t", doublefn, "THRESH", "image threshold", NULL },
    { 1, "SearchRadius", "s", doublefn, "RADIUS", "search radius in cm", NULL },
    { 1, "NumMaxima", "n", intfn, "N", "number of maxima to list (default 10)", NULL },
    { 1, "ShowPower", "p", boolfn, NULL, "whether to list power value", NULL },

    // These are for SAMcoreg
    { 1, "NumVerts", NULL, intfn, "NVERT", "number of vertices to use per pass", NULL },
    { 1, "MinZ", NULL, doublefn, "Z", "only consider vertices above this z-axis cutoff (cm)", NULL },
    { 1, "MinSpan", NULL, doublefn, "SPAN", "minimum eigenvalue span", NULL},
    { 1, "MaxAngle", NULL, doublefn, "ANGLE", "maximum angular error (rad)", NULL },
    { 1, "MaxDot", NULL, doublefn, "THRESH", "threshold to remove vertices with near-radial normals (0..1)", NULL },
    { 1, "MinSNR", NULL, doublefn, "SNR", "minimum signal/noise ratio", NULL },
    { 1, "ShowPasses", "p", boolfn, NULL, "output scatterplot files", NULL },
    { 1, "DeltaTranslate", NULL, doublefn, "TRANS", "line-search & partial derivative translation delta (mm)", NULL },
    { 1, "DeltaRotate", NULL, doublefn, "ANGLE", "line-search & partial derivative rotation delta (degrees)", NULL },
    { 1, "InitDisplacement", NULL, doublefn, "TRANS", "initial translation step-size (mm/mrad)", NULL },
    { 1, "EndDisplacement", NULL, doublefn, "TRANS", "ending translation step-size (mm/mrad)", NULL },
    { 1, "Offset", "o", double3fn, "x y z", "optional translation offset for testing convergence (mm)", NULL },
    { 1, "Damping", NULL, doublefn, "DAMPING", "sequential line-search damping factor", NULL},
    { 1, "LineSteps", NULL, intfn, "NSTEP", "sequential line-search steps (+/-)", NULL},

    { 0, NULL, NULL, NULL, NULL, NULL, NULL }
};

/* Fill in a SAM style PARMINFO structure with default settings. */

void new_params(PARMINFO *Params)
{
    int i;

    // initialize parameters
    Params->DataSetName = NULL;             // argument to -r
    Params->ParmName = NULL;                // basename of parameter file
    Params->NumMark = 0;                    // number of markers
    Params->Marker = NULL;                  // marker info
    Params->BaseStart = -999.;              // baseline window start
    Params->BaseEnd = -999.;                // baseline window end
    Params->DataSegStart = -999.;           // start time of segment relative to marker
    Params->DataSegEnd = -999.;             // end time of segment relative to marker
    Params->FilterType = FFT;               // filter type: FFT (default) or IIR
    Params->CovHP = -999.;                  // covariance high pass frequency
    Params->CovLP = -999.;                  // covariance low pass frequency
    Params->ImageHP = -999.;                // image high pass frequency
    Params->ImageLP = -999.;                // image low pass frequency
    Params->OrientHP = -999.;               // orientation covariance high pass frequency
    Params->OrientLP = -999.;               // orientation covariance low pass frequency
    Params->NoiseHP = -999.;                // noise high pass frequency
    Params->NoiseLP = -999.;                // noise low pass frequency
    Params->SmoothHP = -999.;               // image time-series smoothing high pass frequency
    Params->SmoothLP = -999.;               // image time-series smoothing low pass frequency
    Params->Notch = FALSE;                  // notch filter default is off
    Params->Hz = 60.;                       // mains frequency
    Params->SAMStart[X_] = -999.;           // X - starting image bounds
    Params->SAMStart[Y_] = -999.;           // Y
    Params->SAMStart[Z_] = -999.;           // Z
    Params->SAMEnd[X_] = -999.;             // X - ending image bounds
    Params->SAMEnd[Y_] = -999.;             // Y
    Params->SAMEnd[Z_] = -999.;             // Z
    Params->SAMStep = -999.;                // image step size
    Params->Target = "NULL";                // target file name
    Params->Extent = -999.;                 // radial extent from ROI centroid
    Params->Interpolate = -999.;            // interpolatation step size
    Params->Atlas = "NULL";                 // initialize atlas file name
    Params->MRIDirectory = "mri";           // initialize MRI directory pathname
    Params->Mask = "NULL";                  // initialize mask file name
    Params->ImageFormat = 0;                // image output format: 0=orig, 1=tlrc (default orig)
    Params->ImageRes = 0;                   // output image/weight resolution (tlrc, only!)
    Params->CovType = -1;                   // covariance type (GLOBAL | SUM | ALL)
    Params->Mu = -1.;                       // noise for regularization (-1 means auto)
    Params->Operator = -1;                  // regularization op, add or multiply
    Params->SignSegStart = -999.;           // start time of sign segment
    Params->SignSegEnd = -999.;             // end time of sign segment
    Params->Model = -1;                     // forward model for magnetic fields
    for (i = X_; i <= Z_; i++)              // sphere origin
        Params->Sp[i] = -999.;
    Params->Order = 16;                     // LFU spherical harmonic decomposition order
    Params->ImageMetric = -1;               // image metric (POWER | ENTROPY | MUTUAL INFORMATION | KURTOSIS)
    Params->Tau = -999.;                    // 1/e decay time (sec)
    Params->Lags = 0.;                      // lags between symbolic time-series
    Params->Dims = 4;                       // embedding space dimension (default 4)
    Params->TimeStep = -999.;               // time increment per SAM image
    Params->TimeInt = -999.;                // integration window duration
    Params->RemoveBaseline = NONE;          // default is no removal of baseline
    Params->DirName = "NULL";   // %d/SAM   // initialize output directory path
    Params->BlockName = "NULL";             // initialize marker name for block of trials
    Params->NumPrefix = -1;                 // default uses a delimiter character
    Params->PrefixChar = '_';               // use '_' to delimit the prefix
    Params->MRIPattern = "%M/%P/%s";        // default location of MRI files (see GetFilePath.c)
                                            // defaults for SAMcoreg (see coreg.h)
    Params->MaxVertex = MAX_VERTEX;         // number of verts per sample
    Params->MinZ = MIN_Z;                   // Z axis cutoff (m)
    Params->MinSpan = MIN_SPAN;             // minimum eigenvalue span for solving moment vector
    Params->MaxAngle = MAX_ANGLE;           // max angular error (rad)
    Params->MaxDot = MAX_DOT;               // threshold for excluding near-radial normals
    Params->MinSNR = MIN_SNR;               // minimum signal/noise
    Params->DeltaTrans = DELTA_TRANS;       // line-search & partial derivative translation delta (metres)
    Params->DeltaRotat = DELTA_ROTAT;       // line-search & partial derivative rotation delta (radians)
    Params->InitDisp = INIT_TRANS;          // initial translation step-size (metres)
    Params->EndDisp = END_TRANS;            // ending translation step-size (metres)
    Params->Damping = DAMPING;              // sequential line-search damping factor
    Params->LineSteps = NSTEP;              // sequential line-search steps (+/-)
}
