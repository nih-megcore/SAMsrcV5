// sam_param.h
//
//      Parameter labels:
//          Conditions & Time Windows:
//              NumMarkers <n>                  -- how many markers
//              Marker1 <marker name start end (sec) sum [name]>
//              Marker2 <marker name start end (sec) sum [name]>
//                  "
//              MarkerN <marker name start end (sec) sum [name]>
//              BlockName <marker>              -- block (epoch) of trials
//              Baseline <start end (sec)>      -- relative to Marker1
//              DataSegment <start end (sec)>   -- extended segment for processing RVE, spectra, power, etc.
//          Frequency Bands:
//              CovBand <higpass lowpass> (Hz)
//              ImageBand <highpass lowpass> (Hz)
//              SmoothBand <highpass lowpass> (Hz)
//              Notch <TRUE | FALSE>
//          ROI parameters:
//              XBounds <start end> (cm)
//              YBounds <start end> (cm)
//              ZBounds <start end> (cm)
//              ImageStep <n> (cm)
//              Atlas <atlas file>
//          Covariance Options:
//              CovType <GLOBAL | ALL | SUM
//              Mu <value (ft/rt-Hz)>
//              PropMu <value (unitless)>
//          Model <SPHERE | MULTI | NOLTE>
//          ImageMetric <POWER | RVE | KURTOSIS>
//          TimeStep <n> (sec)
//          TimeInt <n> (sec)
//          RemoveBaseline <NONE | BYVOXEL | GLOBAL>
//          SignSegment <start end> (sec)
//          MinSNR <snr>
//
//  Author: Stephen E. Robinson
//          Core MEG Group
//          NIMH
//

#ifndef H_PARAMS
#define H_PARAMS

// baseline removal
#define NONE        0               // don't remove baseline
#define BYVOXEL     1               // remove means for each voxel
#define GLOBAL      2               // remove global means from each voxel

// FreeSurfer surfaces
#define NO_SURF     -1
#define WM_SURF     0
#define SMWM_SURF   1
#define PIAL_SURF   2

// transform type
#define ORIG    0
#define TLRC    1

// filter type
#define FFT     0
#define IIR     1

// CovType
#define GLOBAL_ 0                   // global covariance
#define SUM_    1                   // sum of covariance for all markers
#define ALL_    10                  // separate covariance for each marker

// misc
#define X_      0
#define Y_      1
#define Z_      2

// 'MARKINFO' -- structure to hold marker conditions & segments

typedef struct {
    char        *MarkName;
    char        *MarkName2;         // NULL or a name to use instead of MarkName
    double      MarkStart;
    double      MarkEnd;
    int         Sum;                // TRUE | FALSE for sum covariance
    int         SegFile;            // TRUE if this is a SegFile
    char        *FileName;          // name of file containing segment times
} MARKINFO;

// Generic time range struct

typedef struct {
    double t0;                      // times in seconds
    double t1;
} TIMERANGE;

// ImageFormat buffer

typedef struct {
    char *fmt;
    double res;
} IMGFORMAT;

// Regularization parameter

typedef struct {
    double mu;                      // mu value
    int op;                         // operator
} MUINFO;

#define ADD_MU  0
#define MPY_MU  1
#define PROP_MU 2

// Head model, sphere center for single sphere, optional order for Nolte

typedef struct {
    int model;
    double sphere[3];
    int order;
} MODELINFO;

// models
#define SSPHERE     0               // single sphere
#define MSPHERE     1               // multiple local spheres
#define NOLTE       2               // Nolte model

// ImageMetric  @@@ this is messy

typedef struct {
    int metric;
    double tau;
    double lags;
    int dims;
} METRICINFO;

// imaging metrics
#define MOMENT      0               // image functions of source moment
#define POWER       1               // image functions of source power
#define HILBERT     2               // image functions of hilbert envelope
#define RV_ENTROPY  3               // image functions of rank vector entropy
#define S_ENTROPY   4               // image functions of spectral entropy
#define KURTOSIS    5               // image functions of excess kurtosis
#define MUT_INFO    6               // image functions of mutual information
#define ST_ENTROPY  7               // image functions of symbolic transfer entropy
#define RVC_ENTROPY 8               // image functions of conditional rank vector entropy

// 'PARMINFO' -- structure to hold analysis parameters

typedef struct {

    // names & directories
    char        *ParmName;          // parameter file name without .param
    char        *DataSetName;       // input dataset
    char        *DirName;           // output directory for image files
    char        *BlockName;         // marker name for block of trials
    int         NumPrefix;          // number of prefix characters
    char        PrefixChar;         // prefix delimiter character (when NumPrefix < 0)

    // conditions
    int         NumMark;            // number of markers
    MARKINFO    *Marker;            // marker info
    double      BaseStart;          // baseline window start
    double      BaseEnd;            // baseline window end
    double      DataSegStart;       // start time of segment relative to marker
    double      DataSegEnd;         // end time of segment relative to marker

    // frequency passbands
    double      CovHP;              // covariance high pass frequency
    double      CovLP;              // covariance low pass frequency
    double      ImageHP;            // image high pass frequency
    double      ImageLP;            // image low pass frequency
    double      OrientHP;           // orientation [Global] covariance high pass frequency
    double      OrientLP;           // orientation [Global] covariance low pass frequency
    double      NoiseHP;            // noise high pass frequency
    double      NoiseLP;            // noise low pass frequency
    double      SmoothHP;           // highpass for bandpass filtering output (0 if smoothing)
    double      SmoothLP;           // lowpass for smoothing imaging time-series
    int         FilterType;         // filter type (FFT | IIR)
    int         Notch;              // notch flag (TRUE | FALSE)
    double      Hz;                 // mains frequency

    // ROI parameters
    double      SAMStart[3];        // starting image bounds
    double      SAMEnd[3];          // ending image bounds
    double      SAMStep;            // image step size
    char        *Target;            // target file full path
    char        *TargetName;        // target file name
    char        *ROIList;           // ROI file full path
    char        *ROIName;           // ROI file name
    char        *SeedChannel;       // name of seed channel
    double      Extent;             // Extent of ROI for patch
    double      Interpolate;        // interpolation value for weights
    char        *Atlas;             // name of atlas containing dipole positions & normals
    char        *AtlasName;         // atlas basename for transform filename
    char 	*XformName;     // optional additional name for transform
    char        *MRIDirectory;      // %M path for MRI directory, see GetFilePath.c
    char        *MRIPattern;        // specification of hull, etc., locations
    char        *Mask;              // mask file name
    char        *HullName;          // hull file name
    int         ImageFormat;        // 0=orig, 1=tlrc
    double      ImageRes;           // output image/weight resolution (mm)

    // covariance options
    int         CovType;            // covariance type: GLOBAL, SUM, ALL
    double      Mu;                 // actually, added noise in fT per root Hz (replaced mu as a multiplier for sigma^2 * I)
    int         Operator;           // ADD_MU: '+' or MPY_MU: 'x'
    double      SignSegStart;       // start time of sign segment
    double      SignSegEnd;         // end time of sign segment

    // model parameters
    int         Model;              // conducting model type
    int         Order;              // model order for Nolte solution
    double      Sp[3];              // local sphere origin

    // coreg parameters
    int         MaxVertex;          // vertex sample size (for freesurfer segmentation)
    double      MinZ;               // z-axis cutoff for mesh vertices
    double      MinSpan;            // minimum eigenvalue ratio
    double      MaxAngle;           // maximum angular error
    double      MaxDot;             // maximum dot product for including vertices with near-radial normals
    double      MinSNR;             // minimum S/N
    double      DeltaTrans;         // delta for translation line-search or partial derivatives (m)
    double      DeltaRotat;         // delta for rotation for line-search or partial derivatives (rad)
    double      InitDisp;           // initial translation step-size (m/rad)
    double      EndDisp;            // stopping translation step-size (m/rad)
    double      Damping;            // sequential line-search damping factor
    int         LineSteps;          // sequential line-search steps (+/-)

    // image metric parameters
    int         ImageMetric;        // e.g., POWER, ENTROPY, KURTOSIS
    double      Tau;                // 1/e decay time (sec)
    double      Lags;               // sample delay between two symbolic time series
    int         Dims;               // embedding space dimension (default 4)

    // time dimension controls
    double      TimeStep;           // time increment per SAM image
    double      TimeInt;            // integration window duration

    // voxel normalization
    int         RemoveBaseline;     // flag for baseline removal (valid for SAMpower, SAMentropy, etc.)

} PARMINFO;

#endif
