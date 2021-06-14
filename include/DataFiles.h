// DataFiles.h -- proposed data file format
//
//      Author: Stephen E. Robinson
//              MEG Core Facility
//              NIMH
//

#ifndef H_DATAFILES
#define H_DATAFILES

#include <stdint.h>
#include <geoms.h>              // contains sensor structure definitions

#define VERSION         1       // this is, after all, Version one...

// order for virual memory data arrays
#define BY_CHANNEL      0
#define BY_EPOCH        1

// 'ChannelType' flag values
#define TYPE_MEG        1
#define TYPE_EEG        2
#define TYPE_REF_MAG    3
#define TYPE_REF_GRAD   4
#define TYPE_TRIGGER    5
#define TYPE_EXTERNAL   6
#define TYPE_OTHER      7

// 'Units' flag values
#define UNIT_TESLA      1
#define UNIT_VOLTS      2
#define UNIT_BITS       3
#define UNIT_OTHER      4

// 'BalInfo' structure for synthetic gradiometers
typedef struct {
    double  Weight;             // balance weight
    int     RefChan;            // reference channel index
} BalInfo;

// 'HeaderInfo' structure is at the beginning of the data file
typedef struct {
    short   Version;            // file version # (=VERSION!)
    char    *DsPath;            // full pathname of dataset
    char    *SetName;           // subj_study_date_run string
    int     NumChannels;        // total number of channels available
    int     NumAcquired;        // number of channels acquired
    int     NumSquid;           // number of SQUID sensors (primary + reference)
    int     *SqIndex;           // SqIndex[NumSquid] -- all primary & reference channel indices
    int     NumPri;             // number of primary MEG channels
    int     *PsIndex;           // PsIndex[NumPri] -- primary sensor channel index
    int     NumRef;             // number of reference channels
    int     *RsIndex;           // RsIndex[NumRef] -- reference sensor channel index
    int     NumEEG;             // number of EEG channels
    int     *PeIndex;           // PeIndex[NumEEG] -- primary EEG electrode index
    int     NumEpochs;          // number of epochs
    int     MaxSamples;         // maximum points per channel per epoch
    int     PreTrigSamples;     // number of pretrigger samples per epoch
    double  SampleRate;         // samples/second (Hz)
    double  PowerFreq;          // local power mains frequency (Hz)
    int     Avg;                // is avg?
    short   DataFormat;         // short, long, float, or double
    int     spillt;             // number of trials before split, or zero
    int32_t *trialbuf;          // temp buffer for one trial of one channel
} HeaderInfo;

// there will be an array, ChannelInfo[NumChannels], following the 'Header'
typedef struct {
    short   ChannelType;        // TYPE_MEG, TYPE_EEG, etc.
    short   Flag;               // flag for editing, etc.
    short   ChannelNumber;      // physical channel in a/d multiplexed list order
    short   Units;              // physical quantity measured -- UNIT_TESLA, UNIT_VOLTS, etc.
    double  Scale;              // multiplier for data (when using integer types, only)
    double  AnalogHpFreq;       // high-pass analog filter (Hz)
    double  AnalogLpFreq;       // low-pass analog filter frequency (Hz)
    double  HeadR2;             // head radius^2
    char    ChannelName[16];    // channel name
    union {
        SENSOR      MEGSensor;  // MEG sensor geometry, including balance coefficients
        ELECTRODE   EEGSensor;  // EEG sensor geometry structure
        OPM         OPMSensor;  // OPM sensor geometry structure
        short       NoGeom;     // 'TRUE' if no sensor geometry (as in either MEG or EEG)
    } Geom;
    short   BalOrder;           // balance order: 0, 1, 2, 3 (for real or ideal); 10, 11, 12, 13 (for adaptive coefficients)
    BalInfo *Balance;           // Balance[NumBal] -- balance information
    int     NumBal;             // number of balance coefficients
} ChannelInfo;

// each trigger marker has one MarkInfo structure (new in April 2009)
typedef struct {
    char    MarkName[32];       // marker name string
    int     MarkCode;           // 32-bit trigger marker code (for CTF and 4D)
    int     MarkSample;         // sample index for rising edge of trigger marker
    double  MarkTime;           // time for rising edge of trigger marker
} MarkInfo;

// there will be one EpochInfo structure for each epoch in the dataset
typedef struct {
    int     NumSamples;         // total number of samples per epoch per channel
    MarkInfo *Marker;           // Marker[N] -- trigger marker information
    int     NumMarkers;         // number of trigger markers
    short   Flag;               // for editing epoch (TRUE or FALSE)
} EpochInfo;

typedef struct {
    double  Weight;
    long    ChanIndex;
} WtItem;

typedef struct {
    FIELD   voxel;              // voxel coordinate & unit current vector
    WtItem  *Wts;               // Wts[M] -- SAM weight vector w/ indices
    int     M;                  // number of weights
} WtInfo;

#endif // H_DATAFILES
