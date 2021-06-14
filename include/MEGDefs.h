/* ==========================================================================
 * $Header: /data/dougmck/code/ctf/ctfutil/include/RCS/MEGDefs.h,v 13.1 2000/06/20 14:45:17 dougmck Exp $
 * ==========================================================================
 *
 * MEGDefs.h
 *
 * Definitions used by the dataset and its clients.
 *
 * Copyright (c) CTF Systems Inc., 1995-1996.  All Rights Reserved.
 * Commercially Confidential Information
 *
 * ==========================================================================
 *
 * $Log: MEGDefs.h,v $
 * Revision 13.1  2000/06/20 14:45:17  dougmck
 * Merged in changes from release 4.12.2 ( HSC )
 *
 * Revision 6.0  1999/11/03 23:01:13  dnm
 * Initial check in for Linux/HP-UX 10 Code
 *
 * Revision 10.6  1999/05/14  21:13:26  sixtus
 * Added VirtualSensors.
 *
 * Revision 10.5  98/03/13  12:33:40  12:33:40  sixtus (Sixtus Lee)
 * Added eSAMSensor and SAMSensor to the enums of SensorType and SensorClass,
 * respectively.
 *
 * Revision 10.4  97/07/30  14:17:59  14:17:59  tmetzger ()
 * Removed a lot of dead code.
 * Enhanced the meg4GeneralResRec to handle DSQ800 and DSQ2000 triggers.
 *
 * Revision 10.3  97/04/18  12:01:31  12:01:31  murray (Murray)
 * added commercially confidential comment
 *
 * Revision 10.2  96/11/14  09:56:47  09:56:47  murray (Murray)
 * add definitions for G1OI
 *
 * Revision 10.1  96/02/01  16:34:15  16:34:15  qabbany (Moustafa Elqabbany)
 * Added standard headers and fixed the point class
 *
 * Revision 1.1  95/11/30  16:52:17  16:52:17  tmetzger (Tom)
 * Initial revision
 *
 *
 * ==========================================================================
 */

#ifndef _H_MegDefs
#define _H_MegDefs

#include <CTFStdTypes.h>
#include <MEGTypes.h>

#ifdef __GNUC__
#define __PACKED__ __attribute__ ((aligned( 2 )))
#else
#define __PACKED__
#endif


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define MAX_COILS			8
#define SENSOR_LABEL		31
#define MAX_NUM_COEFS		50
#define MAX_AVERAGE_BINS	8
#define MAX_BALANCING		MAX_NUM_COEFS
#define GENERALRESID		30000

#define G1BRINDEX			1		// define index for the coefficients
#define G2BRINDEX			2
#define G3BRINDEX			3
#define G2OIINDEX			4
#define G3OIINDEX			5
#define EDDYINDEX			6
#define G1OIINDEX			7

// SensorType
typedef enum {
	eMEGReference,
	eMEGReference1,
	eMEGReference2,
	eMEGReference3,
	eMEGSensor,
	eMEGSensor1,
	eMEGSensor2,
	eMEGSensor3,
	eEEGRef,
	eEEGSensor,
	eADCRef,
	eADCAmpRef = eADCRef,
	eStimRef,
	eTimeRef,
	ePositionRef,
	eDACRef,
	eSAMSensor,
	eVirtualSensor,
	eSystemTimeRef,
	eADCVoltRef,
	eStimAnalog,
	eStimDigital,
	EEGBipolar,
	eEEGAflg,
	eMEGReset,
	eDipSrc,
	eSAMSensorNorm,
	eOtherRef,
	eInvalidType
} SensorType;


// eChType
typedef enum {
	MEGRef,
	MEGSensor,
	EEGRef,
	EEGSensor,
	ADCRef,
	StimRef,
	TimeRef,
	PositionRef,
	DACRef,
	SAMSensor,
	VirtualSensor,
	badMEGSensor,
	badEEGSensor,
	ADCVoltRef,
	SuppRef,
	DipoleRef,
	OtherRef,
	InvalidClass
} eChType;

typedef eChType SensorClass;

// coiltype
typedef enum {
	CIRCULAR,
	SQUARE,
	VOLTMETER = CIRCULAR,
	AMMETER = SQUARE
} coiltype;

typedef coiltype CoilType;

// d3_point
typedef union d3_point {
	struct			{ DDouble x,y,z, junk;} c;
	struct			{ DDouble r,theta,phi, junk ;} s;
	DDouble			point[4];
} d3_point;


// d2_point
typedef union d2_point {
	struct			{ DDouble x,y; } c;
	struct			{ DDouble r,theta; } p;
	DDouble			point[2];
} d2_point;


// d3_point_ext
typedef union d3_point_ext {	// Externally store points as SDouble
	struct			{ SDouble x,y,z, junk ;} c;
	struct			{ SDouble r,theta,phi, junk ;} s;
	SDouble			point[4];
} d3_point_ext;


// CoilRec_ext
typedef struct CoilRec_ext {
	d3_point_ext	position;	// position of coil
	d3_point_ext	orient;		// orientation of coil
	Int16			numturns;	// number of turns making up the coil
	short			reserved1;	// pad out to the next 8 byte boundary
	short			reserved2;
	short			reserved3;
	SDouble			area;		// area of coil
} CoilRec_ext __PACKED__ ;


// CoilRec
typedef struct CoilRec {
	d3_point	position; 		// position of coil
	d3_point	orient;			// orientation of coil
	Int16		numturns;		// number of turns making up the coil
	DDouble		area;			// area of coil
} CoilRec __PACKED__ ;


// Coef_List
typedef struct coef_List {
    Int16		index;
    CChar		name[SENSOR_LABEL];
} Coef_List __PACKED__;

// NewSensorResRec
typedef struct {
	int16_t		sensorTypeIndex;
	int16_t		originalRunNum;
	CoilType	coilShape;
	double		properGain;			// may be corrected
	double		qGain;
	double		ioGain;
	double		ioOffset;
	int16_t		numCoils;
	int16_t		grad_order_no;
	int32_t		stimPolarity;		// new in version 4.2
    CoilRec_ext coilTbl[MAX_COILS];
    CoilRec_ext HdcoilTbl[MAX_COILS];
} NewSensorResRec __PACKED__ ;


typedef struct CoefResRec   /* Making generic resource for coefficients. */
{
    Int16           num_of_coefs;
    CChar           sensor_list[MAX_BALANCING][SENSOR_LABEL];
    SDouble         coefs_list[MAX_BALANCING];
} __PACKED__ CoefResRec, *CoefResRecP, **CoefResRecH;


// meg4FileSetup
typedef struct {
	char_t		nf_run_name[32];
	char_t		nf_run_title[256];
	char_t		nf_instruments[32];
	char_t		nf_collect_descriptor[32];
	char_t		nf_subject_id[32];
	char_t		nf_operator[32];
	char_t		nf_sensorFileName[60];	// !! docs say 56 but Acq says 60
	int32_t		size;					// length of following array
} __PACKED__ meg4FileSetup;

//    Int32   reserved1;                  /* pad out to the next 8 byte boundary */
//    CStrPtr nf_run_descriptor;
//    Int32   reserved2;                        /* can't have a pointer here */
//} __PACKED__ meg4FileSetup ;

typedef enum { CLASSERROR, BUTTERWORTH } classType;
typedef enum { TYPERROR, LOWPASS, HIGHPASS, NOTCH } filtType;

typedef struct
{
    SDouble     freq;
    classType   fClass;
    filtType    fType;
    Int16       numParam;
    SDoubleArr  params;
} __PACKED__ filter;

/*
 * This enum is used by GeneralRsrc to keep track of which
 * part of the trigger format union it will be reading
 */
enum TriggerStructFormat {
	    MEG40_TRIG_FMT,
	    MEG41_TRIG_FMT,
	    MEG42_TRIG_FMT
	};

/*
 * Trigger structure for the meg4 dataset
 */
typedef struct
{
	UCChar primaryTrigger;
	UCChar secondaryTrigger[MAX_AVERAGE_BINS];
	UCChar triggerPolarityMask;
} __PACKED__ meg40TriggerData;

/*
 * Trigger structure for the meg5 dataset
 */
typedef struct
{
	Bit32   primaryTrigger;
	Bit32   triggerPolarityMask;
} __PACKED__ meg41TriggerData;

typedef struct {
	int32_t		no_samples;
	int16_t		no_channels;
	int16_t		reserved1;		// pad out to the next 8 byte boundary
	double		sample_rate;
	double		epoch_time;
	int16_t		no_trials;
	int16_t		reserved2;		// doesn't pad out to the next 8 byte boundary
	int32_t		preTrigPts;
	int16_t		no_trials_done;
	int16_t		no_trials_display;
	CTFBoolean	save_trials;
//        union
//        {
//            meg40TriggerData meg40trig;
//            meg41TriggerData meg41trig;
//        };
	meg41TriggerData meg41trig;	// the union gets padded funny on x86_64, & we never use meg40 anymore.
	int16_t		reserved6;		// pad to 10 bytes after meg41trig
	int16_t		reserved3;		// pad out to the next 8 byte boundary
	int16_t		trigger_mode;
	int16_t		reserved4;		// pad out to the next 8 byte boundary
	CTFBoolean  accept_reject_Flag;
	int16_t		run_time_display;
	int16_t		reserved5;		// pad out to the next 8 byte boundary
	CTFBoolean  zero_Head_Flag;
	CTFBoolean  artifact_mode;
} __PACKED__ new_general_setup_rec_ext;


// meg41GeneralResRec
struct meg41GeneralResRec {
	char_t						appName[256];
	char_t						dataOrigin[256];
	char_t						dataDescription[256];
	int16_t						no_trials_avgd;
	char_t						data_time[255];
	char_t						data_date[255];
	new_general_setup_rec_ext	gSetUp;
	meg4FileSetup				nfSetUp;
} __PACKED__;       // so we dont have alignment problems
		    // padding is by explicit declarations

typedef struct meg41GeneralResRec meg41GeneralResRec;


// SensorCoefResRec
typedef struct {
	char_t		sensorName[32];
	uint32_t	coefType;
	int32_t		reserved1;		// pad out to the next 8 byte boundary
	CoefResRec	coefRec;
} SensorCoefResRec;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _H_MegDefs */
