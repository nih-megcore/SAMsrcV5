// SAMfiles.h -- structures for file headers for SAM
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//
// SAM coefficient files are structured as follows:
//
//	char		Identity[8] = "SAMCOEFF";	// uniquely identifies coefficient file
//	SAM_HDR		SAMHeader;			// SAM header (what else could it be?)
//	int		ChannelIndex[M];		// index of used primary sensor channel numbers
//	double		SAMCoeffs[0][M];		// 1st SAM coefficient set (units = A-m/T)
//	double		SAMCoeffs[1][M];		// 2nd SAM coefficient set
//			"
//			"
//	double		SAMCoeffs[V][M];		// last SAM coefficient set
//
// SAM version 3 coefficient files are structured as follows:
//
//	char		Identity[8] = "SAMCOEFF";	// uniquely identifies coefficient file
//	SAM_HDR		SAMHeader;			// SAM header (what else could it be?)
//	int		ChannelIndex[M];		// index of used primary sensor channel numbers
//
//		if weights are derived from a list of a few coordinates (100 or fewer):
//		labels & coordinates:
//	TARGET		Voxel[0]
//	TARGET		Voxel[1]
//			"
//			"
//	TARGET		Voxel[V]
//
//	double		SAMCoeffs[0][M];		// 1st SAM coefficient set (units = A-m/T)
//	double		SAMCoeffs[1][M];		// 2nd SAM coefficient set
//			"
//			"
//	double		SAMCoeffs[V][M];		// last SAM coefficient set
//
//
// SAM static image files are structured as follows:
//
//	char		Identity[8] = "SAMIMAGE";	// uniquely identifies image file
//	SAM_HDR		SAMHeader;			// SAM header
//	int		ChannelIndex[M];		// index of used primary sensor channel numbers
//	double		Voxel[0];			// 1st SAM voxel (units = A-m, (A-m)^2, Z, T, F, or P)
//	double		Voxel[1];			// 2nd SAM voxel
//			"
//			"
//	double		Voxel[V];			// last SAM voxel
//
//	Coefficients & image voxels are ordered in X,Y,Z sequence, with Z the least
//		significant index (most rapidly changing), Y is next, and then X.
//		Coordinate indices always advance in the positive direction. This implies
//		that Voxel[0] is in the right, posterior, inferior position relative to
//		the region of interest (bounding box of image).
//
// SAM covariance files are structered as follows:
//
//	char		Identity[8] = "SAMCOVAR";	// uniquely identifies covariance file
//	COV_HDR		CovHeader;			// SAM covariance header
//	int		ChannelIndex[M];		// index of used primary sensor channels
//	double		Cov[0][0];			// 1st covariance element
//	double		Cov[0][1];			// 2nd covariance element
//				"
//				"
//	double		Cov[0][M];
//	double		Cov[1][0];
//	double		Cov[1][1];
//				"
//				"
//	double		Cov[M][M];			// last covariance element
//

#ifndef H_SAMFILES
#define H_SAMFILES

#include <stdint.h>

// covariance file versions
#define COV_VERSION	1					// this is version 1 -- got it?
#define SAM_VERSION	1					// this, too!
#define SAM_VERSION_2		2
#define SAM_VERSION_3		3			// I like this one!

// covariance type definitions
#define	GLOBAL_				0			// global covariance
#define SUM_				1			// sum of covariance for all markers
#define ORIENT_				2			// orientation covariance (global)
#define NOISE_				3			// noise covariance (of global)
#define MARKER1_			4			// marker 1
#define MARKER2_			5			// marker 2
#define MARKER3_			6			// marker 3
#define MARKER4_			7			// marker 4
#define MARKER5_			8			// marker 5
#define MARKER6_			9			// marker 6
#define ALL_				10			// separate covariance for each marker

// SAM file types
#define SAM_TYPE_IMAGE		0			// flags file as a SAM static image file
#define SAM_TYPE_WT_ARRAY	1			// flags file as SAM coefficients for regular target array
#define SAM_TYPE_WT_LIST	2			// flags file as SAM coefficients for target list

// define SAM unit types
#define	SAM_UNIT_COEFF	0				// SAM coefficients A-m/T
#define	SAM_UNIT_MOMENT	1				// SAM source (or noise) strength A-m
#define	SAM_UNIT_POWER	2				// SAM source (or noise) power (A-m)^2
#define SAM_UNIT_RVE	3				// SAM source entropy
#define	SAM_UNIT_SPMZ	4				// SAM z-deviate
#define	SAM_UNIT_SPMF	5				// SAM F-statistic
#define	SAM_UNIT_SPMT	6				// SAM T-statistic
#define	SAM_UNIT_SPMP	7				// SAM probability
#define SAM_UNIT_MUSIC	8				// MUSIC metric
#define SAM_UNIT_G2		9				// SAM kurtosis
#define SAM_UNIT_NORM	10				// SAM normalized coefficients

// 'COV_HDR' -- since this is written & read from disk, we must make the integer types portable
typedef struct {
	int32_t		Version;			// file version number
	char		SetName[256];		// name of parent dataset
	char		SpecName[256];		// name of covariance specification file
	int32_t		NumChans;			// number of channels used by SAM
	double		HPFreq;				// highpass frequency (Hz)
	double		LPFreq;				// lowpass frequency (Hz)
	double		BWFreq;				// bandwidth of filters (Hz)
	double		Noise;				// estimated covariance noise power
	int32_t		NumSegments;		// number of time-segments
	int32_t		NumSamples;			// toal number of samples
	int32_t		CovType;			// covariance type
	uint32_t	dummy1;				// alignment
} COV_HDR;

// 'SAM_HDR' is to be used for both SAM coefficients & SAM static images -- since this is written & read from disk, we must make the integer types portable
typedef struct {
	int32_t		Version;			// 32				file version number (1)
	char		SetName[256];		// 288 (32+256)		name of parent dataset
	int32_t		NumChans;			// 292 (288+4)		number of channels used by SAM
	int32_t		NumWeights;			// 296 (292+4)		number of SAM virtual sensors (0=static image)
	uint32_t	reserved1;			// 300 (296+4)		padding
	double		XStart;				// 308 (300+8)		x-start coordinate (m)
	double		XEnd;				// 316 (308+8)		x-end coordinate (m)
	double		YStart;				// 324 (316+8)		y-start coordinate (m)
	double		YEnd;				// 332 (324+8)		y-end coordinate (m)
	double		ZStart;				// 340 (332+8)		z-start coordinate (m)
	double		ZEnd;				// 348 (340+8)		z-end coordinate (m)
	double		StepSize;			// 356 (348+8)		voxel step size (m)
	double		HPFreq;				// 364 (356+8)		highpass frequency (Hz)
	double		LPFreq;				// 372 (364+8)		lowpass frequency (Hz)
	double		BWFreq;				// 380 (372+8)		bandwidth of filters (Hz)
	double		MeanNoise;			// 388 (380+8)		mean primary sensor noise (T)
	char		MriName[256];		// 644 (388+256)	MRI image file name
	int32_t		Nasion[3];			// 656 (644+12)		MRI voxel index for nasion
	int32_t		RightPA[3];			// 668 (656+12)		MRI voxel index for right pre-auricular
	int32_t		LeftPA[3];			// 680 (668+12)		MRI voxel index for left pre-auricular
	int32_t		SAMType;			// 684 (680+4)		SAM file type or total number of trial
	int32_t		SAMUnit;			// 688 (684+4)		SAM units (a bit redundant, but may be useful) or number of active-state trials
	uint32_t	reserved2;			// 672 (688+4)		more padding
} SAM_HDR;

// version 2 SAM header is not what I would like! -- since this is written & read from disk, we must make the integer types portable
typedef struct {
	int32_t		Version;			// file version number (2)
	char		SetName[256];		// name of parent dataset
	int32_t		NumChans;			// number of channels used by SAM
	int32_t		NumWeights;			// number of SAM virtual sensors (0=static image)
	uint32_t	reserved1;			// padding
	double		XStart;				// x-start coordinate (m)
	double		XEnd;				// x-end coordinate (m)
	double		YStart;				// y-start coordinate (m)
	double		YEnd;				// y-end coordinate (m)
	double		ZStart;				// z-start coordinate (m)
	double		ZEnd;				// z-end coordinate (m)
	double		StepSize;			// voxel step size (m)
	double		HPFreq;				// highpass frequency (Hz)
	double		LPFreq;				// lowpass frequency (Hz)
	double		BWFreq;				// bandwidth of filters (Hz)
	double		MeanNoise;			// mean primary sensor noise (T)
	char		MriName[256];		// MRI image file name
	int32_t		Nasion[3];			// MRI voxel index for nasion
	int32_t		RightPA[3];			// MRI voxel index for right pre-auricular
	int32_t		LeftPA[3];			// MRI voxel index for left pre-auricular
	int32_t		SAMType;			// SAM file type
	int32_t		SAMUnit;			// SAM units (a bit redundant, but may be useful)
	uint32_t	reserved2;			// more padding
	double		MegNasion[3];		// MEG dewar coordinates for nasion (m)
	double		MegRightPA[3];		// MEG dewar coordinates for right pre-auricular (m)
	double		MegLeftPA[3];		// MEG dewar coordinates for left pre-auricular (m)
	char		SAMUnitName[32];	// SAM unit name
} SAM_HDR_v2;

// 'SAM_HDR' is to be used for both SAM coefficients & SAM static images -- since this is written & read from disk, we must make the integer types portable
typedef struct {
	int32_t		Version;			// file version number (3)

	// specification of input dataset & output dimensions
	char		SetName[256];		// name of parent dataset
	int32_t		NumChans;			// number of channels used by SAM
	int32_t		NumWeights;			// number of SAM virtual sensors (0=static image)
	uint32_t	reserved1;			// padding

	// specification of ROI, if either 'SAM_TYPE_IMAGE' or 'SAM_TYPE_WT_ARRAY'
	double		XStart;				// x-start coordinate (m)
	double		XEnd;				// x-end coordinate (m)
	double		YStart;				// y-start coordinate (m)
	double		YEnd;				// y-end coordinate (m)
	double		ZStart;				// z-start coordinate (m)
	double		ZEnd;				// z-end coordinate (m)
	double		StepSize;			// voxel step size (m)

	// specification of time window relative to designated marker
	double		TimeStart;			// start of time window (seconds) relative to marker (or start of dataset)
	double		TimeEnd;			// end of time window (seconds) relative to marker (or end of dataset)
	char		MarkerName[128];	// name of marker, or 'Global' if no markers
	
	// specification of frequency range & noise bandwidth of weights or image
	double		HPFreq;				// highpass frequency (Hz)
	double		LPFreq;				// lowpass frequency (Hz)
	double		BWFreq;				// bandwidth of filters (Hz)
	double		MeanNoise;			// mean primary sensor noise (T)

	// specification of MRI & voxel indices of fiducial markers
	char		MriName[256];			// MRI image file name
	int32_t		Nasion[3];			// MRI voxel index for nasion
	int32_t		RightPA[3];			// MRI voxel index for right pre-auricular
	int32_t		LeftPA[3];			// MRI voxel index for left pre-auricular

	// specification of SAM weight or image type & values
	int32_t		SAMType;			// SAM file type
	int32_t		SAMUnit;			// SAM units (a bit redundant, but may be useful)
	uint32_t	reserved2;			// more padding
} SAM_HDR_V3;

// 'TARGET' holds name & location of each voxel in list
typedef struct {
	double		target_pos[3];			// cartesian position
	char		target_label[16];		// name assigned to target
} TARGET;

#endif	// H_SAMFILES
