/* ==========================================================================
 * $Header: CTF_MRI_Headers.h,v 1.2 98/06/06 16:07:56 dougc Exp $
 * ==========================================================================
 *
 * Author's comments
 *
 * Copyright (c) CTF Systems Inc., 1995-1996.  All Rights Reserved.
 * Commercially Confidential Information
 *
 * ==========================================================================
 *
 * $Log:	CTF_MRI_Headers.h,v $
 * Revision 1.2  98/06/06  16:07:56  16:07:56  dougc (Doug)
 * updated MRI file header to version 2.2
 * includes transformation matrix and original slice thickness
 * 
 * Revision 1.1  97/05/20  11:04:19  11:04:19  johnj (John)
 * Initial revision
 * 
 *
 * ==========================================================================
 */
#ifndef	_CTFMRIDEFS
#define _CTFMRIDEFS
// ************************
// 	file CTF_MRI_Headers.h
//
//	new definitions for MRI (CTF) file formats
//
//	D. Cheyne, July, 1996
// 
// ************************

// slice direction defines -- numeric order is arbitrary!
// note that the following relationship is used in 
// coordinate transfomrations in MRIViewer 
//  sagittal = x
//	coronal = y
// 	axial = z
// with the origin at the upper left, anterior corner of the volume, 
// producing a right-handed coordinate system) with
// increaseing values in the x, y, and z axes for typical slice storage 
// formats i.e., starting at upper left corner and scanning down row by row.
//
enum {  coronal = 1,
		sagittal,
		axial
};

// new index for inputOrientation field 
// used to signify slice orientation ....
#define LEFT_ON_RIGHT 	1
#define LEFT_ON_LEFT 	0

// image type ...
enum { Modality_MRI = 0,
		Modality_CT,
		Modality_PET,
		Modality_SPECT,
		Modality_OTHER
};		

#define VERSION_1_STR "CTF_MRI_FORMAT VER 1.0"
#define VERSION_2_STR "CTF_MRI_FORMAT VER 2.0"
#define VERSION_21_STR "CTF_MRI_FORMAT VER 2.1"
#define VERSION_22_STR "CTF_MRI_FORMAT VER 2.2"


// **************************************************************************
//	this structure used to read MRI files and get some general 
// 	information in the headers used by Import routine etc...
// **************************************************************************
typedef struct MRI_File_Params {
	int		fileType;
	int		version;
	char	filePath[256];
	char	listFileName[256];
	long	fileOffset;
	int		useEndOfFile;
	int		imagesPerFile;
	int		numFiles;
	int		swapBytes;
	int		dataSize;
	int		imageSizeRows;
	int		imageSizeCols;
	int		paddingBytes;
	int		packedBytes;
	long	maxFileValue;
	int		clippingValue;
	int		imageOrientation;
	int		flipHorizontal;
	int		flipVertical;
	int		invertSlices;
	double	sliceFOVHigh;
	double	sliceFOVWide;
	double	sliceSpacing;
	double	sliceThickness;
	double	sliceGap;
	int		fitMethod;
	int		scaleFactor;

} MRI_File_Params;


// **************************************************************************
// *** Version 1.0 header -- used for Macintosh program compatibility only
// *** 256 bytes in length ***
// **************************************************************************
typedef struct MRI_Header_Struct {
	char	filetype[32];			
	char	patient_info[32];		// unused
	long	imageResolution;		// 4 bytes
	long	numSlices;				// 4 bytes
	long	inputOrientation;		// 0 = see above 
	float	mmPerPixel_sagital;		// 4 byte IEEE single precision
	float	mmPerPixel_coronal;
	float	mmPerPixel_axial;			
	float	scale;					// unused		
	long	Nasion_Sag;				// pixel coordinates of fiduciary points
	long	Nasion_Cor;
	long	Nasion_Axi;
	long	LeftEar_Sag;
	long	LeftEar_Cor;
	long	LeftEar_Axi;
	long	RightEar_Sag;
	long	RightEar_Cor;
	long	RightEar_Axi;									
	float	SphereParams[4];		// sphere coordinates, x, y, z, Radius (in mm)		
	float	NotUsed[12];				// save for future use	
	float	TranslationMatrix[16];	// translation matrix Head to MRI-- not used ?
} MRI_Header_Struct;



// **************************************************************************
// ********** additional header structs for  CTF MRI version 2.0 file 
// **************************************************************************

typedef struct HeadModel_Info {
	short			Nasion_Sag;				// fiduciary point voxel locations
	short			Nasion_Cor;				// Sag = sagittal direction
	short			Nasion_Axi;				// Cor = coronal direction
	short			LeftEar_Sag;			// Axi = axial direction
	short			LeftEar_Cor;
	short			LeftEar_Axi;
	short			RightEar_Sag;
	short			RightEar_Cor;
	short			RightEar_Axi;									
	float			defaultSphereX;	 		// default sphere parameters in mm
	float			defaultSphereY;			// (in head based coordinate system
	float			defaultSphereZ;		
	float			defaultSphereRadius;	
} HeadModel_Info;

typedef struct Image_Info {					// scan and/or sequence parameters 
	short			modality;				// 0=MRI, 1=CT, 2 = PET, 3=SPECT, 4=OTHER
	char			manufacturerName[64];
	char			instituteName[64];
	char			patientID[32];
	char			dateAndTime[32];
	char			scanType[32];
	char			contrastAgent[32];
	char			imagedNucleus[32];
	float			Frequency;
	float			FieldStrength;
	float			EchoTime;
	float			RepetitionTime;
	float			InversionTime;
	float			FlipAngle;
	short			NoExcitations;
	short			NoAcquisitions;
	char			commentString[256];
	char			forFutureUse[64];
} Image_Info;

typedef struct Version_2_Header {
	char			identifierString[32];		// "CTF_MRI_FORMAT VER 2.x"	
	short			imageSize;					// always = 256
	short			dataSize;					// 1 = 8 bit data, 2 = 16 bit data
	short			clippingRange;				// max. integer value in data
	short			imageOrientation;			// 0 = left on left, 1 = left on right 		
	float			mmPerPixel_sagittal;		// voxel dimensions in mm
	float			mmPerPixel_coronal;			// voxel dimensions in mm
	float			mmPerPixel_axial;			// voxel dimensions in mm
	HeadModel_Info	headModel;					// structure defined above (34 bytes)
	Image_Info		imageInfo;					// structure defined above (638 bytes)
	float			headOrigin_sagittal;		// voxel location of head origin
	float			headOrigin_coronal;			// voxel location of head origin
	float			headOrigin_axial;			// voxel location of head origin
						// ** euler angles to align MR to head coordinate system (angles in degrees!)
	float			rotate_coronal;				// 1. rotate in coronal plane by this angle
	float			rotate_sagittal;			// 2. rotate in sagittal plane by this angle
	float			rotate_axial;				// 3. rotate in axial plane by this angle 
	short			orthogonalFlag;				// true if image is orthogonalized to head frame
						// the following is present only in version 2.2 and greater
	short			interpolatedFlag;			// true if slices were interpolated during conversion
	float			originalSliceThickness;		// 
	float			transformMatrix[4][4];		// 4x4 transformation matrix (head to mri)
	unsigned char	unused[202];				// pad header to 1028 bytes
} Version_2_Header;

#endif
