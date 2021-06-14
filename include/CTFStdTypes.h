/* ==========================================================================
 * $Header: /data/dougmck/code/ctf/ctfutil/include/RCS/CTFStdTypes.h,v 13.0 1999/12/09 19:35:29 klam Exp $
 * ==========================================================================
 *
 * Author's comments
 *
 * Copyright (c) CTF Systems Inc., 1995-1996.  All Rights Reserved.
 * Commercially Confidential Information
 *
 * ==========================================================================
 *
 * $Log: CTFStdTypes.h,v $
 * Revision 13.0  1999/12/09 19:35:29  klam
 * Revised to v13.0.
 *
 * Revision 6.0  1999/11/03 23:01:06  dnm
 * Initial check in for Linux/HP-UX 10 Code
 *
** Revision 10.4  97/04/30  14:39:19  14:39:19  hu (Hu)
** *** empty log message ***
** 
** Revision 10.3  97/04/18  11:39:27  11:39:27  murray (Murray)
** added commercially confidential comment
** 
** Revision 10.2  96/12/12  15:10:07  15:10:07  sixtus (Sixtus Lee)
** Added Boolean definition.
** 
** Revision 10.1  96/12/03  15:34:57  15:34:57  sixtus (Sixtus Lee)
** Redefine CTFBoolean as typedef int (instead of enum).
** Define true and false as 1 and 0, respectively.
** 
** Revision 10.0  96/01/18  16:53:32  16:53:32  qabbany (Moustafa Elqabbany)
** Major Revision!
** 
 * Revision 1.2  95/11/30  16:55:49  16:55:49  tmetzger (Tom)
 * Added
 * 
 * Revision 1.1  95/11/30  16:52:17  16:52:17  tmetzger (Tom)
 * Initial revision
 * 
 *
 * ==========================================================================
 */
/******************************************************************************
// CTFStdTypes.h
//
//	
//
//	Definition file for portable type declarations
//
//	Copyright й 1994 CTF Systems Inc. 1994. All rights reserved.
//
//	Version 1.0		June 28, 1994	
//	Original release
//
//
//	1.1			June 29, 1994	
//	Added define for C or C++ 
//
//
//	1.2			August 2, 1994	
//	Added #defines for xDBL_MAX, xDBL_MIN, xDBL_EPSILON etc. for double 
//	types as well as SDOUBLE_SIZE, LDOUBLE_SIZE and SINGLE_SIZE
//
//
//	1.3			October 12, 1994
//	Added definition of CStrPtrC, which is a pointer to a constant C
//	string.
//
//	1.4			November 9, 1994
//	Made the inclusion of <types.h> conditional.
//
//	Standard C Includes
*/

#ifndef H_CTFStdTypes
#define	H_CTFStdTypes

#include <stdint.h>
#include <math.h> 
#include <limits.h>

/*********************************
*				 *
*  Floating Point Definitions	 *
*				 *
*********************************/
/*
//	Three definitions are provided for floating point operation
//
//	DDouble - used for general floating point operation. Specifies
//	the most efficient size for floats for the particular
//	platform. This type can be any size and is determined 
//	by the compiler.
//
//	SDouble - used for importing/exporting data where data has to be
//	a standard size and requires double precision. This
//	should represent a standard IEEE double precision floating
//	point number.
//
//	Single - used for importing/exporting data where data has to be
//	a standard size and requires single precision.  This
//	should represent a standard IEEE single precision floating
//	point number.
//
//	Note: I have not used the the size in the type names. DDouble
//	length can vary from compiler to compiler. Single and SDouble
//	imply the level of precision, not the size of the varible. 
//		
//
//	for 680x0
//
//			double 		-> 12 bytes
//			short double 	-> 8 bytes
//
//	for PowerPC
//
//			double		-> 8 bytes
//			long double	-> 8 bytes
*/

typedef	double			SDouble;
typedef	double			DDouble;

#define DDBL_MANT_DIG		DBL_MANT_DIG
#define DDBL_DIG			DBL_DIG
#define DDBL_MIN_EXP		DBL_MIN_EXP
#define DDBL_MIN_10_EXP		DBL_MIN_10_EXP
#define DDBL_MAX_EXP		DBL_MAX_EXP
#define DDBL_MAX_10_EXP		DBL_MAX_10_EXP
#define DDBL_MAX			DBL_MAX
#define DDBL_EPSILON		DBL_EPSILON
#define DDBL_MIN			DBL_MIN
/*
//	DNM August 2, 1994
//
//	Definitions for short doubles.
//	
//	These values were copied from <float.h>. Except for 
//	"(*...) __D_xxx" values, these values should be portable because
//	they are for an IEEE 64 bit double.
*/
#define SDBL_MANT_DIG			53
#define SDBL_DIG				15
#define SDBL_MIN_EXP			(-1021)
#define SDBL_MIN_10_EXP			(307)
#define SDBL_MAX_EXP			1024
#define SDBL_MAX_10_EXP			308
#define SDBL_MAX				(* (double *) __D_MAX)
#define SDBL_EPSILON			(* (double *) __D_EPSILON)
#define SDBL_MIN				(* (double *) __D_MIN)

typedef		float	Single;

/*
//	Use the following in function argument lists when arguments are an array
//	of floats
*/

typedef		DDouble*	DDoubleArr;	/* Pointer to aaray of floats */
typedef		SDouble*	SDoubleArr;	/* Pointer to aaray of floats */
typedef		Single*		SingleArr;	/* Pointer to aaray of floats */

enum	{
	SDOUBLE_SIZE = sizeof( SDouble ),
	LDOUBLE_SIZE = sizeof( DDouble ),
	SINGLE_SIZE = sizeof( Single ) };


/************************
*						*
*   C String defines	*
*						*
*************************/
/*
//	Note:	Definitions are prefixed with a 'C' to distinguish strings from
//		Pascal strings used by the Macintosh Toolbox
//
//	The define for 'char_t' should be changed to 'signed' or 'unsigned' as 
//	appropriate for the compiler.
//
*/


/* could be 'unsigned char' in some compilers */
typedef	char	char_t;			

/* 255 character string + NULL terminator */
typedef	char_t		CStr256[256];	

/* 79 character string + NULL terminator */
typedef	char_t		CStr80[80];			

/* 31 character string + NULL terminator */
typedef	char_t		CStr32[32];			

/* 15 character string + NULL terminator */
typedef char_t		CStr16[16];

/* pointer to string */
typedef	char_t*		CStrPtr;			

/* pointer to constant string */
typedef const char_t*	CStrPtrC;


/* pointer to list of strings	ее DNM v1.1 */
typedef	CStrPtr*	CStrList;			

/* single character */
typedef	char_t		CChar;				
typedef unsigned char	UCChar;

/* pointer to constant unsigned string */
typedef const UCChar*	UCStrPtrC;

/************************
*						*
* Integer defines		*
*						*
*************************/
/*
//	Note:	Definitions prefixed with a 'U' are unsigned
//
//	The size of the standard 'int' is compiler dependent and shoud be 
//	avoided.
//
//	IntxxPtr and UIntxxPtr's are to be used to pass arrays of the respective
//	types where no parameters are to be returned
*/

typedef	char	Int8;			/* -128 to +127 */
typedef	short	Int16;			/* -32768 to 32767 */
/* this works on the hp's and under Think C since ints are 4 bytes long */
typedef	long	Int32;			/* -2,147,483,648 to 2,147,483,647 */

typedef	int8_t*		Int8Arr;		/* Pointers to Intxx */
typedef	int16_t*	Int16Arr;				
typedef	int32_t*	Int32Arr;					


// typedef	unsigned char	UInt8;		/* 0 to 255 */
// typedef	unsigned short	UInt16;		/* 0 to 65535 */
/* this works on the hp's and under Think C since ints are 4 bytes long */
// typedef	unsigned long	UInt32;		/* 0 to 4,294,967,295 */

typedef	uint8_t*	UInt8Arr;		/* Pointers to UIntxx */
typedef	uint16_t*	UInt16Arr;			
typedef	uint32_t*	UInt32Arr;			

/*
//	Limit values of the integer types
*/

#define	Int8_MAX	CHAR_MAX	
#define	Int16_MAX	SHRT_MAX	
#define	Int32_MAX	LONG_MAX	

#define	Int8_MIN	CHAR_MIN	
#define	Int16_MIN	SHRT_MIN	
#define	Int32_MIN	LONG_MIN	

#define	UInt8_MAX	UCHAR_MAX	
#define	UInt16_MAX	USHRT_MAX	
#define	UInt32_MAX	ULONG_MAX	



/**********************
*		      *
* Byte defines        *
*		      *
**********************/

typedef uint8_t		Bit8;
typedef uint16_t	Bit16;
typedef uint32_t	Bit32;


/******************************************************************************
*
*
*
*
*	Additional Standard types
*
*
*
*******************************************************************************/
/*
//	CTFBoolean	is provided for values that can	have either True
//	or False values.
*/

typedef int CTFBoolean;

// TEST TEST TEST
//
// The type bool has values of either true or false. The value of true equals 1
// and the value of false equals 0. Therefore, don't need to redefine false/true
// again.
//
//#define false 0
//#define true 1

#ifndef _XtIntrinsic_h
// typedef char    Boolean;
#define FALSE 0
#define TRUE 1
#endif

/* 
//	OSErr is provided as a return type for functions that return an error
//	or status code. 
*/
#if 0			/* OSErr is defined for Think C */
typedef	short	OSErr;	/* Type for error codes returned from functions	 */
#endif



/******************
*		  *
*   End of File	  *
*		  *
******************/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif	/*	H_CTFStdTypes */

