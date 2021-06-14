/* ==========================================================================
* $Header: /data/dougmck/code/ctf/ctfutil/include/RCS/MEGTypes.h,v 13.0 1999/12/09 19:35:29 klam Exp $
* ==========================================================================
*
* Author: Sixtus Lee
*
* Copyright (c) CTF Systems Inc., 1995-1996.  All Rights Reserved.
* Commercially Confidential Information
*
* ==========================================================================
*
* $Log: MEGTypes.h,v $
* Revision 13.0  1999/12/09 19:35:29  klam
* Revised to v13.0.
*
* Revision 6.0  1999/11/03 23:01:14  dnm
* Initial check in for Linux/HP-UX 10 Code
*
 * Revision 10.3  97/06/04  17:15:55  17:15:55  johnj (John)
 * added invalid value stuff
 *
 * Revision 10.2  97/04/25  16:08:29  16:08:29  murray (Murray)
 * added Trigger_t, Marker_t, and TrialClass_t
 *
 * Revision 10.1  97/04/18  12:03:07  12:03:07  murray (Murray)
 * added commercially confidential comment
 *
 * Revision 10.0  96/01/18  16:53:45  16:53:45  qabbany (Moustafa Elqabbany)
 * Major Revision!
 *
 * Revision 1.2  96/01/05  14:51:49  14:51:49  sixtus (Sixtus Lee)
 * No changes.
 *
 * Revision 1.1  95/11/10  19:45:45  19:45:45  sixtus (Sixtus Lee)
 * Initial revision
 *
*
* ==========================================================================*/
/******************************************************************************
// MEGTypes.h
//
//
//
//      Definition file for MEG standard types.
//
//      Copyright 1994 CTF Systems Inc. 1994. All rights reserved.
//
//      Version 1.0             October 12, 1994
//      Original release
//
*/

#ifndef H_MEGTypes
#define H_MEGTypes

#include <limits.h>
#include <CTFStdTypes.h>
#include <CoeffTypes.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

//
// Standard MEG Types
// Sample_t should be the same type as size_t
//
typedef uint32_t        Sample_t;               // Number of samples type
typedef int16_t         Channel_t;              // Number of channels type
typedef int16_t         Trial_t;                // Number or trials type
typedef SDouble         CTFTime;                // Time type
typedef SDouble         SampleRate;             // Sample Rate type
typedef int16_t         Run_t;                  // Run Number type
typedef int16_t         Coil_t;                 // Coil Number type
typedef int16_t         Trigger_t;      // Trigger number type
typedef int16_t         Marker_t;       // Number of markers
typedef int16_t         TrialClass_t;   // Number of trial classifications

//
// Invalid values for standard MEG types
//
#if 0 // @@@ ULONG vs uint32 
const Sample_t			INVALID_SAMPLE_NUMBER		= ULONG_MAX;
#endif
const Channel_t			INVALID_CHANNEL_NUMBER		= SHRT_MAX;
const Trial_t			INVALID_TRIAL_NUMBER		= SHRT_MAX;
const Run_t				INVALID_RUN_NUMBER			= SHRT_MAX;
const Coil_t			INVALID_COIL_NUMBER			= SHRT_MAX;
const Trigger_t			INVALID_TRIGGER_NUMBER		= SHRT_MAX;
const Marker_t			INVALID_MARKER_NUMBER		= SHRT_MAX;
const TrialClass_t		INVALID_TRIALCLASS_NUMBER	= SHRT_MAX;

/* Formats for MEG Data */

enum FormatType
{
	fMEG4,  fMEG3
};

enum TriggerMode
{
	noTrigger,              TrialPerTrigger,        BeginTrialTrigger,
	BeginTrialKeyboard
};

/******************
*                 *
*   End of File   *
*                 *
******************/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /*      H_MEGTypes */

