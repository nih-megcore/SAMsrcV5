/* TimeSegs.h -- structures for SAM covariance integration
 *
 * Author: Stephen E. Robinson
 *
 * Copyright (c) CTF Systems Inc., 1995-2009.  All Rights Reserved.
 * Commercially Confidential Information
 */

#ifndef H_TIMESEGS
#define H_TIMESEGS

#define MAX_MARKERS	100			// number of different segment types

#include <SAMfiles.h>

// minimalistic trial structure
typedef struct {
	int		Epoch;				// epoch
	int		TS;					// trial starting sample
	int		TE;					// trial ending sample
} TRIAL;

// time-ordered list of markers in dataset
typedef struct {
	int		Epoch;				// epoch
	double	Latency;			// time latency to trigger
	char	Name[128];			// marker name
} SAM_MARKS;

// this structure is the basis for the covariance integration list
typedef struct {
	int		Epoch;				// epoch number
	int		TS;					// starting sample
	int		TE;					// ending sample (inclusive)
	int		TT;					// segment sample count
} COV_SEG;

// covariance specs for an active or control state
typedef struct {
	char		Name[128];		// marker name/label
	double		TimeStart;		// start-time -- relative to marker
	double		TimeEnd;		// end-time -- relative to marker
	int			State;			// state flag -- active, control, etc.
} COV_SPEC;

#endif	// H_TIMESEGS
