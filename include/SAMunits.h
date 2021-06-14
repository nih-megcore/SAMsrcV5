/* ==========================================================================
 * $Header: SAMunits.h,v 1.1 99/04/06 14:47:54 ser Exp $
 * ==========================================================================
 *
 * Author's comments
 *
 * SAMunits -- external string definitions of SAM file units
 *
 * Author: Stephen E. Robinson
 * Copyright (c) CTF Systems Inc., 1995-1996.  All Rights Reserved.
 * Commercially Confidential Information
 *
 * ==========================================================================
 *
 * $Log:	SAMunits.h,v $
 * Revision 1.1  99/04/06  14:47:54  14:47:54  ser (Stephen Robinson)
 * Initial revision
 * 
 *
 * ==========================================================================
 */

#ifndef H_SAMUNITS
#define H_SAMUNITS

// SAM units -- these are not used in SAM analysis, but may be filled in for (or by) display programs
extern char	*SAMUnitName[7] = {		// SAM unit names:
	"A-m/T",				// 	SAM coefficients
	"A-m",					//	SAM source (or noise) strength
	"(A-m)^2",				//	SAM source (or noise) power
	"Z",					//	SAM z-deviate
	"F",					//	SAM F-statistic
	"T",					//	SAM T-statistic
	"P"					//	SAM probability
};

#endif	// H_SAMUNITS
