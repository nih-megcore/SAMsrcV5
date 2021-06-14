/* ==========================================================================
 * $Header: /data/dougmck/code/ctf/ctfutil/include/RCS/CoeffTypes.h,v 13.1 2000/06/20 15:15:07 dougmck Exp $
 * ==========================================================================
 *
 * Author's comments
 *
 * Copyright (c) CTF Systems Inc., 1995-1996.  All Rights Reserved.
 * Commercially Confidnetial Information
 *
 * ==========================================================================
 *
 * $Log: CoeffTypes.h,v $
 * Revision 13.1  2000/06/20 15:15:07  dougmck
 * Merged in changes for release 4.12.2 ( HSC )
 *
 * Revision 10.4  20/0./1.  8.:6.:4.  8.:6.:4.  milee ()
 * GxAB and GxAI added to do adaptive balancing.
 * Revision 13.0  1999/12/09 19:35:29  klam
 * Revised to v13.0.
 * 
 * Revision 10.3  10/.0/.1  .1:.0:.4  .1:.0:.4  milee ()
 * Add G1AR, G2AR and G3AR
 * Revision 6.0  1999/11/03 23:01:09  dnm
 * Initial check in for Linux/HP-UX 10 Code
 * 
 * Revision 10.2  97/04/18  11:42:57  11:42:57  murray (Murray)
 * added commercially confidential comment
 * 
 * Revision 10.1  96/11/14  14:28:26  14:28:26  johnj (John)
 * added G1OI
 * 
 * Revision 10.0  96/01/18  16:53:35  16:53:35  qabbany (Moustafa Elqabbany)
 * Major Revision!
 * 
 * Revision 1.2  95/11/30  16:55:51  16:55:51  tmetzger (Tom)
 * Added
 * 
 * Revision 1.1  95/11/30  16:52:17  16:52:17  tmetzger (Tom)
 * Initial revision
 * 
 *
 * ==========================================================================
 */
#ifndef H_COEFFTYPES
#define H_COEFFTYPES

#include <CTFStdTypes.h>


const	Bit32		NOGRAD	= 0x00000000;
const	Bit32		G1BR	= 0x47314252;
const	Bit32		G2BR	= 0x47324252;
const	Bit32		G3BR	= 0x47334252;
const	Bit32		G2OI	= 0x47324f49;
const	Bit32		G3OI	= 0x47334f49;
const	Bit32		G1OI	= 0x47314f49;
const	Bit32		G0AR	= 0x47304152;
const   Bit32		G1AR    = 0x47314152;
const	Bit32		G2AR	= 0x47324152;
const   Bit32       G3AR    = 0x47334152;
const 	Bit32		G0AB	= 0x47304142;
const 	Bit32		G1AB	= 0x47314142;
const	Bit32		G2AB	= 0x47324142;
const	Bit32		G3AB	= 0x47334142;
const	Bit32		G1AI	= 0x47314149;
const	Bit32		G2AI	= 0x47324149;
const	Bit32		G3AI	= 0x47334149;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/* Coefficient gradient orders */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* H_COEFFTYPES */

