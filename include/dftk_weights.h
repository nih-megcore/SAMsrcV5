#ifndef DFTK_WEIGHTS_H
#define DFTK_WEIGHTS_H
/*******************************************************************************

	$Id: dftk_weights.h,v 1.7 2000/01/11 22:48:55 skip Exp $

	Name:		(file) dftk_weights.h

	Purpose:	Defines the classes, structures and user block names for
			the system configuration file's weights.

	Notes for searching for E tables & weight tables in system config files:

			This applies to system config files only.

			The intended use is that there will be several
			dftk_E_table structures in the system config file.
			Each will go by the name, DFTK_E_TABLE_NAME, and
			can be accessed by dftk_config::search_user_block(...)
			and its related functions.  In the E table structure
			is a filter name, filter, which is used to look up
			the related weight tables (concatenate DFTK_WT_TABLE_PREFIX
			and the filter name).

			To find all weight tables, one can also go through all
			user blocks and check the first few characters of the
			name with DFTK_WT_TABLE_PREFIX.  i.e.:

			strncmp(block_name, DFTK_WT_TABLE_PREFIX, strlen(DFTK_WT_TABLE_PREFIX)


	Notes for obtaining E & weight tables from run config files:

			There will be one and only one of each, and these will go
			by the names APPLIED_E_TABLE_NAME and APPLIED_WT_TABLE_NAME.


	Notes for accessing E table entries:

			The "entries" field is the first of an array of total_entries.
			They are indexed by CFO (config-file order / channel number),
			so use the channel number.

	Notes for creating E tables:

			E tables should have the version set to
			CURRENT_DFTK_E_TABLE_VERSION

	Author:		slee

	Creation Date:	12-05-95

	Copyright 1995, Biomagnetic Technologies, inc.

	$Log: dftk_weights.h,v $
	Revision 1.7  2000/01/11 22:48:55  skip
	Extensive changes to support dynamic weight tables.

	Revision 1.6  1998/07/26 04:17:05  scobb
	Added ANY_CHAN_TYPE.

	Revision 1.5  1997/11/20 19:51:28  mfm
	Solaris 2: Changed include for stdtypes.h to types.h.

	Revision 1.4  1996/12/06 01:22:01  slee
	new functions.

# Revision 1.3  1996/08/13  01:21:16  slee
# Added name of the default weight table.
#
# Revision 1.2  1996/06/11  01:35:02  slee
# Added offset_weight_table(), to add a constant analog value to the weight table.
#
# Revision 1.1  1996/05/24  00:49:34  slee
# Initial revision
#

*******************************************************************************/
/* static char dftk_weights_h[] = "$Id: dftk_weights.h,v 1.7 2000/01/11 22:48:55 skip Exp $"; */


#include <dirent.h>		// for MAXNAMLEN
#include <sys/types.h>	// for time_t
#include <tes_con.h>		// For definitions of weight channel constants

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

// Names used in system config file
#define DFTK_E_TABLE_NAME	"B_E_TABLE"	/* user block name for all E tables */
#define DFTK_WT_TABLE_PREFIX	"BWT_"		/* user block prefix for all weight tables */
						/* This is limited to 4 chars so any valid filter */
						/*  name may be used */

// Names used in run config file
#define APPLIED_E_TABLE_NAME	"B_E_table_used"/* user block name for E table applied to data */
#define APPLIED_WT_TABLE_NAME	"B_weights_used"/* user block name for wt table applied to data */


const int ANY_CHAN_TYPE 		= 999;	/* used for allowing afw to use any channel as a reference */

const int DFTK_WEIGHT_ANY_NUMBER	= -1;	// Used to specify any number of reference channel types or entries

#define	DFTK_WEIGHT_ANY_NAME		NULL	// Used as a wildcard for string parameters

// Default weight table name (dftk_weight_table::name)

#define DEFAULT_WEIGHT_TABLE_NAME	"Universal"

// These constants are delimiters and other keywords that are found in the 
// text files used to create both E and weight tables. Table-specific parameters
// are defined in the header files associated with the table type.

const int DFTK_WEIGHT_MAX_INPUT_LINE	= 2000;
#define	DFTK_WEIGHT_DELIMITER_STRING	" \t"
#define	DFTK_WEIGHT_COMMENT_TOKEN	'#'
#define	DFTK_WEIGHT_ENTRY_TOKEN		"Entry"
#define	DFTK_WEIGHT_CH_HEADER_TOKEN	"Channel"

// This constant is used by index lookup functions to indicate either
// the end of a list of index values, or the absense of an index corresponding
// to a channel name.

const long DFTK_WEIGHT_END_OF_LIST      = -1;


// Now include the E and weight table-specific files that use the constants above.

#include <dftk_E_table.h>
#include <dftk_weight_table.h>

#endif
