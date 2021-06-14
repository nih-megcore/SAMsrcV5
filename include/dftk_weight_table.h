#ifndef DFTK_WEIGHT_TABLE_H
#define DFTK_WEIGHT_TABLE_H
/*******************************************************************************

    $Id: dftk_weight_table.h,v 1.3 2004/03/20 00:03:38 zirnhelt Exp $

    Name:  (file) dftk_weight_table.h

    Purpose: Defines the classes, structures and user block names for
            the system configuration file's weight tables.

    Notes for searching for weight tables in system config files:

            This applies to system config files only.

            To find all weight tables, one can also go through all
            user blocks and check the first few characters of the
            name with DFTK_WT_TABLE_PREFIX.  i.e.:

            strncmp(block_name, DFTK_WT_TABLE_PREFIX, strlen(DFTK_WT_TABLE_PREFIX)


    Notes for obtaining the weight table from run config files:

            There will be one and only one, and it will go
            by the name APPLIED_WT_TABLE_NAME.


    Author:  Merle F. McClelland

    Creation Date: 19 August 1998

    Copyright 1998, Biomagnetic Technologies, inc.

    $Log: dftk_weight_table.h,v $
    Revision 1.3  2004/03/20 00:03:38  zirnhelt
    cast character string to get rid of warnings.

    Revision 1.2  2000/04/06 22:11:54  skip
    Added error code for setting names.

    Revision 1.1  1999/03/08 21:38:04  wge
    Initial revision


*******************************************************************************/
/// static char dftk_weight_table_h[] = "$Id: dftk_weight_table.h,v 1.3 2004/03/20 00:03:38 zirnhelt Exp $";

/// This offset must be added to the weight values when they are sent to the DAS

const short DFTK_WEIGHT_ANALOG_OFFSET = 2048;

/// ------------ WH2500 Fixed-size Table Constant Definitions  --------------

/// These constants define the number of different types of weights for the WH2500
/// fixed-size design. Note that unlike previously, the DSP weight array is now a
/// single array, instead of three separate arrays for low, high, and gradiometer
/// weights. This makes the design consistent with the variable-length tables, and
/// is compatible with existing tables (since the single array overlays the same
/// space as the previous 3 arrays in the table).

// const int DFTK_WH2500_NUM_ANALOG_WEIGHTS = 3;
const int DFTK_WH2500_NUM_LOW_GAIN_WEIGHTS = 3;
const int DFTK_WH2500_NUM_HIGH_GAIN_WEIGHTS = 3;
const int DFTK_WH2500_NUM_GRADIOMETER_WEIGHTS = 5;
// const int DFTK_WH2500_NUM_DSP_WEIGHTS  = 11;
// const int DFTK_WH2500_NUM_DSP_WEIGHTS  = (DFTK_WH2500_NUM_LOW_GAIN_WEIGHTS +
//       DFTK_WH2500_NUM_HIGH_GAIN_WEIGHTS +
//      DFTK_WH2500_NUM_GRADIOMETER_WEIGHTS);

/// analog array

const long DFTK_WH2500_W_REF_ANALOG_MxA = 0;
const long DFTK_WH2500_W_REF_ANALOG_MyA = 1;
const long DFTK_WH2500_W_REF_ANALOG_MzA = 2;

/// low_gain_mag channel indexes

const long DFTK_WH2500_W_REF_LOW_MxA  = 0;
const long DFTK_WH2500_W_REF_LOW_MyA  = 1;
const long DFTK_WH2500_W_REF_LOW_MzA  = 2;

/// high_gain_mag channel indexes

const long DFTK_WH2500_W_REF_HIGH_MxaA  = 3;
const long DFTK_WH2500_W_REF_HIGH_MyaA  = 4;
const long DFTK_WH2500_W_REF_HIGH_MzaA  = 5;

/// grad channel indexes

const long DFTK_WH2500_W_REF_GRAD_GxxA  = 6;
const long DFTK_WH2500_W_REF_GRAD_GyyA  = 7;
const long DFTK_WH2500_W_REF_GRAD_GyxA  = 8;
const long DFTK_WH2500_W_REF_GRAD_GzxA  = 9;
const long DFTK_WH2500_W_REF_GRAD_GzyA  = 10;


/// Total Number of available Channels on a 2500 system NOT including the derived chans
const long DFTK_WH2500_TOTAL_ENTRIES = 258;

/// ------------ WH2500 Fixed-size Table Structure Definitions  --------------

/// This structure is used to read and write WH2500-style weight tables

// struct dftk_weight_WH2500
typedef struct {
	int16_t		analog[3];		// Analog / MDAC weights
	int16_t		unused;			// unused but needed for structure padding
	float_t		dsp[11];		// DSP weights
} dftk_weight_WH2500;

/// --------------- Variable-size Table Constant Definitions  --------------

/// These constants define the size of the different arrays in the dftk_weight_table structure
/// Note that for the variable length design, there is only a single, large, DSP weights array, instead of the
/// separate arrays for high/low/gradiometer weights. This makes the array design consistent with
/// the DAS. The interpretation as to which entry is which channel is determined by the analog and dsp
/// channel name lists located after the header. Routines are provided to get index lists for the different
/// types of channels (low gain, high gain, gradiometer).

const int DFTK_NUM_DSP_WEIGHT_TYPES = 3;

/// The assumption is that the constants "MAX_ANALOG_CHANS" and
/// "MAX_REFERENCE_CHANS" are defined in the header files for the
/// real-time system, and thus for each system type, this file will
/// use the current stage's values. Although the variable-sized design
/// doesn't use these constants for array bounds, applications for a
/// particular stage can use these to validate weight tables. When new
/// system types are defined, the constants for system types other than
/// the type of system for the current stage should be added here, using
/// the same naming convention (replace WH3600 with the system design name).
/// The constants DFTK_CURRENT_NUM_ANALOG_WEIGHTS and DFTK_CURRENT_NUM_DSP_WEIGHTS
/// should be set to the values for the system in the current stage.

const int DFTK_WH3600_NUM_ANALOG_WEIGHTS  = MAX_ANALOG_CHANS;
const int DFTK_WH3600_NUM_DSP_WEIGHTS  = MAX_REFERENCE_CHANS;

const int DFTK_CURRENT_NUM_ANALOG_WEIGHTS  = DFTK_WH3600_NUM_ANALOG_WEIGHTS;
const int DFTK_CURRENT_NUM_DSP_WEIGHTS  = DFTK_WH3600_NUM_DSP_WEIGHTS;

/// These constants define the version numbers for the two currently-define table types.
/// Unless the fundamental structure of the weight tables changes, future system types
/// should be able to use the VARIABLE_LENGTH type. The ALL constant is used to specify
/// that all types of tables are to be acted upon in various routines.

const int DFTK_WEIGHT_TABLE_VERSION_ALL   = 0;
const int DFTK_WEIGHT_TABLE_VERSION_WH2500  = 1;
const int DFTK_WEIGHT_TABLE_VERSION_VARIABLE_LENGTH = 2;
const int DFTK_WEIGHT_TABLE_VERSION_LAST  = DFTK_WEIGHT_TABLE_VERSION_VARIABLE_LENGTH;

const int DFTK_WEIGHT_TABLE_VERSION_THIS_SYSTEM       = DFTK_WEIGHT_TABLE_VERSION_VARIABLE_LENGTH;

const int DFTK_WEIGHT_TABLE_NAME_LEN   = 32;
const int DFTK_WEIGHT_TABLE_DESCRIPTION_LEN  = 80;
const int DFTK_WEIGHT_TABLE_RESERVED_LEN  = 72;

/// -------------------- Error Code Definitions  -------------------
///
/// These are constants for various error return status codes
/// If you add or change something here, remember to update the associated
/// array of strings in the dftk_weight_table.C source file!
///

enum
{
    DFTK_WEIGHT_TABLE_NO_ERROR      = 0,
    DFTK_WEIGHT_TABLE_INVALID_VERSION_ERROR,
    DFTK_WEIGHT_TABLE_NO_USER_BLOCK_ERROR,
    DFTK_WEIGHT_TABLE_SMALL_USER_BLOCK_ERROR,
    DFTK_WEIGHT_TABLE_NUM_ENTRIES_CONFIG_MISMATCH_ERROR,
    DFTK_WEIGHT_TABLE_NO_VALID_CHANNELS_IN_CONFIG_ERROR,
    DFTK_WEIGHT_TABLE_NO_VALID_REFERENCES_IN_CONFIG_ERROR,
    DFTK_WEIGHT_TABLE_OPEN_FILE_ERROR,
    DFTK_WEIGHT_TABLE_ALLOCATION_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_NO_HEADER_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_ANALOG_INVALID_NAME_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_ANALOG_DUPLICATE_NAME_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_ANALOG_INVALID_NUMBER_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_ANALOG_NUMBER_OUT_OF_RANGE_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_DSP_INVALID_NAME_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_DSP_DUPLICATE_NAME_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_DSP_INVALID_NUMBER_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_TOO_MANY_WEIGHT_COLUMNS_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_TOO_MANY_HEADER_COLUMNS_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_INVALID_ENTRY_NUMBER_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_CHANNEL_INVALID_NAME_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_CHANNEL_DUPLICATE_NAME_ERROR,
    DFTK_WEIGHT_TABLE_PARSE_TOO_FEW_COLUMNS_OR_ROWS_ERROR,
    DFTK_WEIGHT_TABLE_INVALID_CHANNEL_NAME_ERROR,
    DFTK_WEIGHT_TABLE_INVALID_ANALOG_NAME_ERROR,
    DFTK_WEIGHT_TABLE_INVALID_DSP_NAME_ERROR,
    DFTK_WEIGHT_TABLE_INVALID_NAME_ERROR,
    DFTK_WEIGHT_TABLE_INVALID_DESCRIPTION_ERROR,
    DFTK_WEIGHT_TABLE_CANT_SET_NAMES,
    DFTK_WEIGHT_TABLE_LAST_ERROR          /// Don't put anything below this!
};



/// -------------------- Text File Format Definitions  -------------------
///
/// These definitions are used to define the text file format that must be
/// adhered-to when using the parse member function. Note that is is not the
/// same format as older WH2500 style files. All WH2500 style files must
/// be updated to be used with the class.
///
/// The text weight tables must be in the following format:
///
/// # (lines beginning with "#" are ignored as comments, as are blank lines
///
/// Entry        Channel         <analog ch name>   ...  <dsp ch name>   ...
/// <e>          <ch name>       <weight value>    ... <weight value>  ...
/// .
/// .
/// .
///
/// All fields are white space-delimited (spaces, tabs).
///
/// The first non-comment line is a header consisting of the word "Entry" followed
/// by the word "Channel", and then channel name strings representing the names of
/// the channels associated with the columns of weight values below. The first set
/// of names and associated columns must be the analog weight channels, followed
/// by the set of DSP (digital) weight channels. Note that within these two groups the
/// order of these channel columns is not important, as the parsing function will
/// use the channel name string at the top of the column to locate the channel's
/// entry in the table. The entries in the table are sorted in config file order,
/// both in the weight Channel axis and in the entry (entry channel) axis. Duplicated,
/// extra, or insufficient entries (as determined by the number of channels found
/// in the config file) are flagged as errors. Channels that are not valid
/// reference channels are also flagged as errors. In the case of system that
/// have alternate reference channel sources, the channel names must match the
/// specific source channel name for each reference. For example, for WH3600,
/// which has Left, Right, and Center sources for the Mx, My, and Mz channels,
/// and has a user block in the config file that selects the source for each,
/// the labels might be "MxLA", "MyCA", and "MzRA". Both low and high gain
/// channels must be listed (channel names with and without the "a" modifier).
/// This ensures that the table corresponds to the selection of reference channels
/// in the config file.
///
/// Following the header are entries consisting of an Entry number "e" (for historical
/// purposes - ignored), the Channel Name, and then two sets of weight values. The
/// first set, corresponding to the analog channel names, are short integer values.
/// These values are checked for proper format and range (+32767..-32768).  The second
/// set of weight values correspond to the DSP channel names and are in floating point
/// format. These must be valid floating point numbers but no range checking is
/// performed. The number of shorts and floats must match the number of analog and
/// DSP channel names at in the header.
///
/// The number of entry lines must match the number of entries in the table, which
/// is determined by the type of table and the config file.
/// For WH2500 type tables, the text file must contain entries for all MEG channels,
/// but may contain entries for other channels as well. For variable type tables,
/// only MEG channels are permitted, since the tables contain only MEG channel entries..
///
/// The Entry number must be a valid decimal
/// integer, but also is not checked. The channel name must be valid, and must
/// represent a channel in the associated config file. Duplicated entries are
/// not permitted.
///
/// If a parsing error occurs, the parse routine will return an error code and
/// a line number where the first error was encountered, to assist in debugging.
///
/// For variable-length tables, the table is created with exactly the number of
/// entries found in the list. For WH2500 tables, additional entries will be created
/// to fill-in the required full set of 258 entries, based on the number of channels
/// in the associated config file. The file may contain dummy entries for the
/// non-MEG channels - these are stored but not used by any weight processing
/// applications.
///
/// -------------------- Description of User Block Formats -------------------
///
/// There is a difference between the way the tables are stored in the user blocks
/// and they way they are managed in memory. For WH2500 type tables, we maintain
/// the existing previously-specified user block format, consisting of a header
/// followed by an array of dftk_weight_WH2500 structures. No channel names are
/// stored, but the length of the array is variable. When such a table is read
/// or written, it is converted to/from the variable format at is used for all
/// weight tables when they are resident in memory as objects.
///
/// Likewise, variable-length tables are converted to/from a packed format
/// in user blocks.
///
/// The dftk_weight_table structure is the common header for both WH2500
/// fixed-size arrays and the variable-sized design. It is also the common
/// data structure between the objects in memory and the data stored in the
/// user blocks. The header is always stored as-is when writing a user block.
///
/// In the user block, the header is followed by the weight table data, as
/// described below:
///
/// --------------------- WH2500 Fixed Size User Blocks ----------------------
///
/// For the fixed-size design, the header is immediately followed by an array of
/// "num_entries" dftk_weight_WH2500 structures. Note that the number of entries
/// is normally the total number of channels in the system, with non-MEG channel
/// entries set to zero. No channel name arrays exist in the fixed-size design.
/// The order of channels in the list must match the CFO order of the config file
/// with which the weight table is being used. However, there is no way to verify
/// this.
///
/// ----------------------- Variable Size User Blocks ------------------------
///
/// For the variable-length design, the short value and the float values for each
/// channel are separated into two separate arrays. This allows for correct data
/// alignment regardless of the number of items in the array. Also, the arrays
/// contain entries for only MEG channels. The entry channel names array associates
/// the entries in the list with the channel to which they apply. The use of names
/// instead of CFO numbers permits the use of the tables with any config file,
/// since channel order in the table is not based on a particular config file.
/// Applications that need to access the tables must look up the entries by name,
/// rather than parsing the list sequentially and assuming a particular order.
///
/// Immediately following the header are four arrays:
///
/// 1. Array of Entry Channel Names, "num_entries" entries, each
/// DFTK_CHANNEL_STRLEN in length, in the order in which the weights
/// associated with each channel are stored in the array (in the
/// "num_entries" axis).
/// 2. Array of Analog Channel Names, "num_analog" entries, each
/// DFTK_CHANNEL_STRLEN in length, in the order in which the corresponding
/// entries appear in the Analog Weights array below.
/// 3. Array of DSP Channel Names, "num_dsp" entries, each
/// DFTK_CHANNEL_STRLEN in length, in the order in which the corresponding
/// entries appear in the DSP Weights array below.
/// 4. DSP Weights array, "num_dsp" wide by "num_entries" high, in float format.
/// Column order must match the order of names in the DSP Channel Names
/// array.
/// 5. Analog Weights array, "num_analog" wide by "num_entries" high, in
/// short format. Column order must match the order of names in the
/// Analog Channel Names array.
///
/// Note that arrays 1, 2,  and 3 must be integer multiples of 4 bytes so that
/// the following DSP array is aligned properly in memory (floats must be on a
/// 4 byte boundary). This is assured since DFTK_CHANNEL_STRLEN is 16 bytes.
/// Likewise, this is why the Analog Weights appear after the DSP weights, so
/// that the DSP weights are always aligned, regardless
/// of the number of shorts in the Analog array (which may be 0!).
///
/// --------------------------------------------------------------------------
///
/// When either type of block is read into an object, it is converted into the
/// memory format that is defined by the object. This means that the data in memory
/// is not stored as a single large block.
///
/// Access to both types of arrays (WH2500 fixed and the variable type)
/// should be performed by accessor functions provided by the class.

/// ------------------- Header Structure Format Definitions  -----------------
///
/// The common header for both user blocks and the objects
/*
typedef struct {
    long   version;
    long   entry_size; /// Not used in variable-length design -
                        /// must calculate from num_analog and num_dsp
    long   num_entries;
    char   name[DFTK_WEIGHT_TABLE_NAME_LEN];
    char   description[DFTK_WEIGHT_TABLE_DESCRIPTION_LEN];
    long   num_analog; /// Not used in WH2500 user blocks
    long   num_dsp; /// Not used in WH2500 user blocks
    char   reserved[DFTK_WEIGHT_TABLE_RESERVED_LEN];
} dftk_weight_table_header;
*/
#endif
