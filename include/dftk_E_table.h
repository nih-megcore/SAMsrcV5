#ifndef DFTK_E_TABLE_H
#define DFTK_E_TABLE_H
/*******************************************************************************

    $Id: dftk_E_table.h,v 1.2 2004/03/20 00:03:16 zirnhelt Exp $

    Name:		(file) dftk_E_table.h

    Purpose:	Defines the classes, structures and user block names for
            the system configuration file's E Tables.

    Notes for searching for E tables in system config files:

            This applies to system config files only.

            The intended use is that there will be several
            dftk_E_table structures in the system config file.
            Each will go by the name, DFTK_E_TABLE_NAME, and
            can be accessed by dftk_config::search_user_block(...)
            and its related functions.  In the E table structure
            is a filter name, filter, which is used to look up
            the related weight tables (concatenate DFTK_WT_TABLE_PREFIX
            and the filter name).


    Notes for obtaining E tables from run config files:

            There will be one and only one, and it will go
            by the name APPLIED_E_TABLE_NAME.


    Notes for accessing E table entries:

            All accesses must be via the class member functions.

    Author:		Merle F. McClelland

    Creation Date:	19 August 1998

    Copyright 1998, Biomagnetic Technologies, inc.

    $Log: dftk_E_table.h,v $
    Revision 1.2  2004/03/20 00:03:16  zirnhelt
    cast character string to get rid of warnings.

    Revision 1.1  1999/12/06 21:43:32  skip
    Initial revision


*******************************************************************************/

// static char dftk_E_table_h[] = "$Id: dftk_E_table.h,v 1.2 2004/03/20 00:03:16 zirnhelt Exp $";

/////////////////////////////////// E Table Definitions ////////////////////////////////////

// The following structures and constants are used to overlay user block data from a config file
// in order to read and write such blocks. All application access to E Tables must be through the
// class dftk_E_table, not through these structures, as the class manages memory allocation and
// deallocation as needed.

// These constants define the size of the different arrays in the dftk_E_value_WH2500 structure

const int DFTK_WH2500_NUM_E_TABLE_ENTRIES 	= 6;
const int DFTK_WH2500_NUM_LOW_GAIN_ENTRIES 	= 3;
const int DFTK_WH2500_NUM_HIGH_GAIN_ENTRIES 	= 3;
const int DFTK_WH2500_FIRST_LOW_GAIN_INDEX 	= 0;
const int DFTK_WH2500_FIRST_HIGH_GAIN_INDEX 	= 3;

// Indicies for specific channels

const int DFTK_WH2500_E_REF_MxA 		= 0;
const int DFTK_WH2500_E_REF_MyA 		= 1;
const int DFTK_WH2500_E_REF_MzA 		= 2;
const int DFTK_WH2500_E_REF_MxaA 		= 3;
const int DFTK_WH2500_E_REF_MyaA 		= 4;
const int DFTK_WH2500_E_REF_MzaA 		= 5;

// These constants define version numbers for the two different supported E Table formats
// The ALL constant is used to specify that all type of weight tables are to be acted upon
// by various functions.

const int DFTK_E_TABLE_VERSION_ALL 		= 0;
const int DFTK_E_TABLE_VERSION_WH2500 		= 1;
const int DFTK_E_TABLE_VERSION_VARIABLE_LENGTH	= 2;
// const int DFTK_E_TABLE_VERSION_LAST             = DFTK_E_TABLE_VERSION_VARIABLE_LENGTH;
const int DFTK_E_TABLE_VERSION_LAST             = 2;

// const int DFTK_E_TABLE_VERSION_THIS_SYSTEM	= DFTK_E_TABLE_VERSION_VARIABLE_LENGTH;
const int DFTK_E_TABLE_VERSION_THIS_SYSTEM	= 2;

const int DFTK_E_TABLE_FILTER_NAME_LEN		= 16;
const int DFTK_E_TABLE_RESERVED_LEN		= 28;

/// -------------------- Error Code Definitions  -------------------
///
/// These are constants for various error return status codes
/// If you add or change something here, remember to update the associated
/// array of strings in the dftk_E_table.C source file!
///

enum
{
    DFTK_E_TABLE_NO_ERROR      = 0,
    DFTK_E_TABLE_INVALID_VERSION_ERROR,
    DFTK_E_TABLE_NO_USER_BLOCK_ERROR,
    DFTK_E_TABLE_SMALL_USER_BLOCK_ERROR,
    DFTK_E_TABLE_NUM_ENTRIES_CONFIG_MISMATCH_ERROR,
    DFTK_E_TABLE_NO_VALID_CHANNELS_IN_CONFIG_ERROR,
    DFTK_E_TABLE_NO_VALID_REFERENCES_IN_CONFIG_ERROR,
    DFTK_E_TABLE_OPEN_FILE_ERROR,
    DFTK_E_TABLE_ALLOCATION_ERROR,
    DFTK_E_TABLE_PARSE_NO_HEADER_ERROR,
    DFTK_E_TABLE_PARSE_E_VALUE_NAME_ERROR,
    DFTK_E_TABLE_PARSE_E_VALUE_INVALID_NAME_ERROR,
    DFTK_E_TABLE_PARSE_E_VALUE_DUPLICATE_NAME_ERROR,
    DFTK_E_TABLE_PARSE_E_VALUE_INVALID_NUMBER_ERROR,
    DFTK_E_TABLE_PARSE_TOO_MANY_E_VALUE_COLUMNS_ERROR,
    DFTK_E_TABLE_PARSE_TOO_MANY_HEADER_COLUMNS_ERROR,
    DFTK_E_TABLE_PARSE_INVALID_ENTRY_NUMBER_ERROR,
    DFTK_E_TABLE_PARSE_CHANNEL_INVALID_NAME_ERROR,
    DFTK_E_TABLE_PARSE_CHANNEL_DUPLICATE_NAME_ERROR,
    DFTK_E_TABLE_PARSE_TOO_FEW_COLUMNS_OR_ROWS_ERROR,
    DFTK_E_TABLE_INVALID_CHANNEL_NAME_ERROR,
    DFTK_E_TABLE_INVALID_E_VALUE_NAME_ERROR,
    DFTK_E_TABLE_INVALID_NAME_ERROR,
    DFTK_E_TABLE_LAST_ERROR                    	// Don't put anything below this!

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
/// Entry	Channel 	<ref ch name> 	<ref ch name> 	...
/// <e>		<ch name>	<E value> 	<E value>  	...
/// .
/// .
/// .
///
/// All fields are white space-delimited (spaces, tabs).
///
/// The first non-comment line is a header consisting of the word "Entry" followed
/// by the word "Channel", and then channel name strings representing the names of
/// the channels associated with the columns of E values below. Note that the
/// order of these channel columns is not important, as the parsing function will
/// use the channel name string at the top of the column to locate the channel's
/// entry in the table. The entries in the table are sorted in config file order,
/// both in the E Channel axis and in the entry (entry channel) axis. Duplicated,
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
/// purposes - ignored), the Channel Name, and then an E value for each of the columns
/// in floating point format. The number of these lines must match the number of
/// entries in the table, which is determined by the type of table and the config file.
/// For WH2500 type tables, the text file must contain entries for all MEG channels,
/// but may contain entries for other channels as well. For variable type tables,
/// only MEG channels are permitted, since the tables contain only MEG channel entries..
///
/// The numbers must be valid floating point format, but no other checking is
/// performed (e.g. range checking). The Entry number must be a valid decimal
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
/// followed by an array of arrays of 6 floats each. No channel names are
/// stored, but the length of the array is variable. When such a table is read
/// or written, it is converted to/from the variable format at is used for all
/// weight tables when they are resident in memory as objects.
///
/// Likewise, variable-length tables are converted to/from a packed format
/// in user blocks.
///
/// The dftk_E_table_header structure is the common header for both WH2500
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
/// "num_entries" arrays of 6 floating point values. Note that the number of entries
/// is normally the total number of channels in the system, with non-MEG channel
/// entries set to zero. No channel name arrays exist in the fixed-size design.
/// The order of channels in the list must match the CFO order of the config file
/// with which the weight table is being used. However, there is no way to verify
/// this.
///
/// ----------------------- Variable Size User Blocks ------------------------
///
/// The array contains entries for only MEG channels. The MEG channel names array
/// associates he entries in the list with the channel to which they apply. The use of names
/// instead of CFO numbers permits the use of the tables with any config file,
/// since channel order in the table is not based on a particular config file.
/// Applications that need to access the tables must look up the entries by name,
/// rather than parsing the list sequentially and assuming a particular order.
///
/// Immediately following the header are four arrays:
///
/// 1. Array of Entry Channel Names, "num_entries" entries, each
///	DFTK_CHANNEL_STRLEN in length, in the order in which the E values
///	associated with each channel are stored in the array (in the
///	"num_entries" axis).
/// 2. Array of E Value Channel Names, "num_E_values" entries, each
///	DFTK_CHANNEL_STRLEN in length, in the order in which the corresponding
///	entries appear in the E Values array below.
/// 3. E Values array, "num_E_values" wide by "num_entries" high, in float format.
///	Column order must match the order of names in the E Values Channel Names
///	array.
///
/// Note that arrays 1 and 2 must be integer multiples of 4 bytes so that
/// the following values array is aligned properly in memory (floats must be on a
/// 4 byte boundary). This is assured since DFTK_CHANNEL_STRLEN is 16 bytes.
///
/// --------------------------------------------------------------------------
///
/// When either type of block is read into an object, it is converted into the
/// memory format that is defined by the object. This means that the data in memory
/// is not stored as a single large block.
///
/// Access to both types of arrays (WH2500 fixed and the variable type)
/// should be performed by accessor functions provided by the class.

typedef struct {
	float_t			E[DFTK_WH2500_NUM_E_TABLE_ENTRIES];		// MxA, MyA, MzA, MxaA, MyaA, MzaA
} dftk_E_WH2500;


/// This structure is used for both the WH2500 fixed-size and the variable-sized designs.
///
/// For the variable length design, the entry_size field is set to num_E_values * sizeof(float).
///
/// Access to these arrays should be performed by accessor functions provided by the library
///

typedef struct {
	int32_t			version;								// version of the table
	int32_t			entry_size;								// size of one entry in the table
	int32_t			num_entries;							// # of entries
	char_t			filter[DFTK_E_TABLE_FILTER_NAME_LEN];	// high-pass filter name
	int32_t			num_E_values;							// The number of E values per entry
	char_t			reserved[DFTK_E_TABLE_RESERVED_LEN];
} dftk_E_table_header;

extern char *dftk_E_table_version_strings[];

/// A full set of classes is defined to create and access E tables of either format.
/// There are several constructors provided:
///
/// 1. Given a config file user block user_data_space pointer, create an E table object
/// 	of the appropriate type.
/// 2. Given basic parameters, create an empty E table object
///
/// Member functions are provided to load and save the contents of an E table object
/// into a user block.

class dftk_E_table : public generic, public dftk
{
    private:
    // This data will appear as-is in the beginning of the user block

        dftk_E_table_header	header;

        dftk_config		*config;	// The config file associated with this table
        dftk_user_block		*user_block;	// The user block associated with this table
        char			user_block_name[DFTK_PROCESS_NAMELEN+1];// The name of the current user
                                                // block, or the name of the user block that should be
                                                // created if one does not exist.

        time_t                  timestamp;      // Timestamp from the user block

    // These arrays do not appear in the WH2500 fixed-size user blocks

        char			**entry_channel_name_array;
        char			**E_value_channel_name_array;

    /// In the user block, this array is written after the channel name array
    /// for the variable-length design, or after the header for the WH2500 fixed-size design
        float			**E_values_array;

        /// This is private so users can't pull the rug out from under us by
        /// changing the config file.
        inline void set_config(dftk_config *new_config)
        {
            config = new_config;
        };

        /// These routines perform the actual creation of the two types of in-memory user
        /// block images, and are called from create_user_block_image()

        dftk_E_table_header *create_WH2500_user_block_image(
                    int &status, unsigned long *size_in_bytes = NULL);
        dftk_E_table_header *create_variable_user_block_image(
                     int &status, unsigned long *size_in_bytes = NULL);

        /// This routine is used by the save function to allocate memory for the two types of tables
        /// for storage in user blocks. It allocates a single memory block and fill-in the header
        /// and transfer the object's data into the format appropriate for the table type. The
        /// size_in_bytes parameter, if passed, is filled-in with the size of the allocated block.

        dftk_E_table_header *create_user_block_image(int &status, unsigned long *size_in_bytes = NULL);

        /// These routines perform the actual translation of the two in-memory formats of user
        /// blocks into data compatible with the class. Called by load_user_block_image().

        int load_WH2500_user_block_image(dftk_E_table_header *E_table, unsigned long size_of_block);
        int load_variable_user_block_image(dftk_E_table_header *E_table, unsigned long size_of_block);

        /// This routine will take a memory block pointed-to by the E_table parameter, interpret
        /// it, and if valid, load it into this object's data structures. Existing data will be
        /// deleted if the passed block is valid. The size_of_block parameter is used to validate
        /// the header fields.

        int load_user_block_image(dftk_E_table_header *E_table, unsigned long size_of_block);

        /// Routine to allocate the arrays
        void                    allocate_arrays();

    // Routines to parse the config file for channel information

        long			fill_entry_channel_info(int &status);
        long			fill_E_channel_info(int &status);

    // Routines to allocate memory and do basic setup

        void 			init_WH2500(int &status, char *filter_name);
        void 			init_variable(int &status, char *filter_name);
        void			init(int &status, long version, char *filter_name);

        // Routines to deallocate memory when reconstruction takes place, or the object
        // is destroyed
        void                    null_pointers();
        void                    freeall();

       /// This routine will validate whether the passed parameters are
       /// valid for a variable length E table
       int validate_parameters(long entry, long E_index);

    public:

    ///// CONSTRUCTORS

        /// Copy constructor
        dftk_E_table(dftk_E_table &e);

        /// Create an E Table from user block data
        dftk_E_table(int &status, dftk_config *config, dftk_user_block *block);

        /// Create an E Table from a text file using the parse routine. Errors
        /// are returned in status. line_number will be set to the line number
        /// in the input file where a parsing error occurred.
        dftk_E_table(int &status, dftk_config *config, char *pathname, int &line_number);

        /// Create an E Table from a open file using the parse routine. Errors
        /// are returned in status. line_number will be set to the line number
        /// in the input file where a parsing error occurred.
        dftk_E_table(int &status, dftk_config *config, FILE *fp, int &line_number);

        /// Create an empty table based on the passed version number. The number of
        /// entries and the number of E values will be derived from the passed config file,
        /// as will the channel names.
        dftk_E_table(int &status, dftk_config *config, long version, char *filter_name);

        /// DESTRUCTOR
        ~dftk_E_table();

        /// Assignment operator
        dftk_E_table operator=(dftk_E_table &e);

        /// Saves the table in a user block - will return result in status variable
        /// If the block was loaded from a user block, the same user block is updated.
        /// If no user block exists, a new one is created. In both cases the user block
        /// timestamp is updated. If the flag indicates that the table is to be an
        /// applied table (i.e. one in a run config file), the special user block name
        /// is used. Others`wise, if the object was created from a previous user block,
        /// that name is used, else a new name is automatically generated.

        int save(int applied_table = FALSE);

        /// Loads the table from the passed user block, destroying all existing data
        int load(dftk_user_block *block);

        /// This routine will parse the specified text file and load the contents into this
        /// extisting object, replacing what is already here. The format of the data in the
        /// object (WH2500 fixed or variable-length) will depend entirely on the format of the
        /// contents of the text file. Errors will be returned if the file is invalid. The
        /// line_number variable will be set to the line number inthe input file where the
        /// error occurred.

        int parse(char *pathname, int &line_number);

        /// This routine will parse the specified open text file and load the contents into this
        /// extisting object, replacing what is already here. The format of the data in the
        /// object (WH2500 fixed or variable-length) will depend entirely on the format of the
        /// contents of the text file. Errors will be returned if the file is invalid.

        int parse(FILE *fp, int &line_number);

        /// This routine will parse the standard input stream and load the contents into this
        /// extisting object, replacing what is already here. The format of the data in the
        /// object (WH2500 fixed or variable-length) will depend entirely on the format of the
        /// contents of the input text. Errors will be returned if the input is invalid
        inline int		parse(int &line_number)
        {
            return (parse(stdin, line_number));
        };

        /// This routine will print the E table in a formatted text report to the
        /// specified file. The optional second parameter specifies that the table header
        /// be formatted and printed before the E values. The entry_limit field can be
        /// used to limit the report to the first "n" entries.
        int print_report(FILE *fp, int print_table_header_info = FALSE,
                         long entry_limit = DFTK_WEIGHT_ANY_NUMBER);

        /// This routine will print the E table in a formatted text report to
        /// the specified file name (must be a fully-qualified path to a file).
        int print_report(char *pathname, int print_table_header_info = FALSE,
                                      long entry_limit = DFTK_WEIGHT_ANY_NUMBER);

        /// This routine will print the E table in a formatted text report to
        /// standard output
        inline int print_report(int print_table_header_info = FALSE,
                                long entry_limit = DFTK_WEIGHT_ANY_NUMBER)
        {
            return(print_report(stdout, print_table_header_info, entry_limit));
        };

        /// Data Member accessor functions - note that these cannot be changed
        /// once the object is created
        /// so there are no corresponding "set" functions
        inline long get_version()
        {
            return (header.version);
        };

        inline char *get_version_name()
        {
            return (header.version <= DFTK_E_TABLE_VERSION_LAST ?
                              dftk_E_table_version_strings[header.version - 1] :
                    (char *)"Invalid");
        };
        inline long get_entry_size()
        {
            return (header.entry_size);
        };
        inline long get_num_entries()
        {
            return (header.num_entries);
        }
        inline long get_num_E_values()
        {
            return (header.num_E_values);
        };

        inline time_t get_timestamp()
        {
            return (timestamp);
        };


        // Data Member accessor functions - these members can be modified after the object is created,
        // so there are both get and set versions.

        inline void set_filter(char *filter)
        {
            strncpy(header.filter, filter, DFTK_E_TABLE_FILTER_NAME_LEN);
        };
        inline char *get_filter()
        {
            return (header.filter);
        };

    // All application access to E Table arrays for either fixed or variable-sized arrays should be done
    // through these get and set functions.

        float get_E_value(long entry, long E_index);
        void set_E_value(long entry, long E_index, float value);

    /// Given an entry index, get the entry channel name - returns -1 if invalid index

        char *get_entry_ch_name(long entry);

    /// Given a channel name, get the index into the array for the channel - returns -1 if invalid name

        long get_entry_ch_index(char *channel_name);

    /// Given an entry index, set the associated channel name.
    /// Entry index and channel_name length must be valid
    /// on an error is returned.

        int set_entry_ch_name(long entry, char *channel_name);

    /// Given an E_index value , get the string for the channel
    /// name - returns NULL if invalid E_index

        char *get_E_value_ch_name(long E_index);

    /// Given a channel name, get the index into the E array for the
    /// channel - returns -1 if invalid name

        long get_E_value_ch_index(char *channel_name);

    /// Given an E_index, set the associated channel name. E_index and
    /// channel_name length must be valid on an error is returned.

        int set_E_value_ch_name(long E_index, char *channel_name);

    /// Shortcut routine that will clear all of E values through the entire array

        void clear_E_values();

    /// Offset the E values. The parameter is a pointer to a list of offset
    /// values, and the number of entries in the list must match the number of
    /// E values channels in the table. Each E value is offset by its corresponding
    /// value in the offset_list.

        void offset_E_values(float *offset_list);

    /// This is a routine that will take a status returned by one of the class functions
    /// and return a printable string

        static char *get_error_string(int status);

};

#endif
