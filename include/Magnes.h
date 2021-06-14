#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/types.h>

#define MEGTYPE			1		// primary MEG sensor
#define EEGTYPE			2		// EEG channe;
#define REFERENCE		3		// reference magnetometer or gradiometer
#define	EXTERNAL		4		// external channel
#define	TRIGGER			5		// "TRIGGER" and "RESPONSE" channels
#define	UTILITY			6		// utility such as current
#define	DERIVED			7		// derived channel (could this be for SAM?)
#define	SHORTED			8		// shorted channel

#define FORMAT_SHORT	1
#define FORMAT_LONG		2
#define FORMAT_FLOAT	3
#define FORMAT_DOUBLE	4

#define FILE_TYPE_STRLEN		(5)
#define DFTK_HEADER_STRLEN		(16)
#define DFTK_CONFIG_STRLEN		(32)
#define DFTK_DAP_STRLEN			(16)
#define DFTK_PROCESS_NAMELEN	(20)
#define DFTK_HEADER_STRLEN		(16)
#define DFTK_CHANNEL_STRLEN		(16)
#define DFTK_EVENT_STRLEN		(16)
#define DFTK_CHANNEL_REF_STRLEN	(16)
#define FILE_TYPE_STRLEN		(5)
#define CMINSLICES				(100)	//The number of slices the data buffer is limited to.


typedef struct {		// ok
	int16_t				version;						// software version
	char				sitename[DFTK_CONFIG_STRLEN];	// system location
	char				dap_hostname[DFTK_DAP_STRLEN];	// name of associated dap
	int16_t				sys_type;						// system type in {9003, 9012, ...}
	int32_t				sys_options;					// optional features of the system
	int16_t				supply_freq;					// power supply frequency (50/60)
	int16_t				total_chans;					// number of (FET) channels in the system
	float_t				system_fixed_gain;
	float_t				volts_per_bit;
	int16_t				total_sensors;
	int16_t				total_user_blocks;
	int16_t				next_derived_channel_number;
	int32_t				checksum;
	char				reserved[32];
} dftk_config_data;

typedef struct {
	double_t			m[4][4];
} Xfm;

typedef struct {										// compiler says sizeof(Pnt) = 24 bytes
	double_t			p[3];
} Pnt;

typedef struct {
	double_t			v[3];
} Vec;

typedef struct {		// ok
	int32_t				nbytes;
	char				type[DFTK_PROCESS_NAMELEN];
	int32_t				checksum;
} dftk_process_header;

typedef struct {		// ok
	dftk_process_header	hdr;
	char				user[32];
	int32_t				timestamp;
	int32_t				user_space_size;
	char				reserved[32];		// was 32
} dftk_user_block_data;

typedef struct {		// ok
	char				name[DFTK_CHANNEL_STRLEN];			// channel name
	int16_t				chan_no;							// channel #
	uint16_t			type;
	int16_t				sensor_no;
	float_t				gain;
	float_t				units_per_bit;
	char				yaxis_label[DFTK_CHANNEL_STRLEN];	// y-axis units
	double_t			aar_val;							// auto artifact rejection value
	int32_t				checksum;
	char				reserved[32];
} dftk_channel_data;

typedef struct {		// ok
	Pnt					position;
	Vec					direction;
	double_t			radius;
	double_t			wire_radius;
	int16_t				turns;
	int32_t				checksum;
//	char				reserved[32];
	char				reserved[24];	// changed to obtain alignment
} dftk_loop_data;

typedef struct {		// ok
	int32_t				size;
	int32_t				checksum;
	char				reserved[32];
} dftk_device_header;

typedef struct {		// ok
	dftk_device_header	hdr;
	float_t				inductance;
	char				pad1[4];		// new padding
	Xfm					transform;
	int16_t				xform_flag;
	int16_t				total_loops;
	char				reserved[32];
} dftk_meg_device_data;

typedef struct {		// ok
	dftk_device_header	hdr;
	float_t				impedance;
//	char				pad1[4];		// new padding
	Xfm					transform;
	char				reserved[32];
} dftk_eeg_device_data;

typedef struct {		// ok
	dftk_device_header	hdr;
	int32_t				user_space_size;
	char				reserved[32];
} dftk_external_device_data;

typedef struct {		// ok
	dftk_device_header	hdr;
	int32_t				user_space_size;
	char				reserved[32];
} dftk_trigger_device_data;

typedef struct {		// ok
	dftk_device_header	hdr;
	int32_t				user_space_size;
	char				reserved[32];
} dftk_utility_device_data;

typedef struct {		// ok
	dftk_device_header	hdr;
	int32_t				user_space_size;
	char				reserved[32];
} dftk_derived_device_data;

typedef struct {		// ok
	dftk_device_header	hdr;
	char				reserved[32];
} dftk_shorted_device_data;

typedef struct {		// ok
	int16_t				version;
	char				file_type[FILE_TYPE_STRLEN];
	int16_t				data_format;
	int16_t				acq_mode;
	uint32_t			total_epochs;
	uint32_t			input_epochs;
	uint32_t			total_events;
	uint32_t			total_fixed_events;
	float_t				sample_period;
	char				xaxis_label[DFTK_HEADER_STRLEN];
	uint32_t			total_processes;
	int16_t				total_chans;
	uint32_t			checksum;
	uint32_t			total_ed_classes;
	int16_t				total_associated_files;
	int16_t				last_file_index;
	int32_t				timestamp;
	char				reserved[24];
} dftk_header_data;

typedef struct {		// ok
	int32_t				pts_in_epoch;
	float_t				epoch_duration;
	float_t				expected_iti;
	float_t				actual_iti;
	int32_t				total_var_events;
	int32_t				checksum;
	int32_t				epoch_timestamp;	// start time relative to the start of acqusition in time slices
	char				reserved[28];
} dftk_epoch_data;

typedef struct {		// ok
	char				event_name[DFTK_EVENT_STRLEN];
	float_t				start_lat;
	float_t				end_lat;
	float_t				step_size;
	int16_t				fixed_event;
	int32_t				checksum;
	char				reserved[32];
} dftk_event_data;

typedef struct {		// ok
	char				chan_label[DFTK_CHANNEL_REF_STRLEN];
	int16_t				chan_no;
	int16_t				attributes;
	float_t				scale;
	char				yaxis_label[DFTK_CHANNEL_REF_STRLEN];
	int16_t				valid_min_max;
	char				pad1[6];	//new padding
	double_t			ymin;
	double_t 			ymax;
	int32_t				index;
	int32_t				checksum;
	char				reserved[32];
} dftk_channel_ref_data;

typedef struct {		// ok
	dftk_process_header	hdr;
	char				user[32];
	time_t				timestamp;
	char				filename[256];
	int32_t				total_steps;
	char				reserved[32];
} dftk_proc_data;

typedef struct {		// ok
	dftk_process_header	hdr;
	float_t				frequency;
	char				reserved[32];
} dftk_filter_data;

typedef struct {		// ok
	dftk_process_header	hdr;
	float_t				high_frequency;
	float_t				low_frequency;
	char				reserved[32];
} dftk_band_filter_data;

typedef struct {		// ok
	dftk_process_header	hdr;
	int32_t				user_space_size;
	char				reserved[32];
} dftk_user_process_data;

typedef struct {						// compiler says sizeof(Spi_hs_header) = 136 bytes
	uint32_t		version;
	int32_t			timestamp;
	int32_t			checksum;
	int32_t			num_dg_points;		// number of digitized headshape points FOLLOWING THIS HEADER
	Pnt				ip_pnts[5];			// Cartesian coordinates of digitized index points ( 120 bytes )
} Spi_hs_header;
