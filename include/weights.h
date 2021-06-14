#ifndef	WEIGHTS_INCLUDED
#define	WEIGHTS_INCLUDED
/*****************************************************************

	$Id: weights.h,v 1.1 1996/07/25 15:57:18 skip Exp $

	Name:		weights.h

	Definitions for computing and applying weights.

	Author:		Ken Velarde

	Creation Date:	12 July 1996


	Copyright 1996, Biomagnetic Technologies, inc.

	$Log: weights.h,v $
# Revision 1.1  1996/07/25  15:57:18  skip
# Initial revision
#


*****************************************************************/


// Define a structure to hold needed statistics for each channel
struct ChanStatistics
{
	// These are for channel identification and scaling
	char	name [DFTK_CHANNEL_STRLEN];	// Channel name
	short	type;				// Channel type
	short	chan_no;			// Channel number
	long	in_index;			// Index into the input data file
	long	out_index;			// Index into the output data buffer
	long	file_index;			// Index into the output data file
	double	conversion;			// Factor to convert stored data
						// into real values
	// These hold the DC offset, sum and sum squared
	long	num_cov;			// Number of covariance elements
	double	sum;				// Sum for the epoch or part
	double	sum2;				// Sum of squares for the epoch or part
	double	*cov;				// Covariance elements

	// These hold the variance & covariance
	double	Variance;			// Variance
	double	dc_offset;			// Variance
	double	*Covariance;			// Covariance
	double	*Weights;			// Weight values
};

// Define a structure to hold needed statistics for each channel
struct afw_Weights
{
	// Dimensions
	long	num_signal_channels;
	long	num_reference_channels;
	char	Creator[PDF_NAME_LEN];
	float	HP_Filter;
	time_t  timestamp;

	ChanStatistics	*reference_channels;	// Array of reference channels
	ChanStatistics	*signal_channels;	// Array of signal channels
};

int WriteWeights(afw_Weights *out_weights, char *filename, int version);
ChanStatistics  *FindChannel(ChanStatistics *channels, int n, char *name);
ChanStatistics  *GetChannelArray (dftk_pdf *pdf, ChanStatistics *input_channels, int total_channels,
			dftk_channel_ref **chan_list, short type, long *numchans);
ChanStatistics *GetChannelsToCopy (dftk_pdf *pdf, dftk_channel_ref **chan_list, int *numchans);
afw_Weights *ReadWeights (char *filename);
float GetHP_Filter (dftk_pdf *input_pdf);

#endif
