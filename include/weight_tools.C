/*****************************************************************

    $Id: weight_tools.C,v 1.3 2000/05/18 18:46:39 skip Exp $

    Name:        afw - Apply Fixed Weights

    Purpose:    This application reads a weight file and
            applies the weights to the input data file
            for the specified signal channels.

    Command Line Options:
            The following options are required.  Note
            That the data file must be either specified
            using the PSsrp options or using the -posted
            option.
            -P patient:
                The Patient ID to use for the
                input data file
            -S scan:
                The Scan Name to use for the
                input data file.
            -s session:
                The Session Name to use for the
                input data file.  There must
                be one instance of this option
                for each data file.
            -r run:
                The Run Name to use for the
                input data file.  There must
                be one instance of this option
                for each data file.
            -p datafile:
                The Data File Name to use for the
                input data file.  There must
                be one instance of this option
                for each data file.
            -posted:
                Flag indicating that the input
                data file is to be taken from
                the first posted data file.
            -w weightfile:
                Name of the weights file to read.

            The following options are optional:
            -d:    Debugging flag
            -O:    Option to remove DC offset
            -n prefix:
                The prefix to prepend to the
                output data file.
            -T u<channel-list>:
                List of signal channels which are
                to be used.  Default is to use all
                signal channels.
            -T c<channel-list>:
                List of channels to copy.  This is
                the list of non-signal and non-reference
                channels which is copied to the output
                data file.  Default is to copy all channels.

    Returns:    0    No error
            1    Command line error
            2    Error reading input data file
            3    Error reading input weight file
            4    Error writing output data file

    Inputs:        int BTi_debug:
                Software debugging flag.  When
                true causes detailed debugging
                statements to be printed to the
                standard error stream.
            char *gProgram:
                Pointer to the name of the
                executable.

    Outputs:    Writes the noise reduced data file into
            the same run as the input file and appends
            the ",n" suffix.

    Author:        Ken Velarde

    Creation Date:    18 July 1996


    Copyright 1996, Biomagnetic Technologies, inc.

    $Log: weight_tools.C,v $
    Revision 1.3  2000/05/18 18:46:39  skip
    Include global header file and define weight_debug macro.

    Revision 1.2  2000/01/11 23:09:18  skip
    Updated to latest WHS version.

    Revision 1.1  1998/05/01 16:49:53  scobb
    Initial revision


*****************************************************************/

#ifdef HAVE_CONFIG_H
#    include <config.h>
#endif /* HAVE_CONFIG_H */

#include <dftk_pdf.h>        // for classes dftk_pdf and pdf_id
#include <clparser.h>        // for processing tools functions
#include <strings.h>
#include <weights.h>
#include "MSIparse.h"
#include <dftk_weights.h>

// Macro to turn on weight debugging
#define WEIGHT_DEBUG ((getenv("weight_debug") == NULL) ? 0 : atoi(getenv("weight_debug")))

// Global variables
extern int      BTi_debug;
extern short    gRemoveDC;
extern char    *gProgram;
extern char    *gWeightFile;       // Name of the weight file to read
extern char    *gReferenceName;    // Get Reference channels from this PDF

afw_Weights *ReadWeightsVersionOne (char *filename);
afw_Weights *ReadWeightsOrig (char *filename);

afw_Weights *
ReadWeights (char *filename)
{
    // Local Variables
    FILE        *fp;
    int        i, j;
    char        version_string[256];

    // Open the file
    if (WEIGHT_DEBUG) fprintf (stderr, "Attempting to open %s\n", filename);
    if ((fp = fopen (filename, "r")) == NULL)
    {    fprintf (stderr, "Error opening %s\n", filename);
    perror (filename);
    return (NULL);
    }

    fscanf (fp, "%s", version_string);
    if (!strcmp(version_string, "MSI.WeightTable.Version"))
    {
        // newer (non original weight table) figure out what version and run that
        // Close the file
        if (WEIGHT_DEBUG) printf("Newer version of weights...\n");
        fclose (fp);
        return (ReadWeightsVersionOne(filename));
    }
    else
    {
        // assume old (original version of weight table, and try to read it
        // Close the file
        if (WEIGHT_DEBUG) printf("Old version of weights...\n");
        fclose (fp);
        return (ReadWeightsOrig(filename));
    }
}

afw_Weights *
        ReadWeightsVersionOne (char *filename)
{
    // Local Variables
    FILE        *fp;
    int        i, j;
    ChanStatistics    *cptr;
    afw_Weights        *weights;

    // Open the file
    if (WEIGHT_DEBUG) fprintf (stderr, "Attempting to open %s\n", filename);
    if ((fp = fopen (filename, "r")) == NULL)
    {    fprintf (stderr, "Error opening %s\n", filename);
    perror (filename);
    return (NULL);
    }

    // Allocate the weights structur
    if ((weights = new afw_Weights) == NULL)
    {    fprintf (stderr, "Error allocating weights\n");
    return (NULL);
    }

    long    temp_long;

    if (!msi_file_get_long(fp, "MSI.NumSignalChans", &temp_long))
    {
        weights->num_signal_channels = temp_long;
        if (weights->num_signal_channels < 0) return(NULL);
    }
    else
    {
        if (WEIGHT_DEBUG)
            printf("Error: gould not parse TAG %s in file %s\n",
                   "MSI.NumSignalChans", filename);
        return(NULL);
    }

    if (!msi_file_get_long(fp, "MSI.NumReferenceChans", &temp_long))
    {
        weights->num_reference_channels = temp_long;
        if (weights->num_reference_channels < 0) return(NULL);
    }
    else
    {
        if (WEIGHT_DEBUG)
            printf("Error: gould not parse TAG %s in file %s\n",
                   "MSI.NumReferenceChans", filename);
        return(NULL);
    }

    if (!msi_file_get_long(fp, "MSI.TimeStamp", &temp_long))
    {
        weights->timestamp = temp_long;
    }
    else
    {
        if (WEIGHT_DEBUG)
            printf("Error: gould not parse TAG %s in file %s\n",
                   "MSI.TimeStamp", filename);
        weights->timestamp = (long)0;
    }

    sprintf(weights->Creator, "\0");    // for now leave blank

    if (WEIGHT_DEBUG) printf("%ld signal, %ld reference channels\n",
    weights->num_signal_channels, weights->num_reference_channels);


    // Write the number of signal and reference channels
    // Allocate the reference channels
    weights->reference_channels = new ChanStatistics [weights->num_reference_channels];
    if (!weights->reference_channels)
    {
        fprintf (stderr, "Error allocating reference channel array\n");
        delete weights;
        return (NULL);
    }

    // Allocate the signal channel array
    weights->signal_channels = new ChanStatistics [weights->num_signal_channels];
    if (!weights->signal_channels)
    {
        fprintf (stderr, "Error allocating signal channel array\n");
        delete [] weights->reference_channels;
        delete weights;
        return (NULL);
    }

    //  Skip to the actual tabel in the file
    if (msi_file_skipto (fp, "MSI.WeightTable.Data"))
    {
        printf("Error skipping to TAG %s in file %s\n", "MSI.WeightTable.Data", filename);
        return(NULL);
    }

    // Read in the reference channel names
    for (i = 0, cptr = weights->reference_channels; i < weights->num_reference_channels; i++, cptr++)
    {
        bzero ((char *)cptr, sizeof (ChanStatistics));
        fscanf (fp, "%s", cptr->name);
    }

    // Read in the signal channels
    for (i = 0, cptr = weights->signal_channels; i < weights->num_signal_channels; i++, cptr++)
    {
        // Read in the name
        bzero ((char *)cptr, sizeof (ChanStatistics));
        fscanf (fp, "%s", cptr->name);

        // Allocate the weights array
        if ((cptr->Weights = new double [weights->num_reference_channels]) == NULL)
        {    fprintf (stderr, "Error allocating weights for channel %s\n", cptr->name);
        delete weights;
        return (NULL);
        }

        // Read the weights
        for (j = 0; j < weights->num_reference_channels; j++)
            fscanf (fp, "%lf", &cptr->Weights [j]);

        if ((WEIGHT_DEBUG) && (i == 0))
        {
            printf("First Channel: %s =  %lf\n", cptr->name, cptr->Weights[0]);
        }

    }

    // Close the file
    fclose (fp);
    return (weights);
}


afw_Weights *
        ReadWeightsOrig (char *filename)
{
    // Local Variables
    FILE        *fp;
    int        i, j;
    ChanStatistics    *cptr;
    afw_Weights        *weights;

    // Open the file
    if (WEIGHT_DEBUG) fprintf (stderr, "Attempting to open %s\n", filename);
    if ((fp = fopen (filename, "r")) == NULL)
    {
        fprintf (stderr, "Error opening %s\n", filename);
        perror (filename);
        return (NULL);
    }

    // Allocate the weights structur
    if ((weights = new afw_Weights) == NULL)
    {
        fprintf (stderr, "Error allocating weights\n");
        return (NULL);
    }

    fscanf (fp, "%ld\n", &weights->num_signal_channels);
    fscanf (fp, "%ld\n", &weights->num_reference_channels);

    if (WEIGHT_DEBUG)
        printf("%ld signal, %ld reference channels\n",
        weights->num_signal_channels, weights->num_reference_channels);

    sprintf(weights->Creator, "\0");    // for now leave blank

    // Write the number of signal and reference channels
    // Allocate the reference channels
    weights->reference_channels = new ChanStatistics [weights->num_reference_channels];
    if (!weights->reference_channels)
    {
        fprintf (stderr, "Error allocating reference channel array\n");
        delete weights;
        return (NULL);
    }

    // Allocate the signal channel array
    weights->signal_channels = new ChanStatistics [weights->num_signal_channels];
    if (!weights->signal_channels)
    {
        fprintf (stderr, "Error allocating signal channel array\n");
        delete [] weights->reference_channels;
        delete weights;
        return (NULL);
    }

    // Read in the reference channel names
    for (i = 0, cptr = weights->reference_channels; i < weights->num_reference_channels; i++, cptr++)
    {    bzero ((char *)cptr, sizeof (ChanStatistics));
    fscanf (fp, "%s", cptr->name);
    }

    // Read in the signal channels
    for (i = 0, cptr = weights->signal_channels; i < weights->num_signal_channels; i++, cptr++)
    {
        // Read in the name
        bzero ((char *)cptr, sizeof (ChanStatistics));
        fscanf (fp, "%s", cptr->name);

        // Allocate the weights array
        if ((cptr->Weights = new double [weights->num_reference_channels]) == NULL)
        {    fprintf (stderr, "Error allocating weights for channel %s\n", cptr->name);
        delete weights;
        return (NULL);
        }

        // Read the weights
        for (j = 0; j < weights->num_reference_channels; j++)
            fscanf (fp, "%lf", &cptr->Weights [j]);
    }

    // Close the file
    fclose (fp);
    return (weights);
}

ChanStatistics *
        GetChannelsToCopyAndRef(dftk_pdf *pdf, dftk_channel_ref **chan_list, int *numchans,
                                ChanStatistics *ref_chan_stats, int num_ref_stats)
{
    // Local variables
    dftk_channel_ref    **chan_ref_ptr, *chan_ref;
    dftk_channel_ref    **ref_chan_ref_ptr, *ref_chan_ref;
    list_position    lp;
    ChanStatistics    *channels, *out_chan_ptr, *ref_ptr;
    short        format = pdf->get_data_format ();
    int i;

    // First allocate enough room for the channels
    (*numchans) = 0;
    if ((channels = new ChanStatistics [pdf->get_total_channels ()]) == NULL)
    {
	fprintf (stderr, "Error allocating channels array\n");
	return (NULL);
    }

    // If the chan_list is non-NULL, only get those channels
    out_chan_ptr = channels;
    if (chan_list)
    {
	// Build up the output array from all non-MEGTYPE and non-REFERENCE channels
	// in the input list.
        for (chan_ref_ptr = chan_list; (*chan_ref_ptr); chan_ref_ptr++)
        {
            // Only deal with the correct type of channel
            chan_ref = (*chan_ref_ptr);
            if (chan_ref->get_type () == MEGTYPE ||
                chan_ref->get_type () == TRIGGER) continue;

            if (chan_ref->get_type () == REFERENCE)
            {
        	// now find any ignored REFERENCE channels
        	// if this channel is NOT in the ref_channel list,
        	// also copy the channel
                int found_ref = FALSE;
                for (i = 0,  ref_ptr = ref_chan_stats; i < num_ref_stats; i++, ref_ptr++)
                {
                    if (!strcmp(chan_ref->get_name(), ref_ptr->name))
                    {
                        found_ref = TRUE;
                    }
                }
                if (found_ref)
                {
                    continue;
                }
                else
                {
            	    // didn't find ref, therefore need to copy
                    if (WEIGHT_DEBUG) fprintf(stderr, "copying ref: %s\n", chan_ref->get_name());
                }
            }

            // Set the information about this channel
            bzero ((char *)out_chan_ptr, sizeof (ChanStatistics));
            strcpy (out_chan_ptr->name, chan_ref->get_name ());
            out_chan_ptr->type       = chan_ref->get_type ();
            out_chan_ptr->chan_no    = chan_ref->get_chan_no ();
            out_chan_ptr->in_index   = chan_ref->get_index ();
            out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
            out_chan_ptr++;
            (*numchans)++;
        }
    }
    else    // Get all non-MEGTYPE and non-REFERENCE channels
    {
	for (lp = pdf->first_channel_ref_position (); lp; lp = pdf->next_channel_ref_position (lp))
	{
            // Get the channel ref from the pdf
            chan_ref = pdf->fetch_channel_ref (lp);

            if (chan_ref->get_type () == MEGTYPE ||
                chan_ref->get_type () == TRIGGER) continue;

            if (chan_ref->get_type () == REFERENCE)
            {
                // now find any ignored REFERENCE channels
                // if this channel is NOT in the ref_channel list,
                // also copy the channel
                int found_ref = FALSE;
                for (i = 0,  ref_ptr = ref_chan_stats; i < num_ref_stats; i++, ref_ptr++)
                {
                    if (!strcmp(chan_ref->get_name(), ref_ptr->name))
                    {
                        found_ref = TRUE;
                    }
                }

                if (found_ref)
                {
                    continue;
                }
                else
                {
                    // didn't find ref, therefore need to copy
                    if (WEIGHT_DEBUG) fprintf(stderr, "copying ref: %s\n", chan_ref->get_name());
                }
            }

            // Set the information about this channel
            bzero ((char *)out_chan_ptr, sizeof (ChanStatistics));
            strcpy (out_chan_ptr->name, chan_ref->get_name ());
            out_chan_ptr->type       = chan_ref->get_type ();
            out_chan_ptr->chan_no    = chan_ref->get_chan_no ();
            out_chan_ptr->in_index   = chan_ref->get_index ();
            out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
            out_chan_ptr++;
            (*numchans)++;
        }
    }

    // Always copy the trigger channels
    for (lp = pdf->first_channel_ref_position (); lp; lp = pdf->next_channel_ref_position (lp))
    {
        // Get the channel ref from the pdf
        chan_ref = pdf->fetch_channel_ref (lp);

        if (chan_ref->get_type () == TRIGGER)
        {
            // Set the information about this channel
            bzero ((char *)out_chan_ptr, sizeof (ChanStatistics));
            strcpy (out_chan_ptr->name, chan_ref->get_name ());
            out_chan_ptr->type       = chan_ref->get_type ();
            out_chan_ptr->chan_no    = chan_ref->get_chan_no ();
            out_chan_ptr->in_index   = chan_ref->get_index ();
            out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
            out_chan_ptr++;
            (*numchans)++;
        }
    }

    return (channels);
}

ChanStatistics *
        GetChannelsToCopy (dftk_pdf *pdf, dftk_channel_ref **chan_list, int *numchans)
{
    // Local variables
    dftk_channel_ref    **chan_ref_ptr, *chan_ref;
    list_position        lp;
    ChanStatistics        *channels, *out_chan_ptr;
    short            format = pdf->get_data_format ();

    // First allocate enough room for the channels
    (*numchans) = 0;
    if ((channels = new ChanStatistics [pdf->get_total_channels ()]) == NULL)
    {    fprintf (stderr, "Error allocating channels array\n");
    return (NULL);
    }

    // If the chan_list is non-NULL, only get those channels
    out_chan_ptr = channels;
    if (chan_list)
    {
        // Build up the output array from all non-MEGTYPE and non-REFERENCE channels
        // in the input list.
        for (chan_ref_ptr = chan_list; (*chan_ref_ptr); chan_ref_ptr++)
        {
            // Only deal with the correct type of channel
            chan_ref = (*chan_ref_ptr);
            if (chan_ref->get_type () == MEGTYPE ||
                chan_ref->get_type () == REFERENCE ||
                chan_ref->get_type () == TRIGGER) continue;

            // Set the information about this channel
            bzero ((char *)out_chan_ptr, sizeof (ChanStatistics));
            strcpy (out_chan_ptr->name, chan_ref->get_name ());
            out_chan_ptr->type       = chan_ref->get_type ();
            out_chan_ptr->chan_no    = chan_ref->get_chan_no ();
            out_chan_ptr->in_index         = chan_ref->get_index ();
            out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
            out_chan_ptr++;
            (*numchans)++;
        }
    }
    else    // Get all non-MEGTYPE and non-REFERENCE channels
    {
        for (lp = pdf->first_channel_ref_position (); lp; lp = pdf->next_channel_ref_position (lp))
        {
            // Get the channel ref from the pdf
            chan_ref = pdf->fetch_channel_ref (lp);

            if (chan_ref->get_type () == MEGTYPE ||
                chan_ref->get_type () == REFERENCE ||
                chan_ref->get_type () == TRIGGER) continue;

            // Set the information about this channel
            bzero ((char *)out_chan_ptr, sizeof (ChanStatistics));
            strcpy (out_chan_ptr->name, chan_ref->get_name ());
            out_chan_ptr->type       = chan_ref->get_type ();
            out_chan_ptr->chan_no    = chan_ref->get_chan_no ();
            out_chan_ptr->in_index   = chan_ref->get_index ();
            out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
            out_chan_ptr++;
            (*numchans)++;
        }
    }


    // Always copy the trigger channels
    for (lp = pdf->first_channel_ref_position (); lp; lp = pdf->next_channel_ref_position (lp))
    {
        // Get the channel ref from the pdf
        chan_ref = pdf->fetch_channel_ref (lp);

        if (chan_ref->get_type () == TRIGGER)
        {
            // Set the information about this channel
            bzero ((char *)out_chan_ptr, sizeof (ChanStatistics));
            strcpy (out_chan_ptr->name, chan_ref->get_name ());
            out_chan_ptr->type       = chan_ref->get_type ();
            out_chan_ptr->chan_no    = chan_ref->get_chan_no ();
            out_chan_ptr->in_index   = chan_ref->get_index ();
            out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
            out_chan_ptr++;
            (*numchans)++;
        }
    }

    return (channels);
}

ChanStatistics *
        GetRefChannelsToCopy (dftk_pdf *pdf, int *numchans)
{
    // Local variables
    dftk_channel_ref    **chan_ref_ptr, *chan_ref;
    list_position        lp;
    ChanStatistics        *channels, *out_chan_ptr;
    short            format = pdf->get_data_format ();

    // First allocate enough room for the channels
    (*numchans) = 0;
    if ((channels = new ChanStatistics [pdf->get_total_channels ()]) == NULL)
    {
        fprintf (stderr, "Error allocating channels array\n");
        return (NULL);
    }

    // If the chan_list is non-NULL, only get those channels
    out_chan_ptr = channels;

    // Now get only REFERENCE channels
    for (lp = pdf->first_channel_ref_position (); lp; lp = pdf->next_channel_ref_position (lp))
    {
        // Get the channel ref from the pdf
        chan_ref = pdf->fetch_channel_ref (lp);

        if (chan_ref->get_type() != REFERENCE) continue;

        // Set the information about this channel
        bzero ((char *)out_chan_ptr, sizeof (ChanStatistics));
        strcpy (out_chan_ptr->name, chan_ref->get_name ());
        out_chan_ptr->type       = chan_ref->get_type ();
        out_chan_ptr->chan_no    = chan_ref->get_chan_no ();
        out_chan_ptr->in_index         = chan_ref->get_index ();
        out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
        out_chan_ptr++;
        (*numchans)++;
    }

    return (channels);
}

ChanStatistics  *
        GetChannelArray (dftk_pdf *pdf, ChanStatistics *input_channels, int total_channels,
                         dftk_channel_ref **chan_list, short type, long *numchans)
{
    // Function prototypes
    ChanStatistics        *FindChannel (ChanStatistics *list, int n, char *name);

    // Local variables
    ChanStatistics        *channels, *in_chan_ptr,
    *out_chan_ptr;
    list_position        lp;
    dftk_channel_ref    *chan_ref, **chan_ref_ptr;
    int            i;
    short            format = pdf->get_data_format ();

    // First allocate enough room for the channels
    (*numchans) = 0;
    if ((channels = new ChanStatistics [pdf->get_total_channels ()]) == NULL)
    {
        fprintf (stderr, "Error allocating channels array\n");
        return (NULL);
    }

    // If the chan_list is non-NULL, only get those channels
    out_chan_ptr = channels;
    if (chan_list)
    {
        // Build up the output array from both the input array and the channel list
        for (chan_ref_ptr = chan_list; (*chan_ref_ptr); chan_ref_ptr++)
        {
            // Only deal with the correct type of channel
            if ((type != ANY_CHAN_TYPE) && ((*chan_ref_ptr)->get_type() != type)) continue;

            // Find this channel in the input array
            in_chan_ptr = FindChannel (input_channels, total_channels, (*chan_ref_ptr)->get_name ());
            if (!in_chan_ptr)
            {
                if ((*chan_ref_ptr)->get_type() == MEGTYPE)
                {
                    fprintf (stderr, "Warning: channel %s not in weight file\n",
                             (*chan_ref_ptr)->get_name ());
                    continue;
                }
                else
                {
                    fprintf (stderr, "Error: channel %s not in weight file\n",
                             (*chan_ref_ptr)->get_name ());
                    fprintf (stderr, "Data file must contain all REFERENCE channels\n");
                    fprintf (stderr, "used to create the weight file.\n");
                    delete [] channels;
                    return (NULL);
                }
            }

            (*out_chan_ptr) = (*in_chan_ptr);
            out_chan_ptr->in_index = (*chan_ref_ptr)->get_index ();
            if (format == SHORT || format == LONG)
                out_chan_ptr->conversion = (*chan_ref_ptr)->get_conversion_factor ();
            else    out_chan_ptr->conversion = 1.0;
            out_chan_ptr++;
            (*numchans)++;
        }
    }
    else    // Only get the channels from the input_channels array
    {
        for (i = 0, in_chan_ptr = input_channels;
             i < total_channels; in_chan_ptr++, i++)
        {
                // Get the channel ref from the pdf
            if ((lp = pdf->search_channel_ref (in_chan_ptr->name)) == NULL)
            {
                if (type == MEGTYPE)
                {
                    fprintf (stderr,
                             "Warning: channel %s not in data file\n",
                             in_chan_ptr->name);
                    continue;
                }
                else
                {
                    fprintf (stderr, "Error: channel %s not in weight file\n",
                            in_chan_ptr->name);
                    fprintf (stderr,
                             "Data file must contain all REFERENCE channels\n");
                    fprintf (stderr,
                             "used to create the weight file.\n");
                    delete [] channels;
                    return (NULL);
                }
            }

                // Get a pointer to this channel
            chan_ref = pdf->fetch_channel_ref (lp);

            // Copy the information to the returned array
            (*out_chan_ptr) = (*in_chan_ptr);
            out_chan_ptr->in_index = chan_ref->get_index ();
            if (format == SHORT || format == LONG)
                out_chan_ptr->conversion = chan_ref->get_conversion_factor ();
            else    out_chan_ptr->conversion = 1.0;
            out_chan_ptr++;
            (*numchans)++;
        }
    }

    return (channels);
}

// finds the channel named 'name' in the array of channels
ChanStatistics  *
        FindChannel (ChanStatistics *channels, int n, char *name)
{
    // Local variables
    ChanStatistics    *cptr;
    int        i;

    for (i = 0, cptr = channels; i < n; i++, cptr++)
        if (strcmp (cptr->name, name) == 0) return (cptr);

    return (NULL);
}

int
        WriteWeights (afw_Weights *out_weights, char *filename, int version)
{
    // Local Variables
    FILE        *fp;
    int        i, j;
    ChanStatistics    *cptr;

    // Open the file
    if (WEIGHT_DEBUG) fprintf (stderr, "Attempting to open %s\n", filename);
    if ((fp = fopen (filename, "w+")) == NULL)
    {    fprintf (stderr, "Error opening %s\n", filename);
    perror (filename);
    return (1);
    }

    if (version)
    {
        fprintf(fp, "MSI.WeightTable.Version\t%d\t\n", version);
    }

    if (!version)
    {
        fprintf(fp, "%ld\n", out_weights->num_signal_channels);
        fprintf(fp, "%ld\n", out_weights->num_reference_channels);
    }
    else
    {
        fprintf(fp, "MSI.NumSignalChans\t%ld\t\n", out_weights->num_signal_channels);
        fprintf(fp, "MSI.NumReferenceChans\t%ld\t\n", out_weights->num_reference_channels);
        if (strlen(out_weights->Creator) > 1) fprintf(fp, "MSI.Creator\t%s\t\n", out_weights->Creator);
        if(out_weights->HP_Filter >= 0.0) fprintf(fp, "MSI.HPFilter\t%3.1f\t\n", out_weights->HP_Filter);
        fprintf(fp, "MSI.TimeStamp\t%ld\t\n", out_weights->timestamp);
        fprintf(fp, "MSI.WeightTable.Data\n");
    }

    // Write the reference channel names
    for (i = 0, cptr = out_weights->reference_channels; i < out_weights->num_reference_channels; i++, cptr++)
    {
        fprintf (fp, "\t%s\t", cptr->name);
    }

    fprintf (fp, "\t\n");

    // Write the signal channels
    for (i = 0, cptr = out_weights->signal_channels; i < out_weights->num_signal_channels; i++, cptr++)
    {
        // Write the name
        fprintf(fp, "%s\t", cptr->name);

        // Write the weights
        for (j = 0; j < out_weights->num_reference_channels; j++)
            fprintf(fp, "%3.6lf\t", cptr->Weights [j]);

        fprintf(fp, "\n");
    }

    // Close the file
    fclose (fp);
    return (0);
}

int get_applied_weights()
{
    return(0);
}

float
GetHP_Filter (dftk_pdf *input_pdf)
{
    char      *ETableSearchName, WeightTableBlockName[20];
    int        retVal;
    int        foundETable = FALSE;
    int        foundWeightTable = FALSE;
    list_position    EPos = NULL;
    list_position    FPos = NULL;
    dftk_user_block    *EBlock;
    dftk_user_block    *FBlock;
    time_t        timeStamp;
    dftk_config    *config;
    char        *EFilterName;
    char        *WeightName = APPLIED_WT_TABLE_NAME;

    config = (dftk_config *)input_pdf->get_config();

    long num_procs = input_pdf->get_total_procs();

    // first find acquisition High pass filter setting
    double        filt_freq = -0.9;
    EFilterName = "-0.9";    // if no filter set

    if ((num_procs))
    {
        if (WEIGHT_DEBUG) fprintf(stderr, "got %ld config procs\n", num_procs);
        list_position proc_lp = input_pdf->first_proc_position();
        while(num_procs--)
        {
            dftk_proc *proc = input_pdf->fetch_proc(proc_lp);
            if (proc)
            {
                if (!strcmp(proc->get_type(), "B_hw_filt"))
                {
                    if (WEIGHT_DEBUG) fprintf(stderr, "found good B_hw_filt\n");
                    long num_steps;
                    list_position step_lp = proc->first_step_position();
                    int found_filt = FALSE;
                    if(num_steps = proc->get_total_steps() > 0)
                    {
                        while((num_steps--) && (!found_filt))
                        {
                            dftk_gen_process *gen_proc = proc->fetch(step_lp);
                            if (!strcmp(gen_proc->get_type(), DFTK_FILT_HP))
                            {
                                found_filt = TRUE;
                                filt_freq = ((dftk_filter *)gen_proc)->get_frequency();
                                if (WEIGHT_DEBUG) fprintf(stderr, "Found HP filter %g\n", filt_freq);
                                double filt_delta = 0.001;    // to use for if tests below
                                if((filt_freq <= (0.0 + filt_delta)) && (filt_freq >= (0.0 - filt_delta)))
                                    EFilterName = "DC";
                                else if ((filt_freq <= (1.0 + filt_delta)) && (filt_freq >= (1.0 - filt_delta)))
                                    EFilterName = "1.0";
                                else if ((filt_freq <= (0.1 + filt_delta)) && (filt_freq >= (0.1 - filt_delta)))
                                    EFilterName = "0.1";
                                else
                                    EFilterName = "-0.7";  // AN ERROR CONDITION THAT SHOULD NEVER HAPPEN
                            }
                        }

                    }
                    else
                        if (WEIGHT_DEBUG) fprintf(stderr, "No steps\n");
                }
                else
                    if (WEIGHT_DEBUG) fprintf(stderr, "found %s\n", proc->get_type());
            }
            else
            {
                fprintf(stderr, "Bad Proc\n");
            }
            proc_lp = input_pdf->next_proc_position(proc_lp);
        }
    }
    else
    {
        fprintf(stderr, "got no config procs\n");
    }
    return(filt_freq);
}

