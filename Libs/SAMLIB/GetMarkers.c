// GetMarkers() -- read CTF 'MarkerFile.mrk' file, including modifications
//	by A Moissiev for safeguarding against memory overwrite by marker names
//	longer than 15 characters
//
//	Author:	Stephen E. Robinson
//			MEG Core Group
//			NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <TimeSegs.h>
#include <lfulib.h>

#define TRUE	1
#define FALSE	0

void	GetMarkers(
	char		*MarkerFile,
	SAM_MARKS	**MarkerOut,
	int			*NumMarkers
)

{
	static SAM_MARKS	*Marker;		// Marker[N]
	double				a;
	double				b;
	int					NumClasses;		// number of marker classes
	int					NumSamples;		// number of markers within a class
	int					i;				// class-index
	int					j;				// sample-endex
	int					k;				// marker-index
	int					done;
	const int			MRK_MAX_LEN = sizeof(Marker->Name) - 1;
	char				LabelField[64];
	char				MarkerName[64];
	FILE				*fp;

	// open 'MarkerFile.mrk' & read parameters req'd to allocate 'Marker' structure
	if((fp = fopen(MarkerFile, "r")) == NULL)
		cleanup("can't open 'MarkerFile.mrk' for read");

	do {
		if(fgets(LabelField, 64, fp) == NULL)
			cleanup("can't read 'LabelField' string");
	} while(strncmp(LabelField, "NUMBER OF MARKERS:", 18) != 0);
	if(fscanf(fp, "%d", &NumClasses) == EOF)
		cleanup("can't read 'NumClasses'");
	if(NumClasses < 1)
		cleanup("'NumClasses < 1");
	for(i=*NumMarkers=0; i<NumClasses; i++) {
		do {
			if(fgets(LabelField, 64, fp) == NULL)
				cleanup("can't read 'LabelField' string");
		} while(strncmp(LabelField, "NUMBER OF SAMPLES:", 18) != 0);
		if(fscanf(fp, "%d", &NumSamples) == EOF)
			cleanup("can't read 'NumSamples'");
		*NumMarkers += NumSamples;
	}
	fclose(fp);
	if((Marker = (SAM_MARKS *)malloc((size_t)(*NumMarkers * sizeof(SAM_MARKS)))) == NULL)
		allocfailed("Marker[]");
	*MarkerOut = Marker;

	// reopen 'MarkerFile.mrk' & read parameters into 'Marker' structure
	if((fp = fopen(MarkerFile, "r")) == NULL)
		cleanup("can't reopen 'MarkerFile.mrk' for read");

	for(i=k=0; i<NumClasses; i++) {

		do {	// search for class name
			if(fgets(LabelField, 64, fp) == NULL)
				cleanup("can't read 'LabelField'");
		} while(strncmp(LabelField, "NAME:", 5) != 0);
		if(fscanf(fp, "%s", MarkerName) == EOF)
			cleanup("can't read 'MarkerName'");
		if(strlen(MarkerName) > MRK_MAX_LEN) {
			fprintf(stderr, "\nWARNING: Name of the marker '%s' is too long. Truncated to '%s'\n", MarkerName, Marker[k].Name);
			fflush(stderr);
		}

		do {	// search for number of samples in this class
			if(fgets(LabelField, 64, fp) == NULL)
				cleanup("can't read 'LabelField'");
		} while(strncmp(LabelField, "NUMBER OF SAMPLES:", 18) != 0);
		if(fscanf(fp, "%d", &NumSamples) == EOF)
			cleanup("can't read 'NumSamples'");

		do {	// synch up to line before marker list
			if(fgets(LabelField, 64, fp) == NULL)
				cleanup("can't read 'LabelField'");
		} while(strncmp(LabelField, "TRIAL NUMBER", 12) != 0);

		for(j=0; j<NumSamples; j++, k++) {
			memset(&(Marker[k].Name), 0, sizeof(Marker[k].Name));	// clear marker name field for current marker
			strcpy(Marker[k].Name, MarkerName);
			if(fscanf(fp, "%d%lf", &Marker[k].Epoch, &Marker[k].Latency) == EOF)
				cleanup("can't read marker epoch & latency");
		}
	}
	fclose(fp);

	// sort the markers in time-order
	do {
		for(i=1, done=TRUE; i<*NumMarkers; i++) {
			a = 1000. * Marker[i-1].Epoch + Marker[i-1].Latency;
			b = 1000. * Marker[i].Epoch + Marker[i].Latency;
			if(a > b) {	// classic bubble-sort
				j = Marker[i-1].Epoch;
				Marker[i-1].Epoch = Marker[i].Epoch;
				Marker[i].Epoch = j;
				a = Marker[i-1].Latency;
				Marker[i-1].Latency = Marker[i].Latency;
				Marker[i].Latency = a;
				strcpy(MarkerName, Marker[i-1].Name);
				strcpy(Marker[i-1].Name, Marker[i].Name);
				strcpy(Marker[i].Name, MarkerName);
				done = FALSE;
			}
		}
	} while(done == FALSE);
}
