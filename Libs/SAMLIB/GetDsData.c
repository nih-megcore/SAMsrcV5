// GetDsData() -- read one channel/epoch from CTF .meg4 data file
//
//      Author: SE Robinson
//              MEG Core Facility
//              NIMH
//

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <math.h>
#include <byteswap.h>
#include <string.h>
#include <samlib.h>
#include <DataFiles.h>

FILE *OpenDsFile(
	HeaderInfo	*Header,
	char		*ext
)
{
	int			i;
	static char	*buf;
	FILE		*f;

	i = strlen(Header->DsPath) + strlen(Header->SetName) * 2 + strlen(ext) + 7;
	if((buf = (char *)malloc((size_t)i)) == NULL)
		allocfailed("buf");
	sprintf(buf, "%s/%s.ds/%s.%s", Header->DsPath, Header->SetName, Header->SetName, ext);
	f = fopen(buf, "r");
	free(buf);
	return f;
}

static FILE	*OpenMeg4File(
	HeaderInfo	*Header,
	int			spill
)
{
	int			i;
	static char	*buf;
	FILE		*f;

	i = strlen(Header->DsPath) + strlen(Header->SetName) * 2 + 20; // + 13 or 14 at least
	if((buf = (char *)malloc((size_t)i)) == NULL)
		allocfailed("buf");
    if(spill > 0) {
		sprintf(buf, "%s/%s.ds/%s.%d_meg4", Header->DsPath, Header->SetName, Header->SetName, spill);
	} else {
		sprintf(buf, "%s/%s.ds/%s.meg4", Header->DsPath, Header->SetName, Header->SetName);
	}
	f = fopen(buf, "r");
	free(buf);
	return f;
}


void	GetDsData(
	HeaderInfo	*Header,		// header information
	int			e,				// epoch index
	int			m,				// channel index
	double		scale,			// integer->Tesla conversion factor
	double		*Signal			// Signal[T] -- one channel/epoch of time series
)
{
	int32_t			i;			// sample value
	int				t;			// sample index
	off_t			offset;		// number of bytes offset
	static int		spill;		// 0 for .meg4, 1 for .1_meg4, etc.
	static int32_t	*xint;		// xint[T] -- raw time series buffer
	FILE			*meg4File;	// data file pointer
	static int		M = 0;		// number of channels
	static int		T;			// number of samples per trial

	// initialize
	if(M == 0) {
		M = Header->NumChannels;
		T = Header->MaxSamples;
		if((xint = (int32_t *)malloc((size_t)T * sizeof(int32_t))) == NULL)
			allocfailed("xint");
	}

	// initial conditions
	spill = 0;
	if(Header->spillt) {
		spill = e / Header->spillt;
		e = e % Header->spillt;
	}

	// open file
	if((meg4File = OpenMeg4File(Header, spill)) == NULL)
        cleanup("'meg4File' open failed");

	// offset data to epoch & read data
	if(Signal != NULL) {
		offset = (off_t)((T * 4) * (m + e * M) + 8);
		if(fseek(meg4File, offset, SEEK_SET) != 0)
			cleanup("'meg4File' fseek failed");
		if(fread((void *)xint, sizeof(int32_t), T, meg4File) != T)
			cleanup("'meg4File' read failed");
		for(t=0; t<T; t++) {
			i = bswap_32(xint[t]);
			Signal[t] = scale * (double)i;
		}
	}

	// done
	fclose(meg4File);
}
