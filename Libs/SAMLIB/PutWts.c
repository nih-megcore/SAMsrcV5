// PutWts() -- write legacy weight to SAM subdirectory for compatibility with DataEditor
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <byteswap.h>
#include <samlib.h>
#include <DataFiles.h>
#include <SAMfiles.h>

#define RETRO_VERSION	1		// this is a legacy weight file

void	PutWts(
	char		*WgtPath,       // full weight file path
	gsl_matrix	*Wgts,          // Wgts - VxM beamformer coefficients
	HeaderInfo	Header,         // SAM header information
	COV_HDR		CovHdr          // covariance header (input)
)

{
	SAM_HDR		WgtHdr;			// weight header (output)
	double		tmp;			// single weight element
	double		Dtmp;			// temporary double for flipping endian
	int			m;				// channel index
	int			n;				// alternate channel index
	int			v;				// vector/voxel index
	int			M;				// number of primary channels
	int			V;				// number of voxels
	char		CoeffID[] =		// SAM image file ID
		"SAMCOEFF";
	static FILE	*fp;			// output file file pointer

	// fetch dimensions of Wgts
	M = Header.NumPri;
	V = Wgts->size1;
	if(M != Wgts->size2)
		cleanup("PutWts: number of channels in Header differs from that in Wgt");

	// fill in header
	memset((void *)&WgtHdr, 0, sizeof(SAM_HDR));	// clear entire header
	WgtHdr.Version = bswap_32(RETRO_VERSION);		// file version number
	sprintf(WgtHdr.SetName, "%s", Header.SetName);	// name of parent dataset
	WgtHdr.NumChans = bswap_32(M);					// number of channels used by SAM
	WgtHdr.NumWeights = bswap_32(V);				// number of SAM virtual sensors
	Dtmp = 0.;
	WgtHdr.XStart = Dflip(&Dtmp);					// x-start coordinate (m)
	Dtmp = 0.;
	WgtHdr.XEnd = Dflip(&Dtmp);						// x-end coordinate (m)
	Dtmp = 0.;
	WgtHdr.YStart = Dflip(&Dtmp);					// y-start coordinate (m)
	Dtmp = 0.;
	WgtHdr.YEnd = Dflip(&Dtmp);						// y-end coordinate (m)
	Dtmp = 0.;
	WgtHdr.ZStart = Dflip(&Dtmp);					// z-start coordinate (m)
	Dtmp = 0.;
	WgtHdr.ZEnd = Dflip(&Dtmp);						// z-end coordinate (m)
	Dtmp = 0.;
	WgtHdr.StepSize = Dflip(&Dtmp);					// voxel step size (m)
	Dtmp = CovHdr.HPFreq;
	WgtHdr.HPFreq = Dflip(&Dtmp);					// highpass frequency (Hz)
	Dtmp = CovHdr.LPFreq;
	WgtHdr.LPFreq = Dflip(&Dtmp);					// lowpass frequency (Hz)
	Dtmp = CovHdr.BWFreq;
	WgtHdr.BWFreq = Dflip(&Dtmp);					// bandwidth of filters (Hz)
	Dtmp = CovHdr.Noise;
	WgtHdr.MeanNoise = Dflip(&Dtmp);				// mean primary sensor noise (T^2)
	memcpy((void *)WgtHdr.MriName, "NONE", 4);		// a place-holder until we decide whether to read the default.hdm file
	for(v=0; v<3; v++) {
		WgtHdr.Nasion[v] = bswap_32(0);
		WgtHdr.RightPA[v] = bswap_32(0);
		WgtHdr.LeftPA[v] = bswap_32(0);
	}
	WgtHdr.SAMType = bswap_32(SAM_TYPE_WT_LIST);
	WgtHdr.SAMUnit = bswap_32(SAM_UNIT_COEFF);

	// open legacy weight file & write ID, header, and channel index
	if((fp = fopen(WgtPath, "w")) == NULL)
		cleanup("can't open weight file for write");
	if(fwrite((void *)CoeffID, 8, 1, fp) != 1)
		cleanup("can't write ID to weight file");
	if(fwrite((void *)&WgtHdr, sizeof(SAM_HDR), 1, fp) != 1)
		cleanup("can't write header to weight file");
	for(m=0; m<M; m++) {
		n = bswap_32(Header.PsIndex[m]);
		if(fwrite((void *)&n, sizeof(int), 1, fp) != 1)
			cleanup("can't write channel index to weight file");
	}

	// write weight array until file closed
	for(v=0; v<V; v++)
		for(m=0; m<M; m++) {
			Dtmp = gsl_matrix_get(Wgts, v, m);
			tmp = Dflip(&Dtmp);
			if(fwrite((void *)&tmp, sizeof(double), 1, fp) != 1)
				cleanup("can't write weight to file");
		}
	fclose(fp);
}
