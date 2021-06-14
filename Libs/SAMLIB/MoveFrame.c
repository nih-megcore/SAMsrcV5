// MoveFrame() -- Move an array of all sensors (primary & reference)
//	to a new frame.
//
//	This version uses a 6-parameter transform: translation in x, y, z
//	& rotation using Euler angle (in unique order)
//
//      Author: Stephen E. Robinson
//              MEG Core Facility
//              NIMH
//

#include <geoms.h>
#include <samlib.h>
#include <siglib.h>


void	MoveFrame(
	HeaderInfo	*Header,			// dataset header info
	ChannelInfo	*InChannel,			// input channel array
	ChannelInfo	*OutChannel,		// output (transformed) channel array
	char		*TransFile			// full path name for transform file
)

{
	XFORM		Where;				// where to move all SQUID sensors
	double		Vr[3];				// rotation around y-axis (radians)
	int			*SQUIDIndex;		// SQUIDIndex[S] -- index of channels for rotation/translation
	int			S;					// total SQUID sensors
	int			i;					// channel index
	int			j;					// 'nuther index
	FILE		*fp;				// generic file pointer

	// determine S
	for(i=S=0; i<Header->NumChannels; i++)
		if(InChannel[i].ChannelType == TYPE_MEG
		|| InChannel[i].ChannelType == TYPE_REF_GRAD
		|| InChannel[i].ChannelType == TYPE_REF_MAG)
			S++;

	// create SQUID index
	if((SQUIDIndex = (int *)malloc((size_t)S * sizeof(int))) == NULL)
		allocfailed("SQUIDIndex[]");
	for(i=j=0; i<Header->NumChannels; i++)
		if(InChannel[i].ChannelType == TYPE_MEG
		|| InChannel[i].ChannelType == TYPE_REF_GRAD
		|| InChannel[i].ChannelType == TYPE_REF_MAG) {
			SQUIDIndex[j] = i;
			j++;
		}

	// move all MEG primary and reference sensors
	if((fp = fopen(TransFile, "r")) == NULL)
		cleanup("can't open transform");
	if(fscanf(fp, "%lf%lf%lf%lf%lf%lf", &Where.Translate[X_], &Where.Translate[Y_], &Where.Translate[Z_], &Vr[PSI_], &Vr[THETA_], &Vr[PHI_]) != 6)
		cleanup("can't read transform");
	fclose(fp);
	EtoR(Vr, Where.Rotate);
	MoveMEG(InChannel, OutChannel, SQUIDIndex, S, &Where);
}

