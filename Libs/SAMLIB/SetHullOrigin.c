// SetHullOrigins() -- sets Hull origin from *.hdm file
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>		// get axes definitions
#include <samlib.h>
#include <lfulib.h>
#include <DataFiles.h>
#include <sam_param.h>


void	SetHullOrigin(
    char        *DSname,
    HeaderInfo  *Header,
    ChannelInfo *Channel,
    HULL        *Hull
) {
    double  Vo[3];  // sphere origin
	int     i;
	int		v;
	char	Label[256];
	char	Name[256];
	FILE	*fp;

	// open head model file
	sprintf(Name, "%s.ds/default.hdm", DSname);
	if((fp = fopen(Name, "r")) == NULL)
		cleanup("'default.hdm' open failed");

	// read single sphere values
    do {    // synch up to line containing 'ORIGIN_X'
        if (fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'ORIGIN_X' label");
    } while (strncmp(Label, "ORIGIN_X:", 9) != 0);
    if (fscanf(fp, "%lf", &Vo[X_]) != 1)    // read sphere x-origin
        cleanup("can't read sphere x-origin");

    do {    // synch up to line containing 'ORIGIN_Y
        if (fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'ORIGIN_Y' label");
    } while (strncmp(Label, "ORIGIN_Y:", 9) != 0);
    if (fscanf(fp, "%lf", &Vo[Y_]) != 1)     // read sphere y-origin
        cleanup("can't read sphere y-origin");

    do {    // synch up to line containing 'ORIGIN_Z
        if (fscanf(fp, "%s", Label) != 1)
            cleanup("can't read 'ORIGIN_Z' label");
    } while (strncmp(Label, "ORIGIN_Z:", 9) != 0);
    if (fscanf(fp, "%lf", &Vo[Z_]) != 1)     // read sphere z-origin
        cleanup("can't read sphere z-origin");

    for (v=X_; v<=Z_; v++)       // set Hull origin to single sphere
        Hull->Vo[v] = Vo[v];

    // bonus! set sensor origins to single sphere
    for (i=0; i<Header->NumChannels; i++)
        if (Channel[i].ChannelType == TYPE_MEG
         || Channel[i].ChannelType == TYPE_REF_MAG
         || Channel[i].ChannelType == TYPE_REF_GRAD
         && Channel[i].Flag == TRUE)
            for (v=X_; v<=Z_; v++)
                Channel[i].Geom.MEGSensor.LocalSphere[v] = Vo[v] * 0.01;  // convert to metres

	// exit stage left...
	fclose(fp);
}
