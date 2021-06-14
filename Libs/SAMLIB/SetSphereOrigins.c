// SetSphereOrigins() -- sets sphere origins from *.hdm file
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>      // get axes definitions
#include <samlib.h>
#include <lfulib.h>
#include <DataFiles.h>


void SetSphereOrigins(
    char        *name,  // default.hdm path
    HeaderInfo  *Header,
    ChannelInfo *Channel,
    int         Multi
) {
    double  Vo[3];  // sphere origin
    double  r;
    int     i;
    int     n;
    int     v;
    int     count;
    char    Name[256];
    char    Label[256];
    FILE    *fp;

    // open head model file

    fp = fileopen(name, "r");

    // read single (or mean) local sphere values
    if (!Multi) {   // single sphere
        do {                                                    // synch up to line containing 'ORIGIN_X'
            if (fscanf(fp, "%s", Label) != 1)
                cleanup("can't read 'ORIGIN_X' label");
        } while (strncmp(Label, "ORIGIN_X:", 9) != 0);
        if (fscanf(fp, "%lf", &Vo[X_]) != 1)     // read sphere x-origin
            cleanup("can't read sphere x-origin");
        do {
            if (fscanf(fp, "%s", Label) != 1)
                cleanup("can't read 'ORIGIN_Y' label");
        } while (strncmp(Label, "ORIGIN_Y:", 9) != 0);
        if (fscanf(fp, "%lf", &Vo[Y_]) != 1)
            cleanup("can't read sphere y-origin");
        do {
            if (fscanf(fp, "%s", Label) != 1)
                cleanup("can't read 'ORIGIN_Z' label");
        } while (strncmp(Label, "ORIGIN_Z:", 9) != 0);
        if (fscanf(fp, "%lf", &Vo[Z_]) != 1)
            cleanup("can't read sphere z-origin");
        for (i=0; i<Header->NumChannels; i++)
            if (Channel[i].ChannelType == TYPE_MEG
             || Channel[i].ChannelType == TYPE_REF_MAG
             || Channel[i].ChannelType == TYPE_REF_GRAD
             && Channel[i].Flag == TRUE)
                for (v=X_; v<=Z_; v++)
                    Channel[i].Geom.MEGSensor.LocalSphere[v] = Vo[v] * 0.01;  // convert to metres

    } else {    // multisphere

        // synch up to start of multisphere list
        do {
            if (fscanf(fp, "%s", Label) != 1)
                cleanup("can't read 'Radius' label");
        } while (strncmp(Label, "Radius", 6) != 0);     // this is the last word in header line

        // read values for each channel in multisphere list
        do {
            count = fscanf(fp, "%s%lf%lf%lf%lf", Name, &Vo[X_], &Vo[Y_], &Vo[Z_], &r);
            if (count == 5) {
                n = strlen(Name);
                for (i=0; i<Header->NumChannels; i++) {
                    if (!strncmp(Name, Channel[i].ChannelName, n-1)) {   // ignore trailing colon
                        for (v=X_; v<=Z_; v++)
                            Channel[i].Geom.MEGSensor.LocalSphere[v] = 0.01 * Vo[v];
                        break;
                    }
                }
            }
        } while(count == 5);
    }

    // exit stage left...
    fclose(fp);
}
