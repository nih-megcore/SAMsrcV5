// MoveMEG() -- copy all SQUID sensors to new frame
//
//  All SQUID sensors (reference & primary sensors) are copied from the original
//  reference frame into a new frame. The original sensor structure is kept intact.
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <geoms.h>
#include <samlib.h>

void MoveMEG(
    ChannelInfo     *Input,     // Input[S] -- input ChannelInfo structure
    ChannelInfo     *Output,    // Output[S] -- output ChannelInfo structure
    int             *Index,     // list of all SQUID sensor indices
    int             S,          // number of sensors
    XFORM           *where      // where the probe is to be moved
) {
    int             c;          // coil index
    int             i;          // channel index
    register int    s;          // sensor index
    register int    v;          // vector index

    // move all SQUID sensors
    for (s=0; s<S; s++) {
        i = Index[s];           // get channel index in array

        // copy input sensor information to output
        Output[i].Geom.MEGSensor.NumCoils = Input[i].Geom.MEGSensor.NumCoils;
        for (v=X_; v<=Z_; v++) {     // copy sensor origin & local sphere
            Output[i].Geom.MEGSensor.origin.p[v] = Input[i].Geom.MEGSensor.origin.p[v];
            Output[i].Geom.MEGSensor.origin.v[v] = Input[i].Geom.MEGSensor.origin.v[v];
            Output[i].Geom.MEGSensor.LocalSphere[v] = Input[i].Geom.MEGSensor.LocalSphere[v];
        }
        for (c=0; c<Input[i].Geom.MEGSensor.NumCoils; c++) {
            for (v=X_; v<=Z_; v++) {
                Output[i].Geom.MEGSensor.Coil[c].origin.p[v] = Input[i].Geom.MEGSensor.Coil[c].origin.p[v];
                Output[i].Geom.MEGSensor.Coil[c].origin.v[v] = Input[i].Geom.MEGSensor.Coil[c].origin.v[v];
            }
            Output[i].Geom.MEGSensor.Coil[c].Evaluated = FALSE;
        }

        // translate & rotate sensor to 'where'
        MoveSensor(&Output[i].Geom.MEGSensor, where);
    }
}
