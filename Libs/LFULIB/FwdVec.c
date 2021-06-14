// FwdVec() -- compute forward solution vector for a gradiometer
//  (with reference compensation) due to current dipole in homogeneously
//  conducting sphere. This subroutine supports multiple local spheres.
//
//  This version features programmable integration order from 1st to 7th
//  for SQUID sensor flux transformers. If order=0, use OPM forward solution.
//
//  Author: Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <geoms.h>
#include <samlib.h>

void FwdVec(
    FIELD       *Q6,        // location & orientation of current dipole (input)
    double      *Bp,        // Bp[M] -- forward solution vector (output to externally allocated vector)
    HeaderInfo  *Header,    // header information, including primary & reference sensor indices (input)
    ChannelInfo *Channel,   // Channel[M] -- channel structure with primary and reference sensors (input)
    int         Order       // integration order: 1, 2, 3, 4, 5, & 7 are implemented (input)
) {
    double      *Br;        // Br[R] -- reference sensor field array
    int         i;          // index
    int         j;          // another index
    int         k;          // yet another index
    int         g;          // still going...
    int         M;          // number of primary sensors
    int         R;          // number of reference sensors
    int         firstRef;   // index of first reference

    // allocate memory for reference sensor field Br
    M = Header->NumPri;
    R = Header->NumRef;
    Br = new_arrayE(double, R, "Br[]");
    firstRef = Header->RsIndex[0];

    // compute reference channel forward solution for Target source
    for (i=0; i<R; i++) {
        j = Header->RsIndex[i];
        if (Order == 0)
            Br[i] = FwdOPM(Q6, &Channel[j].Geom.OPMSensor);
        else
            Br[i] = FwdSensor(Q6, &Channel[j].Geom.MEGSensor, Order);
    }

    // compute primary channel (balanced) forward solution for Target source
    for (i=0; i<M; i++) {
        j = Header->PsIndex[i];
        if (Order == 0)
            Bp[i] = FwdOPM(Q6, &Channel[j].Geom.OPMSensor);
        else
            Bp[i] = FwdSensor(Q6, &Channel[j].Geom.MEGSensor, Order);
        for (g=0; g<Channel[j].NumBal; g++) {
            k = Channel[j].Balance[g].RefChan - firstRef;
            Bp[i] -= Channel[j].Balance[g].Weight * Br[k];
        }
    }

    free(Br);
}
