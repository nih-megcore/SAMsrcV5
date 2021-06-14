// MoveSensor() -- move a sensor to specified position & rotation
//  relative to starting frame
//
//  Author: Stephen E. Robinson
//          Neuromagnetism Laboratory
//          Dept.of Neurology
//          Henry Ford Hospital
//          Detroit, MI
//          Copyright (c) 2007
//

#include <stdio.h>
#include <geoms.h>
#include <samlib.h>

void MoveSensor(
    SENSOR  *Sensor,        // Sensor structure
    XFORM   *where          // new origin vectors
) {
    register int c;         // coil index

    // move sensor origin vectors
    MoveField(&(Sensor->origin), where);

    // apply rotation matrix to each coil origin FIELD
    for (c=0; c<Sensor->NumCoils; c++) {
        MoveField(&(Sensor->Coil[c].origin), where);
        Sensor->Coil[c].Evaluated = FALSE;
    }
}
