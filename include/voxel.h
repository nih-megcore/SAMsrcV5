// voxel.h --definition of a voxel for dipole position & moment vector
//  Flags are used to control 1) linear solution for moment vector
//  & 2) whether voxel is within ROI bounds.
//
//  Author: SE Robinson
//          MEG Core Facility
//          NIMH
//

#ifndef H_VOXEL
#define H_VOXEL

typedef struct {
    double      p[3];       // vertex coordinate
    double      v[3];       // vertex normal
    int         Solve;      // use linear solution -- TRUE | FALSE
    int         ROI;        // voxel inside ROI -- TRUE | FALSE
} VOXELINFO;

#endif  // H_VOXEL
