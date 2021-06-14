// PutNIFTIWts() -- writes the header & SAM coefficients as a NIFTI file
//  This version passes the beamformer coefficients as a GSL matrix structure
//
// Author:  Stephen E. Robinson
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>
#include <SAMfiles.h>
#include <gsl/gsl_matrix.h>
#include <nifti1.h>
#include "samlib.h"
#include "samutil.h"
#include "sam_param.h"

#define ETA 1.0e-15


void    PutNIFTIWts(
    char        *WgtFilePath,       // full path for writing NIFTI file
    PARMINFO    *Params,            // analysis parameters
    gsl_matrix  *Wgt,               // Wgt - VxM SAM coefficients
    float       Noise               // sensor mean noise floor
)

{
    nifti_1_header  Header;         // output header
    float       Vc[3];
    float       XStart;
    float       XEnd;
    float       YStart;
    float       YEnd;
    float       ZStart;
    float       ZEnd;
    float       Step;
    float       voxel;
    int         i;
    int         j;
    int         k;
    int         M;                  // number of primary sensors in Wgt
    int         V;
    int         m;                  // channel index
    int         v;                  // voxel index
    char        extension[4];       // 4-byte extension
    FILE        *fp;                // file pointer

    // get dimensions from Wgt
    M = Wgt->size2;
    
    // explicitly zero header & extension
    memset((void *)&Header, 0, 348);
    memset((void *)extension, 0, 4);

    // conversion to mm -- for AFNI's sake
    XStart = 1000. * (float)Params->SAMStart[X_];
    XEnd = 1000. * (float)Params->SAMEnd[X_];
    YStart = 1000. * (float)Params->SAMStart[Y_];
    YEnd = 1000. * (float)Params->SAMEnd[Y_];
    ZStart = 1000. * (float)Params->SAMStart[Z_];
    ZEnd = 1000. * (float)Params->SAMEnd[Z_];
    Step = 1000. * (float)Params->SAMStep;

    // count number of voxels in ROI (in mm)
    for(Vc[X_]=XStart, k=0; Vc[X_]<=(XEnd+ETA); Vc[X_]+=Step, k++);
    for(Vc[Y_]=YStart, j=0; Vc[Y_]<=(YEnd+ETA); Vc[Y_]+=Step, j++);
    for(Vc[Z_]=ZStart, i=0; Vc[Z_]<=(ZEnd+ETA); Vc[Z_]+=Step, i++);
    V = i * j * k;
    if(V != Wgt->size1)
        cleanup("number of voxels in Wgt does not agree with computed from dimensions");

    // fill in NIFTI header
    Header.sizeof_hdr = 348;                        // fixed nift1 header size
    //  Header.data_type was cleared, above
    //  Header.db_name was cleared, above
    Header.extents = 0;
    Header.session_error = 0;
    Header.regular = (char)0;
    Header.dim_info = NIFTI_SLICE_SEQ_INC;
    Header.intent_p1 = NIFTI_INTENT_NONE;
    Header.intent_p2 = NIFTI_INTENT_NONE;
    Header.intent_p3 = NIFTI_INTENT_NONE;
    Header.intent_code = NIFTI_INTENT_NONE;
    Header.dim[0] = 4;                              // this is 3d+coefficients
    Header.dim[1] = i;                              // slowest changing index
    Header.dim[2] = j;                              // next
    Header.dim[3] = k;                              // most rapidily changing index
    Header.dim[4] = M;                              // time or seed number (must be 1 if no time/seed dimension)
    Header.datatype = DT_FLOAT;
    Header.bitpix = 32;                             // redundant with float datatype
    Header.slice_start = 0;                         // first slice index
    Header.pixdim[0] = 1.;                          // qfac
    Header.pixdim[1] = Step;                        // i-axis SAM step size in mm
    Header.pixdim[2] = Step;                        // j-axis "
    Header.pixdim[3] = Step;                        // k-axis "
    Header.pixdim[4] = 1.;                          // weight step (index step)
    Header.vox_offset = 352.;                       // offset into .nii file (why did they make this float?)
    Header.scl_slope = 0.;                          // data scaling: slope
    Header.scl_inter = 0.;                          // data scaling: offset
    Header.slice_end = M - 1;                       // last weight index
    Header.slice_code = NIFTI_SLICE_SEQ_INC;        // voxels are sequential in xyzt
    Header.xyzt_units = (NIFTI_UNITS_MM | NIFTI_UNITS_UNKNOWN);
    Header.cal_max = 1000.;                         // max display intensity
    Header.cal_min = -1000.;                        // min display intensity
    Header.slice_duration = 1.;                     // time for 1 slice
    Header.toffset = 0.;                            // time axis shift
    Header.noise = Noise * 1.0e+28;                 // re-purposed glmax!
    Header.glmin = 0;                               // unused
    // Header.descrip was set to zero, above
    // Header.aux_file was set to zero, above
    Header.qform_code = NIFTI_XFORM_SCANNER_ANAT;   // use method 2 for transform
    Header.sform_code = NIFTI_XFORM_UNKNOWN;        // don't use method 3
    // Header.descrip was set to zero, above
    // Header.aux_file was set to zero, above
    Header.qform_code = NIFTI_XFORM_UNKNOWN;        // don't use method 2
    Header.sform_code = NIFTI_XFORM_SCANNER_ANAT;   // use method 3 for transform
    // Header.quatern_b, Header.quatern_c, Header.quatern_d were set to zero, above
    // Header.qoffset_x, Header.qoffset_y, Header.qoffset_z were set ti zero, above

    // use method 3 -- affine transform
    //  +x = Right  +y = Anterior  +z = Superior
    Header.srow_x[0] = 0.;      Header.srow_x[1] = -Step;   Header.srow_x[2] = 0.;      Header.srow_x[3] = -YStart;
    Header.srow_y[0] = 0.;      Header.srow_y[1] = 0.;      Header.srow_y[2] = Step;    Header.srow_y[3] = XStart;
    Header.srow_z[0] = Step;    Header.srow_z[1] = 0.;      Header.srow_z[2] = 0.;      Header.srow_z[3] = ZStart;

    bcopy((void *)"SAM Coefficients", (void *)Header.intent_name, 16);
    bcopy("n+1\0", Header.magic, 4);

    // open NIFTI file for write
    if((fp = fopen(WgtFilePath, "w")) == NULL)
        cleanup("can't open .nii file for write");

    // write header & extension
    if(fwrite((void *)&Header, sizeof(nifti_1_header), 1, fp) != 1)
        cleanup("can't write nifti header to .nii file");
    if(fwrite((void *)extension, 4, 1, fp) != 1)
        cleanup("can't write 'extension' to .nii file");

    // write the SAM coefficients -- note the order for each channel as a '3d image' this is compatible with AFNI
    for(m=0; m<M; m++) {
        for(v=0; v<V; v++) {
            voxel = (float)gsl_matrix_get(Wgt, v, m);
            if(fwrite((void *)&voxel, sizeof(float), 1, fp) != 1)
                cleanup("can't write SAM coefficient to .nii file");
        }
        fflush(fp);
    }

    // done
    fclose(fp);
}
