// PutNIFTIImage -- writes a 3d NIFTI image file

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>
#include <SAMfiles.h>
#include <nifti1.h>
#include "samutil.h"
#include "sam_param.h"

#define ETA	1.0e-12


void	PutNIFTIImage(
	float			*Image,			// Image[V]
	char			*fpath,			// full pathname
	PARMINFO		*Params,		// parameters
	char			*intent			// intent of data (i.e., "mean" or "variance")
)

{
	nifti_1_header	NiiHdr;			// NIFTI header
//	unsigned char	*ExtHdr;		// NiiHdr[N] -- extended NIFTI header
	char			extension[4];	// extension
	double			Vc[3];
	float			XStart;
	float			XEnd;
	float			YStart;
	float			YEnd;
	float			ZStart;
	float			ZEnd;
	float			Step;
	int				i;
	int				j;
	int				k;
//	int				N;
	int				V;
	size_t			len;
	FILE			*fp;			// file pointer

	// explicitly zero header
	memset(&NiiHdr, 0, 348);							// explicitly zero header
	memset(extension, 0, 4);							// for now, ignore the fact that we may have an extended header

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

	// fill in NIFTI header
	NiiHdr.sizeof_hdr = 348;								// fixed nift1 header size
	//	NiiHdr.data_type was cleared, above
	//	NiiHdr.db_name was cleared, above
	NiiHdr.extents = 0;
	NiiHdr.session_error = 0;
	NiiHdr.regular = (char)0;
	NiiHdr.dim_info = NIFTI_SLICE_SEQ_INC;
	NiiHdr.intent_p1 = NIFTI_INTENT_NONE;
	NiiHdr.intent_p2 = NIFTI_INTENT_NONE;
	NiiHdr.intent_p3 = NIFTI_INTENT_NONE;
	NiiHdr.intent_code = NIFTI_INTENT_NONE;
	NiiHdr.dim[0] = 3;										// 3D volume
	NiiHdr.dim[1] = i;										// slowest changing index:			(x-axis)
	NiiHdr.dim[2] = j;										// next:							(y-axis)
	NiiHdr.dim[3] = k;										// most rapidily changing index:	(z-axis)
	NiiHdr.dim[4] = 1;										// time or seed number (must be 1 if no time/seed dimension)
	NiiHdr.datatype = DT_FLOAT;								// floats should be sufficient
	NiiHdr.bitpix = 32;										// redundant with float datatype
	NiiHdr.slice_start = 0;									// first slice index
	NiiHdr.pixdim[0] = 1.;									// qfac
	NiiHdr.pixdim[1] = Step;								// i-axis SAM step size
	NiiHdr.pixdim[2] = Step;								// j-axis "
	NiiHdr.pixdim[3] = Step;								// k-axis "
	NiiHdr.pixdim[4] = 1.;									// time step (else 1.) for seed images
	NiiHdr.vox_offset = 352.;								// offset into .nii file (why did they make this float?)
	NiiHdr.scl_slope = 0.;									// data scaling: slope
	NiiHdr.scl_inter = 0.;									// data scaling: offset
	NiiHdr.slice_end = V - 1; 								// last slice index
	NiiHdr.slice_code = NIFTI_SLICE_SEQ_INC;				// voxels are sequential in xyzt
	NiiHdr.xyzt_units = NIFTI_UNITS_MM;						// units of pixdim[1..4]
	NiiHdr.cal_max = 1000.;									// max display intensity
	NiiHdr.cal_min = -1000.;									// min display intensity
	NiiHdr.slice_duration = 1.;								// time for 1 slice
	NiiHdr.toffset = 0.;									// time axis shift
//	NiiHdr.glmax = 0;										// unused
	NiiHdr.glmin = 0;										// unused
	// NiiHdr.descrip was set to zero, above
	// NiiHdr.aux_file was set to zero, above
	NiiHdr.qform_code = NIFTI_XFORM_SCANNER_ANAT;			// use method 2 for transform
	NiiHdr.sform_code = NIFTI_XFORM_UNKNOWN;				// don't use method 3
	// NiiHdr.descrip was set to zero, above
	// NiiHdr.aux_file was set to zero, above
	NiiHdr.qform_code = NIFTI_XFORM_UNKNOWN;				// don't use method 2
	NiiHdr.sform_code = NIFTI_XFORM_SCANNER_ANAT;			// use method 3 for transform
	// NiiHdr.quatern_b, Header.quatern_c, Header.quatern_d were set to zero, above
	// NiiHdr.qoffset_x, Header.qoffset_y, Header.qoffset_z were set to zero, above

	// use method 3 -- affine transform
	NiiHdr.srow_x[0] = 0.;		NiiHdr.srow_x[1] = -Step;	NiiHdr.srow_x[2] = 0.;		NiiHdr.srow_x[3] = -YStart;
	NiiHdr.srow_y[0] = 0.;		NiiHdr.srow_y[1] = 0.;		NiiHdr.srow_y[2] = Step;	NiiHdr.srow_y[3] = XStart;
	NiiHdr.srow_z[0] = Step;	NiiHdr.srow_z[1] = 0.;		NiiHdr.srow_z[2] = 0.;		NiiHdr.srow_z[3] = ZStart;

	// set intent
	len = strlen((char *)intent);
	bcopy((void *)intent, (void *)NiiHdr.intent_name, len);
	bcopy("n+1\0", NiiHdr.magic, 4);

	// open output file
	if((fp = fopen(fpath, "w")) == NULL)
		cleanup("can't open .nii file for write");

	// write header & extension
	if(fwrite((void *)&NiiHdr, sizeof(nifti_1_header), 1, fp) != 1)
		cleanup("can't write nifti header to .nii file");
	if(fwrite((void *)extension, 4, 1, fp) != 1)
		cleanup("can't write 'extension' to .nii file");
#if 0
	if((int)NiiHdr.vox_offset != 352 && ExtHdr != NULL) {
		N = (int)NiiHdr.vox_offset - 352;
		if(fwrite((void *)ExtHdr, N, 1, fp) != 1)
			cleanup("can't write extended header to .nii file");
	}
#endif
	if(fwrite((void *)Image, sizeof(float), V, fp) != V)
		cleanup("can't write data image file");
	fclose(fp);
}
