// GetNIFTIWts() -- reads a NIFTI SAM coefficient file file
//	This version returns a pointer to a GSL matrix that is
//	statically allocated, containing the beamformer coefficients.
//
// Author:	Stephen E. Robinon
//			MEG Core Facility
//			NIMH
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>
#include <SAMfiles.h>
#include <nifti1.h>
#include <lfulib.h>

gsl_matrix	*GetNIFTIWts(
	char			*WgtFilePath,	// full pathname for NIFTI coefficient file (input)
	nifti_1_header	*Header,		// .nii header
	char			Ext[4],			// 4-character extension
	unsigned char	**ExtHeader		// optional extended header
)

{
	static gsl_matrix	*Wgt;		// Wgt - VxM beamformer weights (static to persist for multiple calls)
	unsigned char	*NiiExt;		// NiiExt[N] -- optional extended header
	float			voxel;			// coefficient
	int				i;
	int				j;
	int				k;
	int				m;				// channel index
	int				v;				// voxel index
	int				M;				// number of channels
	int				V;				// number of coefficients per channel
	int				offset;			// offset to start of voxels
	size_t			N;				// number of characters in NIFTI extended header`
	FILE			*fp;			// file pointer

	// open .nii file for read
	if((fp = fopen(WgtFilePath, "r")) == NULL)
		cleanup("can't open .nii file for read");

	// read .nii header & extension
	if(fread((void *)Header, sizeof(nifti_1_header), 1, fp) != 1)
		cleanup("can't read .nii header");
	if(fread((void *)Ext, 4, 1, fp) != 1)
		cleanup("can't read .nii extension");
	if(Header->datatype != DT_FLOAT)
		cleanup("only float datatype supported");
	offset = (int)Header->vox_offset;
	if(offset > 352) {
		N = (size_t)(offset - 352);
		if((NiiExt = (unsigned char *)malloc(N * sizeof(unsigned char))) == NULL)
			allocfailed("NiiExt[]");
		if(fread((void *)NiiExt, N, 1, fp) != 1)
			cleanup("can't read extended NIFTI header");
		*ExtHeader = NiiExt;
	} else {
		*ExtHeader = NULL;
	}

	// get dimensions
	i = Header->dim[1];
	j = Header->dim[2];
	k = Header->dim[3];
	M = Header->dim[4];				// number of channels
	V = i * j * k;					// number of coefficients per channel

	// read coefficients
	Wgt = gsl_matrix_alloc(V, M);
	for(m=0; m<M; m++)
		for(v=0; v<V; v++) {
			if(fread((void *)&voxel, sizeof(float), 1, fp) != 1)
				cleanup("can't read coefficient from .nii file");
			gsl_matrix_set(Wgt, v, m, (double)voxel);
		}

	// done
	fclose(fp);
	return(Wgt);
}
