// GetNIFTIMask() -- reads a NIFTI Mask file

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>
#include <SAMfiles.h>
#include <nifti1.h>
#include <samlib.h>

void	GetNIFTIMask(
	char			*MaskPath,		// full pathname for NIFTI mask file
	unsigned char	**Mask,			// mask data
	int				*nv				// number of voxels in mask
)

{
	nifti_1_header	Header;			// .nii header
	static unsigned char *Data;		// Image
	unsigned char	*NiiExt;		// NiiExt[N] -- optional extended header
	int				i;
	int				j;
	int				k;
	int				V;				// number of voxels per image
	int				offset;			// offset to start of voxels
	size_t			N;				// number of characters in NIFTI extended header`
	char			extension[4];	// 4-byte extension
	FILE			*fp;			// file pointer

	// open .nii file for read
	if((fp = fopen(MaskPath, "r")) == NULL)
		cleanup("can't open nifti file for read");

	// read .nii header & extension
	if(fread((void *)&Header, sizeof(nifti_1_header), 1, fp) != 1)
		cleanup("can't read nifti header");
	if(fread((void *)extension, 4, 1, fp) != 1)
		cleanup("can't read nifti extension");
	if(Header.datatype != DT_UNSIGNED_CHAR)
		cleanup("only unsigned char datatype supported");

	// if req'd, read extended header
	offset = (int)Header.vox_offset;
	if(offset > 352) {
		N = (size_t)(offset - 352);
		if((NiiExt = (unsigned char *)malloc(N * sizeof(unsigned char))) == NULL)
			allocfailed("NiiExt[]");
		if(fread((void *)NiiExt, N, 1, fp) != 1)
			cleanup("can't read extended NIFTI header");
	}

	// get dimensions
	i = Header.dim[1];
	j = Header.dim[2];
	k = Header.dim[3];
	V = i * j * k;					// number of voxels per image

	// allocate memory for mask data
	if((Data = (unsigned char *)malloc((size_t)V * sizeof(char))) == NULL)
		allocfailed("Data[]");

	// read image data
	if(fread((void *)Data, sizeof(char), V, fp) != V)
		cleanup("can't read image data");
	fclose(fp);

	// output pointers
	*Mask = Data;
	*nv = V;
}
