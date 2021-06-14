// GetNIFTI() -- reads a NIFTI file

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <geoms.h>
#include <SAMfiles.h>
#include <nifti1.h>
#include <lfulib.h>

void	GetNIFTI(
	char			*ImageName,		// full pathname for NIFTI image file
	int				*nt,			// number of time or seed images
	float			***Data,		// image data
	nifti_1_header	*Header			// .nii header
)

{
	static float	**Image;		// Image
	int				d;				// 3d image index
	int				i;
	int				j;
	int				k;
	int				l;
	int				D;				// number of 3d images
	int				V;				// number of voxels per image
	int				offset;			// offset to start of voxels
	unsigned char	*NiiExt;		// optional extended header
	char			extension[4];	// 4-byte extension
	size_t			N;				// number of characters in NIFTI extended header`
	FILE			*fp;			// file pointer

	// open .nii file for read
	if((fp = fopen(ImageName, "r")) == NULL)
		cleanup("can't open nifti file for read");

	// read .nii header & extension
	if(fread((void *)Header, sizeof(nifti_1_header), 1, fp) != 1)
		cleanup("can't read nifti header");
	if(fread((void *)extension, 4, 1, fp) != 1)
		cleanup("can't read nifti extension");
	if(Header->datatype != DT_FLOAT)
		cleanup("only float datatype supported");
	offset = (int)Header->vox_offset;
	if(offset > 352) {
		N = (size_t)(offset - 352);
		if((NiiExt = (unsigned char *)malloc(N * sizeof(unsigned char))) == NULL)
			allocfailed("NiiExt[]");
		if(fread((void *)NiiExt, N, 1, fp) != 1)
			cleanup("can't read extended NIFTI header");
	}

	// get dimensions
	i = Header->dim[1];
	j = Header->dim[2];
	k = Header->dim[3];
	l = Header->dim[4];				// number of time or seed images
	D = (l == 0)? 1: l;
	V = i * j * k;					// number of voxels per image

	// allocate memory for image data
	if((Image = (float **)malloc((size_t)D * sizeof(float *))) == NULL)
		allocfailed("Image[]");
	for(d=0; d<D; d++)
		if((Image[d] = (float *)malloc((size_t)V * sizeof(float))) == NULL)
			allocfailed("Image[][]");

	// read image data
	for(d=0; d<D; d++)
		if(fread((void *)Image[d], sizeof(float), V, fp) != V)
			cleanup("can't read image data");
	fclose(fp);

	// output pointers
	*Data = Image;
	*nt = D;
}
