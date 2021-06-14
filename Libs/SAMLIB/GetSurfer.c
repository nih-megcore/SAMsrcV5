#include <stdio.h>
#include <stdlib.h>
#include <geoms.h>
#include <voxel.h>
#include "samutil.h"

void	GetSurfer(
	char		*PathName,			// full path name for surfer file (input)
	VOXELINFO	**VertexList,		// address of surfer vertex array (output)
	int			*N					// number of vertices (output)
)
{
	static VOXELINFO	*Vertex;
	int					n;
	int					code;		// unused
	char				Line[256];	// space for one line of coordinates
	FILE				*fp;

	// read number of vertices from *.norm
	if((fp = fopen(PathName, "r")) == NULL)
		Cleanup("can't open surface mesh with norms'");
	n = 0;
	while(fgets(Line, 256, fp) != NULL)		// count lines
 		n++;
	if(n == 0)
		cleanup("empty vertex file");
	*N = n;

	// allocate space for vertices
	if((Vertex = (VOXELINFO *)malloc((size_t)*N * sizeof(VOXELINFO))) == NULL)
		allocfailed("VertexList[]");
	*VertexList = Vertex;

	// reset to start of vertex file
	rewind(fp);

	// read target list into array
	for(n=0; n<(*N); n++) {
		if(fscanf(fp, "%lf%lf%lf%lf%lf%lf%d", &Vertex[n].p[X_], &Vertex[n].p[Y_], &Vertex[n].p[Z_], &Vertex[n].v[X_], &Vertex[n].v[Y_], &Vertex[n].v[Z_], &code) != 7)
			cleanup("error reading vertex");
		Vertex[n].Solve = FALSE;
		Vertex[n].ROI = TRUE;
	}
	fclose(fp);
}
