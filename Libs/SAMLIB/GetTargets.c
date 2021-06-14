// GetTargets() -- reads list of target coordinates (cm) into array
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <geoms.h>
#include <voxel.h>
#include "samutil.h"

void	GetTargets(
	char		*PathName,			// full path name for target file (input)
	VOXELINFO	**TargetList,		// address of target array (output)
	int			*N					// number of targets (output)
	)

{
	static VOXELINFO	*Target;
	double				Vc[3];		// cartesian coordinates
	double              val;
	int					n;
	int					v;
	char				Line[128];	// space for one line of coordinates
	FILE				*fp;

	// open target file for read
	if((fp = fopen(PathName, "r")) == NULL)
		cleanup("target file open failed");

	// read number of uncommented lines in file
	n = 0;
	while(fgets(Line, 128, fp) != NULL) {
		if(Line[0] == '#')
			continue;
 		n++;
 	}
	if(n == 0)
		cleanup("empty target file");
	*N = n;

	// allocate space for target coordinates
	if((Target = (VOXELINFO *)malloc((size_t)n * sizeof(VOXELINFO))) == NULL)
		allocfailed("Target[]");
	*TargetList = Target;

	// reset to start of target file
	rewind(fp);
	
	// read target list into array
	n = 0;
	while(fgets(Line, 128, fp) != NULL) {
		if(sscanf(Line, "%lf%lf%lf%le", &Vc[X_], &Vc[Y_], &Vc[Z_], &val) != 4)
			cleanup("can't read target coordinates");
		for(v=X_; v<=Z_; v++) {
			Target[n].p[v] = 0.01 * Vc[v];
			Target[n].v[v] = 0.;
		}
		Target[n].Solve = TRUE;
		Target[n].ROI = TRUE;
		n++;
	}
	fclose(fp);
}
