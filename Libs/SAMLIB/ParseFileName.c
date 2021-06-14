// ParseFileName() -- parse a file name from its directory pathname
//	name from a full path name
//
//	Author:	Stephen E. Robinson
//			MEG Core Facility
//			NIMH
//

#include <stdio.h>
#include <string.h>

#define FALSE	0
#define TRUE	1

void	ParseFileName(
	char	*PathName,			// input full path name for dataset
	char	*FileName			// file name (last token)
)

{
	int		found;				// delimiter flag
	int		i;					// character index
	int		j;					// string index
	int		n;					// length of PathName
	int		last;				// index of last delimiter

	// find the last '/' delimiter in the PathName string
	n = strlen(PathName);
	for(i=last=0, found=FALSE; i<n; i++)	// find last occurance of '/'
		if(PathName[i] == '/') {
			last = i;
			found = TRUE;
		}

	// get FileName
	if(found == FALSE) {					// we already have the file name
		strcpy(FileName, PathName);
	} else {								// we have the delimiter -- copy the following characters
		i = n - last - 1;
		for(j=0; j<i; j++, last++)
			FileName[j] = PathName[last+1];
		FileName[j] = '\0';
	}
}
