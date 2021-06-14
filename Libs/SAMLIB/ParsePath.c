// ParsePath() -- replacement for 'samUtils() to handle both absolute
//  & relative paths
//
// Author:      Stephen E. Robinson (improved by Tom Holroyd)
//                      Neuromagnetism Laboratory
//                      Dept. of Neurology
//                      Henry Ford Hospital
//                      Detroit, MI
//                      Copyright (c) 2009
//

#include <stdlib.h>
#include <string.h>

int ParsePath(
    char *Path,                 // input path, replaced by setname
    char *Dir                   // output dirname
)
{

    int i;                      // a convenient index
    int j;                      // and another
    int n;                      // number of characters in string
    int last;                   // last '/' in path name
    int found;                  // flag if '/' is found in 'Path'
    char *s;
    char *root;
    char pathname[256];
    char setname[256];

    // clear string buffers
    memset((void *)pathname, 0, 256);
    memset((void *)setname, 0, 256);

    // make sure something was passed
    if (Path == NULL || Path[0] == '\0' || Dir == NULL) {
	return -1;
    }
    // first, determine if this is absolute or relative path
    root = NULL;
    if (Path[0] != '/') {       // an absolute path begins with '/'
	if ((s = getenv("PWD")) == NULL)
	    return -1;
	root = (char *)malloc((size_t) strlen(s) + 2);
	strcpy(root, s);
	strcat(root, "/");      // append a slash to root path
    }
    // next, if req'd, split the dataset name from the rest of the path string
    n = strlen(Path);
    if (Path[n - 1] == '/') {   // if there is a terminal '/', remove it
	Path[n - 1] = 0;
	n--;
    }
    for (i = last = found = 0; i < n; i++)      // find last occurance of '/'
	if (Path[i] == '/') {
	    last = i;
	    found = 1;
	}
    // if req'd, separate path & dataset name at index of last '/'
    if (found == 1) {
	if (root != NULL) {     // relative pathname
	    strcpy(pathname, root);     // start with root
	}
	strncat(pathname, Path, last);  // add dir part

	for (i = 0, j = last + 1; i < n; i++, j++)      // and the set name goes here
	    setname[i] = Path[j];
    } else {                    // we are in the directory where the dataset resides
	strcpy(setname, Path);
	strcpy(pathname, root);
	n = strlen(pathname);
	pathname[n - 1] = 0;    // remove the offending terminal slash!
    }
    if (root)
	free(root);

    // if present, split '.ds' from setname, & copy to 'Path' (output)
    if ((n = strlen(setname)) == 0)     // can't have an empty setname, can we?
	return -1;
    if (!strncmp(&(setname[n - 3]), ".ds", 3)) {
	for (i = n - 3; i < n; i++)
	    setname[i] = 0;
    }
    strcpy(Path, setname);

    // lastly, we return the pathname, without setname
    strcpy(Dir, pathname);

    return 0;
}
