// version.h -- version string
//
//      Author: SE Robinson
//                      MEG Core Facility
//                      NIMH
//

#ifndef H_VERSION
#define H_VERSION

#define SAM_REV 3
#define _VERSION_STR "5.0"

#if __SIZEOF_POINTER__ == 4
#define PRG_REV "Version " _VERSION_STR " (32-bit)"
#else
#define PRG_REV "Version " _VERSION_STR " (64-bit)"
#endif

#endif  // H_VERSION
