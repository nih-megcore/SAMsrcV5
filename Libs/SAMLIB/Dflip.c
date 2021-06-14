// Dflip() -- flip endian order for a double
//
//	Author:	Stephen E. Robinson
//		MEG Core Facility
//		NIMH
//

#include <stdint.h>
#include <byteswap.h>


double	Dflip(
	union {
		double dval;
		unsigned char *cval;
	} value
)

{
	unsigned char buffer[8];

	buffer[0] = value.cval[7];
	buffer[1] = value.cval[6];
	buffer[2] = value.cval[5];
	buffer[3] = value.cval[4];
	buffer[4] = value.cval[3];
	buffer[5] = value.cval[2];
	buffer[6] = value.cval[1];
	buffer[7] = value.cval[0];

	return *((double *)&buffer);
}
