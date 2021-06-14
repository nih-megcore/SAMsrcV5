// Mac OS X definitions for missing byteswap functions found in Linux

#include <stdint.h>

#ifndef bswap_16
# define bswap_16(x) ( \
	((uint16_t)(x) & 0x00ffU) << 8 | \
	((uint16_t)(x) & 0xff00U) >> 8)
#endif

#ifndef bswap_32
# define bswap_32(x) ( \
	((uint32_t)(x) & 0x000000ffU) << 24 | \
	((uint32_t)(x) & 0x0000ff00U) << 8 | \
	((uint32_t)(x) & 0x00ff0000U) >> 8 | \
	((uint32_t)(x) & 0xff000000U) >> 24)
#endif

#ifndef bswap_64
# define bswap_64(x) ( \
	((uint64_t)(x) & 0x00000000000000ffULL) << 56 | \
	((uint64_t)(x) & 0x000000000000ff00ULL) << 40 | \
	((uint64_t)(x) & 0x0000000000ff0000ULL) << 24 | \
	((uint64_t)(x) & 0x00000000ff000000ULL) << 8 | \
	((uint64_t)(x) & 0x000000ff00000000ULL) >> 8 | \
	((uint64_t)(x) & 0x0000ff0000000000ULL) >> 24 | \
	((uint64_t)(x) & 0x00ff000000000000ULL) >> 40 | \
	((uint64_t)(x) & 0xff00000000000000ULL) >> 56)
#endif
