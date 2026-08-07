#include <stdint.h>
#include "sam_platform.h"

#ifdef _WIN32
static uint64_t sam_random_state = UINT64_C(0x4d595df4d0f33173);

static uint64_t sam_next_random(void)
{
    uint64_t value = sam_random_state;
    value ^= value >> 12;
    value ^= value << 25;
    value ^= value >> 27;
    sam_random_state = value;
    return value * UINT64_C(2685821657736338717);
}

void sam_srandom(unsigned int seed)
{
    sam_random_state = seed ? seed : UINT64_C(0x4d595df4d0f33173);
}

long sam_random(void)
{
    return (long)(sam_next_random() & UINT64_C(0x7fffffff));
}

void sam_srand48(long seed)
{
    sam_srandom((unsigned int)seed);
}

double sam_drand48(void)
{
    return (sam_next_random() >> 11) * (1.0 / 9007199254740992.0);
}
#endif
