#include <stdlib.h>
#include <gsl/gsl_version.h>

int main(int a, char **v)
{
#if GSL_MINOR_VERSION == 6
    exit(0);
#else
    exit(1);
#endif
}
