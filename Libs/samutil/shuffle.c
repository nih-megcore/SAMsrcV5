#include <stdlib.h>

/* Create a random permutation by making n - 1 random exchanges (Knuth). */

void shuffle(int *array, int n)
{
    int i, j, t;

    for (i = n - 1; i > 0; --i) {

        /* Select a random array element, from 0 to i. */

        j = random() % (i + 1);

        /* Exchange. */

        t = array[i];
        array[i] = array[j];
        array[j] = t;
    }
}
