#include <stdio.h>
#include "hilbert.h"
 
long int evalhilbert_(bitmask_t const disc[], int *nbits, int *ndims)
{
    return hilbert_c2i(*ndims ,*nbits, disc);
}
