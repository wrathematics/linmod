#include <stdint.h>
#include <string.h>


/* double */

void r_set_na_real_(double *val)
{
  int64_t x = 0x7FF00000000007A2LL;
  memcpy((void *) val, (void *) &x, 8);
}

void r_set_nan_real_(double *val)
{
  int64_t x = 0x7FF0000000000000LL;
  memcpy((void *) val, (void *) &x, 8);
}



/* int */

void r_set_na_int_(int *val)
{
  *val = INT32_MIN;
}

