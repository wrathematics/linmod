#include "linmod.h"
#include <math.h>

#define MAX(m,n) m<n?n:m


SEXP make_lmfit_default_rownames(const int n)
{
  R_INIT;
  int i;
  int buflen;
  SEXP ret;
  
  buflen = (int) (ceil(log10((double)n)) + 1.);
  char *buf = malloc(buflen * sizeof(buf));
  buf[0] = 'x';
  
  newRlist(ret, n);
  
  for (i=0; i<n; i++)
  {
    sprintf(buf+1, "%d", i+1);
    buflen = (int) (ceil(log10((double)i+2)) + 1.);
    buflen = MAX(buflen, 2);
    SET_VECTOR_ELT(ret, i, mkCharLen(buf, buflen));
  }
  
  free(buf);
  
  R_END;
  return ret;
}

