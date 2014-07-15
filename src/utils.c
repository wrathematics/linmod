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



SEXP make_lmfit_default_effectnames(const int m, const int n, const int *pvt)
{
  R_INIT;
  int i, j;
  int buflen;
  int maxpvt = n;
  SEXP ret;
  
  for (i=1; i<n; i++)
  {
    if (pvt[i] < pvt[i-1])
    {
      maxpvt = i;
      break;
    }
  }
  
  buflen = (int) (ceil(log10((double)n)) + 1.);
  char *buf = malloc(buflen * sizeof(buf));
  buf[0] = 'x';
  
  newRlist(ret, m);
  
  for (i=0; i<maxpvt; i++)
  {
    j = pvt[i];
    sprintf(buf+1, "%d", j);
    buflen = (int) (ceil(log10((double)j+1)) + 1.);
    buflen = MAX(buflen, 2);
    SET_VECTOR_ELT(ret, i, mkCharLen(buf, buflen));
  }
  
  for (i=maxpvt; i<m; i++)
    SET_VECTOR_ELT(ret, i, mkChar(""));
  
  free(buf);
  
  R_END;
  return ret;
}

