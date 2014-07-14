/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, Schmidt


#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
  #include <omp.h>
  #if _OPENMP >= 201307
    #define _OPENMP_SUPPORT_SIMD
  #endif
#endif



static void transpose_square(int m, int n, double *x) // m == n
{
  int i, j;
  int ind, tind;
  double tmp;
  
  
  #pragma omp parallel for default(shared) private(i, j, tmp) if(n*n>5000)
  for (i=0; i<n; i++)
  {
    #if defined(_OPENMP_SUPPORT_SIMD)
    #pragma omp for simd
    #endif
    for (j=0; j<i; j++)
    {
      ind = n*i + j;
      tind = i + n*j;
      
      tmp = x[ind];
      x[ind] = x[tind];
      x[tind] = tmp;
    }
  }
  #pragma omp end parallel
}


