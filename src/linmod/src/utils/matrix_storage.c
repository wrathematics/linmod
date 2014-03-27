/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, Schmidt


#include <stdio.h>
#include <stdlib.h>

#define BLOCKSIZE 8

static void transpose_square_naive(int m, int n, double *x) // m == n
{
  int i, j;
  int ind, tind;
  double tmp;


  #pragma omp parallel for default(shared) private(i, j, tmp) if(n*n>5000)
  for (i=0; i<n; i++)
  {
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


void transpose_square(int m, int n, double *x)
{
  //
  if (m*n < 20000)
    transpose_square_naive(m, n, x);


}


// Binary GCD implementation; reference is one of Lehmer's papers I think, don't remember
static int gcd_binary(int a, int b)
{
  int k = 0;
  int tmp;


  // Trivial cases
  if (a == 0)
    return abs(b);

  if (b == 0)
    return abs(a);

  a = abs(a);
  b = abs(b);

  if (a == 1 || b == 1)
    return 1;

  if ( a==b )
    return a;


  // General case

  // Remove greatest common power k of 2
  while ( ((a | b) & 1) == 0 )
  {
    a >>= 1;
    b >>= 1;
    k++;
  }

  // Remove any remaining power of 2
  if ((a & 1)==0)
    a >>= 1;
  else if ((b & 1)==0)
    b >>= 1;

  while (a != 0 && a != b)
  {
    while ((a & 1)==0)
      a >>= 1;

    // gcd(a, b) = gcd( (a-b)/2, b )
    if (a < b)
    {
      tmp = a;
      a = b;
      b = tmp;
    }

    a = a-b;
  }

  // gcd = b*2^k
  return b << k;
}

// -----------------------------------------------------
/*
  Rectangular (non-square) transpose, using a max(m, n)
  buffer.  Based on:
    Catanzaro, B., Keller, A., and Garland, M. "A Decomposition
    for In-place Matrix Transposition".
*/
// -----------------------------------------------------

static inline int eqn_23(int i, int j, int b, int m)
{
  return (i + ((int) j/b)) % m;
}

static inline int eqn_24(int i, int j, int b, int m, int n)
{
  return ((i + ((int) j/b)) % m + j*m) % n;
}

static inline int eqn_26(int i, int j, int a, int m, int n)
{
  return (j + i*n - ((int) i/a)) % m;
}

// column to row index
static void c2r_work(const int m, const int n, double *x)
{
  int i, j;

  const int gcd = gcd_binary(m, n);
  const int a = m/gcd;
  const int b = n/gcd;

  const int mn = (m>n ? m : n);

  double *tmp;
  tmp = malloc(mn * sizeof(*tmp));


  #pragma omp parallel private(i, j) if(m*n>5000)

  // Rotate columns
  if (gcd > 1)
  {
    #pragma omp for
    for (i=0; i<n; i++)
    {
      for (j=0; j<m; j++)
        tmp[i] = x[eqn_23(i, j, b, m) + m*j]; // Gather

      for (i=0; i<m; i++)
        x[i + m*j] = tmp[i];
    }
    #pragma omp end for
  }

  // Row shuffle
  #pragma omp for
  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
      tmp[eqn_24(i, j, b, m, n)] = x[i + m*j]; // Scatter

    for (j=0; j<n; j++)
      x[i + m*j] = tmp[j];
  }
  #pragma omp end for

  // Column shuffle
  #pragma omp for
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
      tmp[i] = x[eqn_26(i, j, a, m, n) + m*j]; // Gather

    for (i=0; i<m; i++)
      x[i + m*j] = tmp[i];
  }
  #pragma omp end for

  #pragma omp end parallel

  free(tmp);
}


void c2r(const int m, const int n, double *x)
{
  if (m == 1 || n == 1)
    return;
  else if (m == n)
    transpose_square(m, n, x);
  else
    c2r_work(m, n, x);
}







// column to row
static void r2c_work(const int m, const int n, double *x)
{
  int i, j;

  const int gcd = gcd_binary(m, n);
  const int a = m/gcd;
  const int b = n/gcd;

  const int mn = (m>n ? m : n);

  double *tmp;
  tmp = malloc(mn * sizeof(*tmp));


  #pragma omp parallel private(i, j) if(m*n>5000)

  // Column shuffle
  #pragma omp for
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
      tmp[i] = x[n*eqn_26(i, j, a, m, n) + j]; // Scatter

    for (i=0; i<m; i++)
      x[n*i + j] = tmp[i];
  }
  #pragma omp end for

  // Row shuffle
  #pragma omp for
  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
      tmp[eqn_24(i, j, b, m, n)] = x[n*i + j]; // Gather

    for (j=0; j<n; j++)
      x[i + m*j] = tmp[j];
  }
  #pragma omp end for

  // Rotate columns
  if (gcd > 1)
  {
    #pragma omp for
    for (i=0; i<n; i++)
    {
      for (j=0; j<m; j++)
        tmp[i] = x[eqn_23(i, j, b, m) + m*j]; // Scatter

      for (i=0; i<m; i++)
        x[n*i + j] = tmp[i];
    }
    #pragma omp end for
  }

  #pragma omp end parallel

  free(tmp);
}







void r2c(const int m, const int n, double *x)
{
  if (m == 1 || n == 1)
    return;
  else if (m == n)
    transpose_square(m, n, x);
  else
    r2c_work(m, n, x);
}










// -----------------------------------------------




void show_matrix(double *x, int w, int h)
{
  int i, j;
  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++)
        printf("%2g ", x[i * w + j]);
    putchar('\n');
  }
}

void rect_test()
{
  int i;
  double m[15];
  for (i = 0; i < 15; i++) m[i] = i + 1;

  puts("Row major:");
  show_matrix(m, 3, 5);

  c2r(3, 5, m);

  puts("\nCol major:");
  show_matrix(m, 3, 5);

  r2c(3, 5, m);

  puts("\nRow major:");
  show_matrix(m, 3, 5);
}

void square_test()
{
  int i;
  double m[16];
  for (i = 0; i < 16; i++) m[i] = i + 1;

  puts("Row major:");
  show_matrix(m, 4, 4);

  c2r(4, 4, m);

  puts("\nCol major:");
  show_matrix(m, 4, 4);
}

int main(void)
{
  rect_test();

/*  square_test();*/

  return 0;
}
