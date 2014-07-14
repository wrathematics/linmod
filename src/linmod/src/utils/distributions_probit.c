/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013-2014, Schmidt


#include <math.h>
#include <stdlib.h>

#define ROOT2    1.414213562373095
#define ROOTEPS  1.490116119384766e-8

const double erfinv_coefficients[25] = {
  0.886226925452758, 0.232013666534654, 0.127556175305598,
  0.086552129241548, 0.064959617745385, 0.051731281984616,
  0.042836720651797, 0.036465929308532, 0.031689005021605,
  0.027980632964995, 0.025022275841198, 0.022609863318898,
  0.020606780379059, 0.018918217250779, 0.017476370562857,
  0.016231500987685, 0.015146315063248, 0.014192316002510,
  0.013347364197421, 0.012594004871332, 0.011918295936392,
  0.011308970105923, 0.010756825303318, 0.010254274081853,
  0.009795005770071
};


// inverse of the error function erf; uses taylor series
double erfinv(const double x)
{
  int i;
  double ret, retold;
  double xsq, tmp;
  
  
  if (x == 1) return INFINITY;
  else if (x > 1) return NAN;
  else if (x == 0) return 0; // why not?
  else if (x == -1) return -INFINITY;
  else if (x < -1) return NAN;
  
  xsq = x*x;
  tmp = x;
  
  ret = x*erfinv_coefficients[0];
  
  for (i=1; i<25; i++)
  {
    tmp *= xsq;
    ret += tmp*erfinv_coefficients[i];
    
    if (fabs(ret - retold) < ROOTEPS)
      break;
  }
  
  return ret;
}

void erfinv_(double *x, double *ret)
{
  *ret = erfinv(*x);
}


// FIXME this is probably not very smart...
double erfinvc(const double x)
{
  return erfinv(1. - x);
}

void erfinvc_(double *x, double *ret)
{
  *ret = erfinv(1. - (*x));
}


double probit(const double x)
{
  return -ROOT2 * erfinvc(2. * x);
}

void probit_(double *x, double *ret)
{
  *ret = probit(*x);
}

