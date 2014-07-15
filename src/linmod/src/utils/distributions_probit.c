/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, Schmidt

// qnorm implementation is from "An algorithm for computing the
// inverse normal cumulative distribution function", 
// http://home.online.no/~pjacklam/notes/invnorm/#FAQ
// This implementation is reasonably quick, and has relative error
// (in norm) less than 1.15e-9.  "Special cases" (dumb inputs) 
// conform to R's qnorm().

#include <math.h>
#include <stdlib.h>


#define P_LOW 0.02425
#define P_HIGH 0.97575

#define ROOT2 1.414213562373095


// Coefficients in rational approximations.
static const double a[6] =
{
  -3.969683028665376e+01, 2.209460984245205e+02, 
  -2.759285104469687e+02, 1.383577518672690e+02, 
  -3.066479806614716e+01, 2.506628277459239e+00
};

static const double b[5] =
{
  -5.447609879822406e+01, 1.615858368580409e+02,
  -1.556989798598866e+02, 6.680131188771972e+01, 
  -1.328068155288572e+01
};

static const double c[6] =
{
  -7.784894002430293e-03, -3.223964580411365e-01, 
  -2.400758277161838e+00, -2.549732539343734e+00, 
  4.374664141464968e+00, 2.938163982698783e+00
};

static const double d[4] =
{
  7.784695709041462e-03, 3.224671290700398e-01, 
  2.445134137142996e+00, 3.754408661907416e+00
};

double qnorm(const double p)
{
  int i;
  double q, r;
  double numerator, denominator;
  
  // Special cases
  if (p < 0 || p > 1)
    return NAN;
  else if (p == 0)
    return -INFINITY;
  else if (p == 1)
    return INFINITY;
  // Rational approximation for lower region
  else if (p < P_LOW)
  {
    q = sqrt(-2.*log(p));
    
    numerator = ((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5];
    denominator = (((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1;
  }
  // Rational approximation for central region
  else if (p < P_HIGH)
  {
    q = p - 0.5;
    r = q*q;
    
    numerator = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q;
    denominator = (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }
  // Rational approximation for upper region
  else
  {
    q = sqrt(-2.*log(1. - p));
    
    numerator = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]);
    denominator = ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
  
  
  return numerator/denominator;
}



double erfinv(const double x)
{
  return qnorm((1. + x)/2.) / ROOT2;
}



double erfinvc(const double x)
{
  return erfinv(1. - x);
}



double probit(const double x)
{
  return -ROOT2 * erfinvc(2. * x);
}

