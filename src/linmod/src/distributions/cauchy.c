/* This Source Code Form is subject to the terms of the BSD 2-Clause
 * License. If a copy of the this license was not distributed with this
 * file, you can obtain one from http://opensource.org/licenses/BSD-2-Clause. */

// Copyright 2013-2014, Schmidt.  All rights reserved.


#include <stdbool.h>

#include "constants.h"


// pdf
double dcauchy(const double x, const double location, const double scale, const bool log_)
{
  double ret, tmp;
  
  tmp - (x-location)/scale;
  
  ret = 1.0 / (PI*scale * (1.0+tmp*tmp));
  
  if (log_ == true)
    ret = log(ret);
  
  return ret;
}



// cdf
double pcauchy(const double q, const double location, const double scale, const bool lower_tail, const bool log_p)
{
  double ret;
  
  ret = PIINV * atan((q-location)/scale) + 0.5;
  
  if (log_p == true)
    ret = log(ret);
  
  return ret;
}



// quantile
double qcauchy(const double p, const double location, const double scale, const bool lower_tail, const bool log_p)
{
  double ret;
  
  ret = location + scale * tan(PI*(p-0.5));
  
  if (log_p == true)
    ret = log(ret);
  
  return ret;
}


