/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, Schmidt


#include <stdlib.h>
#include "lapack.h"

int lmfit_classic(char trans, int m, int n, int nrhs, double *a, double *b)
{
  int info = 0;

  double *work, tmpwork;
  int lwork = -1;

  // Workspace query
  dgels_(&trans, &m, &n, &nrhs, a, &m, b, &m, &tmpwork, &lwork, &info);

  lwork = (int) tmpwork;
  work = malloc(lwork * sizeof(*work));

  // Fit y~x
  dgels_(&trans, &m, &n, &nrhs, a, &m, b, &m, work, &lwork, &info);

  return info;
}



int lmfit_R(char trans, int m, int n, int nrhs, double *a, double *b,
  double tol, double **eff, double **ft, double **rsd, double **tau,
  int rank)
{
  int info = 0;

  double *work, tmpwork;
  int lwork = -1;


  // Workspace query
  dgels_(&trans, &m, &n, &nrhs, a, &m, b, &m, &tmpwork, &lwork, &info);

  lwork = (int) tmpwork;
  work = malloc(lwork * sizeof(*work));

  *eff = malloc(m*nrhs * sizeof(**eff));
  *ft = malloc(m*nrhs * sizeof(**ft));
  *rsd = malloc(m*nrhs * sizeof(**rsd));
  *tau = malloc(m*nrhs * sizeof(**tau));


  // Fit y~x
  rdgels_(&trans, &m, &n, &nrhs, a, &m, b, &m, work, &lwork, &info, &tol, *eff, *ft, *rsd, *tau, &rank);

  return info;
}

