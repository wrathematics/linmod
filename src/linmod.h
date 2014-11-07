#ifndef __LINMOD_H__
#define __LINMOD_H__


#include <RNACI.h> 
#include <stdbool.h>
#include "linmod/src/c_interface/linmod.h"
#include "linmod/src/c_interface/lapack.h"


// utils.c
SEXP make_lmfit_default_rownames(const int n);
SEXP make_lmfit_default_effectnames(const int m, const int n, const int rank, const int *pvt);


#endif
