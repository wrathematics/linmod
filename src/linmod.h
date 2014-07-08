#ifndef __LINMOD_H__
#define __LINMOD_H__


#include <RNACI.h> 
#include "linmod/src/c_interface/linmod.h"
#include "linmod/src/c_interface/lapack.h"

SEXP make_lmfit_default_rownames(const int n);


#endif
