// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2014, Schmidt


#ifndef __GLM_GLM_DEFINES_H_
#define __GLM_GLM_DEFINES_H_


// Bad inputs and errors
#define GLM_oom = -2147483647

#define GLM_BADINPUT_STOPRULE   -4
#define GLM_BADINPUT_N          -5
#define GLM_BADINPUT_P          -6
#define GLM_BADINPUT_FAMILY     -8
#define GLM_BADINPUT_MAXITER   -13
#define GLM_BADINPUT_TOL       -14

// Families
#define GLM_FAMILY_UNSUPPORTED     -1
#define GLM_FAMILY_GAUSSIAN         1
#define GLM_FAMILY_BINOMIAL         2
#define GLM_FAMILY_POISSON          3
#define GLM_FAMILY_GAMMA            4
#define GLM_FAMILY_INVERSEGAUSSIAN  5

#define GLM_FAMILY_BADRESPONSE -8

// Links
#define GLM_LINK_UNSUPPORTED   -2
#define GLM_LINK_CLOGLOG        1
#define GLM_LINK_IDENTITY       2
#define GLM_LINK_INVERSE        3
#define GLM_LINK_LOG            4
#define GLM_LINK_LOGIT          5
#define GLM_LINK_SQRT           6
#define GLM_LINK_PROBIT         7
#define GLM_LINK_CAUCHIT        8
#define GLM_LINK_INVERSESQUARE  9

// Intercept
#define GLM_INTERCEPT_NO  0
#define GLM_INTERCEPT_YES 1

// Stoprules
#define GLM_STOPRULE_MAXITER  1
#define GLM_STOPRULE_COEFS    2
#define GLM_STOPRULE_DEVIANCE 3

// Convergence
#define GLM_CONVERGENCE_NOCONVERGENCE  -1
#define GLM_CONVERGENCE_CONVERGED       1
#define GLM_CONVERGENCE_INFPARAMS       2
#define GLM_CONVERGENCE_NOCHANGE        3

#endif
