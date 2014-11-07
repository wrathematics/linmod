// Copyright 2014, Schmidt

#include <string.h>
#include "linmod.h"
#include "linmod/src/c_interface/glm_defines.h"


static int glm_fit_family_val(char *family)
{
  if (strcmp(family, "gaussian") == 0)
    return GLM_FAMILY_GAUSSIAN;
  else if (strcmp(family, "binomial") == 0)
    return GLM_FAMILY_BINOMIAL;
  else if (strcmp(family, "poisson") == 0)
    return GLM_FAMILY_POISSON;
  else if (strcmp(family, "gamma") == 0)
    return GLM_FAMILY_GAMMA;
  else if (strcmp(family, "inverse.gaussian") == 0)
    return GLM_FAMILY_INVERSEGAUSSIAN;
  else
    return GLM_FAMILY_UNSUPPORTED;
}

static int glm_fit_link_val(char *link)
{
  if (strcmp(link, "cloglog") == 0)
    return GLM_LINK_CLOGLOG;
  else if (strcmp(link, "identity") == 0)
    return GLM_LINK_IDENTITY;
  else if (strcmp(link, "inverse") == 0)
    return GLM_LINK_INVERSE;
  else if (strcmp(link, "log") == 0)
    return GLM_LINK_LOG;
  else if (strcmp(link, "logit") == 0)
    return GLM_LINK_LOGIT;
  else if (strcmp(link, "sqrt") == 0)
    return GLM_LINK_SQRT;
  else if (strcmp(link, "probit") == 0)
    return GLM_LINK_PROBIT;
  else if (strcmp(link, "cauchit") == 0)
    return GLM_LINK_CAUCHIT;
  else if (strcmp(link, "1/mu^2") == 0)
    return GLM_LINK_INVERSESQUARE;
  else
    return GLM_LINK_UNSUPPORTED;
}

static int glm_fit_stoprule_val(char *stoprule)
{
  if (strcmp(stoprule, "maxiter") == 0)
    return GLM_STOPRULE_MAXITER;
  else if (strcmp(stoprule, "coefficients") == 0)
    return GLM_STOPRULE_COEFS;
  else if (strcmp(stoprule, "deviance") == 0)
    return GLM_STOPRULE_DEVIANCE;
  else
    return GLM_BADINPUT_STOPRULE;
}

SEXP R_glm_fit_family_val(SEXP family)
{
  R_INIT;
  SEXP ret;
  newRvec(ret, 1, "int");
  
  INT(ret) = glm_fit_family_val(STR(family));
  
  R_END;
  return ret;
}

SEXP R_glm_fit_link_val(SEXP link)
{
  R_INIT;
  SEXP ret;
  newRvec(ret, 1, "int");
  
  INT(ret) = glm_fit_link_val(STR(link));
  
  R_END;
  return ret;
}

SEXP R_glm_fit_stoprule_val(SEXP stoprule)
{
  R_INIT;
  SEXP ret;
  newRvec(ret, 1, "int");
  
  INT(ret) = glm_fit_stoprule_val(STR(stoprule));
  
  R_END;
  return ret;
}

