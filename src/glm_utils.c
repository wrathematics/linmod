#include <string.h>
#include "linmod.h"
#include "linmod/src/c_interface/glm_defines.h"


int glm_fit_family(char *family)
{
  if (strcmp(family, "gaussian") == 0)
    return GLM_FAMILY_GAUSSIAN;
  else if (strcmp(family, "binomial") == 0)
    return GLM_FAMILY_BINOMIAL;
  else if (strcmp(family, "poisson") == 0)
    return GLM_FAMILY_POISSON;
  else if (strcmp(family, "gamma") == 0)
    return GLM_FAMILY_GAMMA;
  else
    return GLM_FAMILY_UNSUPPORTED;
}

int glm_fit_link(char *link)
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
  else
    return GLM_LINK_UNSUPPORTED;
}

SEXP R_glm_fit_family(SEXP family)
{
  R_INIT;
  SEXP ret;
  newRvec(ret, 1, "int");
  
  INT(ret) = glm_fit_family(STR(family));
  
  R_END;
  return ret;
}

SEXP R_glm_fit_link(SEXP link)
{
  R_INIT;
  SEXP ret;
  newRvec(ret, 1, "int");
  
  INT(ret) = glm_fit_link(STR(link));
  
  R_END;
  return ret;
}

