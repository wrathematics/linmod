string_c2f <- function(str, len)
{
  tmp <- len - length(unlist(strsplit(str, split="")))
  
  if (tmp < 0)
    stop("'str' contains more than 'len' chars")
  else if (tmp == 0L)
    return(str)
  else
    return( paste(str, paste(rep(" ", tmp), collapse="", sep=""), collapse="", sep="") )
}


glm_fit <- function(x, y, family, maxiter=25, tol=1e-8, offset=rep(0.0, nobs), intercept=TRUE, stoprule="deviance")
{
  n <- dim(x)[1L]
  p <- dim(x)[2L]
  
  nobs <- NROW(y)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  if (!is.double(offset))
    storage.mode(offset) <- "double"
  
  storage.mode(intercept) <- "integer"
  
  stoprule <- match.arg(tolower(stoprule), c("maxiter", "coefficients", "deviance"))
  stoprule <- .Call("R_glm_fit_stoprule_val", stoprule)
  
  link <- .Call("R_glm_fit_link_val", tolower(family$link))
  family <- .Call("R_glm_fit_family_val", tolower(family$family))
  
  if (stoprule < 0)
    stop("Invalid stoprule")
  if (link < 0)
    stop("Invalid link")
  if (family < 0)
    stop("Invalid family")
  
  fit <- .Call("R_GLM_FIT", 
               family, link, intercept, stoprule, 
               as.integer(n), as.integer(p),
               x, y, offset, as.integer(maxiter), as.double(tol),
               PACKAGE = "linmod")
  
  return(fit)
}

