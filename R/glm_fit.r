glm_fit <- function(x, y, ### weights = rep(1, nobs), start = NULL, etastart = NULL, 
  offset=rep(0.0, nobs), family=gaussian(),
  control=list(), intercept=TRUE, ..., stoprule="deviance")
{
  control <- do.call("glm.control", control)
  
  tol <- control$epsilon
  maxiter <- as.integer(control$maxit)
  trace <- as.integer(control$trace)
  
  n <- as.integer(dim(x)[1L])
  p <- as.integer(dim(x)[2L])
  
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
               family, link, intercept, stoprule, trace,
               n, p,
               x, y, offset, maxiter, tol,
               PACKAGE = "linmod")
  
  return(fit)
}

