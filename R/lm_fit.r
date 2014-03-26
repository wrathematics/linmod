lm_fit <- function(x, y)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  
  fit <- .Call("R_LM_FIT", x, y, PACKAGE="linmod")
  
  return(fit)
}


lm_fit_R <- function(x, y, tol=1e-7)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  
  fit <- .Call("R_LM_FIT_R", x, y, tol, PACKAGE="linmod")
  
  return(fit)
}
