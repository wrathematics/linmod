lm_fit <- function(x, y, tol=1e-7, checkrank=TRUE)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  storage.mode(checkrank) <- "integer"
  
  fit <- .Call("R_LM_FIT", x, y, tol, checkrank, PACKAGE="linmod")
  
  return(fit)
}
