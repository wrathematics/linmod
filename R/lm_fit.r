lm_fit <- function(x, y, tol=1e-7, check.rank=TRUE)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  storage.mode(check.rank) <- "integer"
  
  fit <- .Call(R_LM_FIT, x, y, tol, check.rank)
  
  attr(fit$qr, "class") <- "qr"
  
  return(fit)
}
