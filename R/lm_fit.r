lm_fit <- function(x, y, offset=NULL, tol=1e-7,  check.rank=TRUE, ...)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  storage.mode(check.rank) <- "integer"
  
  if (!is.null(offset)) ### 'numeric of length n'
  {
    if (is.matrix(y)) nrhs <- ncol(y) else nrhs <- 1L
    if (length(offset) != length(y)/nrhs)
      stop("incompatible dimensions")
    
    if (!is.double(offset))
      storage.mode(offset) <- "double"
  }
  
  fit <- .Call(R_LM_FIT, x, y, offset=offset, tol, check.rank)
  
  attr(fit$qr, "class") <- "qr"
  
  return(fit)
}
