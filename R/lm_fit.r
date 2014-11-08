lm_fit <- function(x, y, offset=NULL, tol=1e-07, singular.ok=TRUE, check.rank=TRUE, ...)
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
  
  fit <- .Call(R_LM_FIT, x, y, offset=offset, tol, as.integer(singular.ok), check.rank)
  attr(fit$qr, "class") <- "qr"
  
  if (!is.null(rownames(x)))
  {
    rownames(fit$qr$qr) <- rownames(x)
    names(fit$residuals) <- rownames(x)
    names(fit$fitted.values) <- rownames(x)
  }
  
  ### FIXME don't fit the other ones in the first place...
  if (!is.null(colnames(x)))
  {
    colnames(fit$qr$qr) <- colnames(x)
    names(fit$coefficients) <- colnames(x)
    names(fit$effects) <- colnames(x)
  }
  
  if (!is.null(attr(x, "assign")))
  {
    attr(fit$qr$qr, "assign") <- attr(x, "assign")
    fit$assign <- attr(x, "assign")
  }
  
  return(fit)
}



### FIXME this isn't quite right, check the return
.lm_fit <- function(x, y, tol=1e-07, check.rank=TRUE)
{
  lm_fit(x=x, y=y, offset=NULL, tol=tol, singular.ok=TRUE, check.rank=check.rank)
}

