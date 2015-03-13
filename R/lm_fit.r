set_lm_names <- function()
{
  
}



#' Linear Model Fitter
#' 
#' A basic linear model fitter.
#' 
#' The outputs are identical to those of R's \code{lm.fit()}, though the
#' internals differe significantly.
#' 
#' By default, \code{lm_fit()} will behave exactly as R's \code{lm.fit()}
#' though with a very different backend.  R uses the very old and deprecated
#' LINPACK (library, not the benchmark) routines for fitting a linear model,
#' whereas linmod uses more recent LAPACK routines.  The former exclusively
#' uses level 1 BLAS, while the latter makes use of level 3 BLAS, which allows
#' for better performance by better utilizing CPU cache.
#' 
#' @param x 
#' The input data matrix.
#' @param y 
#' The vector/matrix of independent variable(s).
#' @param offset 
#' A vector to be included in the predictors.
#' @param tol 
#' Numerical tolerance for the QR.
#' @param singular.ok 
#' logical; if \code{FALSE}, then the function will error if
#' a singular (rank-degenerate) model is detected.
#' @param check.rank 
#' logical; if \code{TRUE}, then the rank-revealing algorithm
#' will be used.  Otherwise the model is assumed to be full rank.
#' @param ... 
#' Extra arguments; ignored with a warning, as in \code{lm.fit()}
#' 
#' @return 
#' A list containing the elements: \tabular{l}{ coefficients \cr
#' residuals \cr effects \cr rank \cr fitted.values \cr assign \cr qr \cr
#' df.residual \cr } where \code{qr} is a list consisting of the elements:
#' \tabular{l}{ qr \cr qraux \cr pivot \cr tol \cr rank \cr }
#' 
#' @examples
#' \dontrun{
#' library(linmod)
#' 
#' n <- 10
#' p <- 3
#' x <- matrix(rnorm(n*p), n, p)
#' y <- rnorm(n)
#' 
#' lm_fit(x, y)
#' }
#' 
#' @export
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
  
  hasnames <- !is.null(colnames(x))
  
  fit <- .Call(R_LM_FIT, x, y, offset=offset, tol, as.integer(singular.ok), check.rank, hasnames)
  attr(fit$qr, "class") <- "qr"
  
  if (!is.null(rownames(x)))
  {
    rownames(fit$qr$qr) <- rownames(x)
    names(fit$residuals) <- rownames(x)
    names(fit$fitted.values) <- rownames(x)
  }
  
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

