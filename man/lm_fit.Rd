\name{lm_fit}
\alias{lm_fit}

\title{Linear Model Fitter}

\description{
  A basic linear model fitter.  
  
  The outputs are identical to those of R's \code{lm.fit()}, though
  the internals differe significantly.
}
\usage{
lm_fit(x, y, offset=NULL, tol=1e-07, singular.ok=TRUE, check.rank=TRUE, ...)
}
\arguments{
  \item{x}{The input data matrix.}
  \item{y}{The vector/matrix of independent variable(s).}
  \item{offset}{A vector to be included in the predictors.}
  \item{tol}{Numerical tolerance for the QR.}
  \item{singular.ok}{logical; if \code{FALSE}, then the function will
  error if a singular (rank-degenerate) model is detected.}
  \item{check.rank}{logical; if \code{TRUE}, then the rank-revealing
  algorithm will be used.  Otherwise the model is assumed to be full
  rank.}
  \item{...}{Extra arguments; ignored with a warning, as in 
  \code{lm.fit()}}
}
\details{
  By default, \code{lm_fit()} will behave exactly as R's \code{lm.fit()}
  though with a very different backend.  R uses the very old and
  deprecated LINPACK (library, not the benchmark) routines for 
  fitting a linear model, whereas linmod uses more recent LAPACK
  routines.  The former exclusively uses level 1 BLAS, while
  the latter makes use of level 3 BLAS, which allows for better
  performance by better utilizing CPU cache.
}
\value{
  A list containing the elements:
\tabular{l}{
  coefficients \cr
  residuals \cr
  effects \cr
  rank \cr
  fitted.values \cr
  assign \cr
  qr \cr
  df.residual \cr
}
  where \code{qr} is a list consisting of the elements:
\tabular{l}{
  qr \cr
  qraux \cr
  pivot \cr
  tol \cr
  rank \cr
}
}
\examples{
\dontrun{
library(linmod)

n <- 10
p <- 3
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)

lm_fit(x, y)
}
}
