library(linmod)

set.seed(1234)
nrhs <- 2



test <- function(m, n, nrhs, offset=NULL, whichprint=NULL)
{
  x <- matrix(rnorm(m*n), m, n)
  y <- matrix(rnorm(m*nrhs), m, nrhs)
  
  linmod:::lmfit_test(x, y, offset=NULL, TRUE, whichprint)
}

#x[, 5] <- x[, 4]
#x[, 3] <- x[, 2]
#x[, 4] <- x[, 3] <- x[, 2]
#x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
#x[, 5] <- x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]



test(m=10, n=5, nrhs=2, offset=NULL, whichprint="coefficients")
test(m=2000, n=100, nrhs=2, offset=NULL, whichprint=NULL)

