library(linmod)

m <- 10
n <- 5

#n <- 10
#m <- 5

bigtest <- function(seed, verbose, check.rank)
{
  NRHS <- 1:2
  RUNS <- 1:5
  
  for (nrhs in NRHS)
  {
    for (run in RUNS)
    {
      set.seed(seed)
      x <- matrix(rnorm(m*n), m, n)
      y <- matrix(rnorm(m*nrhs), m, nrhs)
      
      if (run == 2)       x[, 5] <- x[, 4]
      else if (run == 3)  x[, 3] <- x[, 2]
      else if (run == 4)  x[, 4] <- x[, 3] <- x[, 2]
      else if (run == 5)  x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
#      else if (run == 6) x[, 5] <- x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
      
      print(linmod:::lmfit_test(x=x, y=y, check.rank=check.rank, verbose=verbose))
    }
  }
  
  invisible()
}

bigtest(seed=1234, verbose=FALSE, check.rank=TRUE)
