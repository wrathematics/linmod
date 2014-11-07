library(linmod)


printlmtest <- function(nrhs, rank, test)
{
  printtest <- if (test) "TRUE" else "FALSE"
  cat(sprintf("%s:  nrhs=%d, rank=%d\n", printtest, nrhs, rank))
}

bigtest <- function(m, n, seed, verbose, check.rank)
{
  NRHS <- 1:2
  RUNS <- 1:5
  offsets <- list(NULL, rep(1, m))
  
  for (offset in offsets)
  {
    for (nrhs in NRHS)
    {
      ### Assume full rank
      set.seed(seed)
      x <- matrix(rnorm(m*n), m, n)
      y <- matrix(rnorm(m*nrhs), m, nrhs)
      test <- linmod:::lmfit_test(x=x, y=y, offset=offset, check.rank=FALSE, verbose=verbose)
      printlmtest(nrhs, -1, test)
      
      
      ### Don't (and might not be)
      for (run in RUNS)
      {
        if (run == 2)       x[, 5] <- x[, 4]
        else if (run == 3)  x[, 3] <- x[, 2]
        else if (run == 4)  x[, 4] <- x[, 3] <- x[, 2]
        else if (run == 5)  x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
  #      else if (run == 6) x[, 5] <- x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
        
        test <- linmod:::lmfit_test(x=x, y=y, offset=offset, check.rank=TRUE, verbose=verbose)
        printlmtest(nrhs, n-run+1, test)
      }
    }
  }
  
  invisible()
}


m <- 10
n <- 5

#n <- 10
#m <- 5

bigtest(m=m, n=n, seed=1234, verbose=FALSE)
