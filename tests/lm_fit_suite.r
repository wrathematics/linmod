library(linmod)


printlmtest <- function(nrhs, rank, test, offset)
{
  cat(sprintf("%s:\tnrhs=%d,\t rank=%d,\t offset=%s\n", as.character(test), nrhs, rank, as.character(offset)))
}

bigtest <- function(m, n, seed, verbose, check.rank)
{
  cat(sprintf("#-------------------------------------------#\nm=%d n=%d\n", m, n))
  NRHS <- 1:2
  RUNS <- 1:5
  offsets <- list(NULL, rep(1, m))
  
  for (nrhs in NRHS)
  {
    for (offset in offsets)
    {
      ### Assume full rank
      set.seed(seed)
      x <- matrix(rnorm(m*n), m, n)
      y <- matrix(rnorm(m*nrhs), m, nrhs)
      test <- linmod:::lmfit_test(x=x, y=y, offset=offset, check.rank=FALSE, verbose=verbose)
      printlmtest(nrhs, -1, test, !is.null(offset))
      
      
      ### Don't (and might not be)
      for (run in RUNS)
      {
        if (run == 2)       x[, 5] <- x[, 4]
        else if (run == 3)  x[, 3] <- x[, 2]
        else if (run == 4)  x[, 4] <- x[, 3] <- x[, 2]
        else if (run == 5)  x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
  #      else if (run == 6) x[, 5] <- x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
        
        test <- linmod:::lmfit_test(x=x, y=y, offset=offset, check.rank=TRUE, verbose=verbose)
        printlmtest(nrhs, n-run+1, test, !is.null(offset))
      }
    }
  }
  
  invisible()
}


## TODO
#n <- 10
#m <- 5

bigtest(m=10, n=5, seed=1234, verbose=FALSE)
bigtest(m=10, n=10, seed=1234, verbose=FALSE)
