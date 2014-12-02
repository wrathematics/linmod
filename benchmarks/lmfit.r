library(linmod)
#library(RcppArmadillo)
library(RcppEigen)
library(rbenchmark)


burnin <- function(reps=10)
{
  x <- matrix(rnorm(30), 10)
  y <- rnorm(10)
  
  replicate(fastLm(X=x, y=y), n=reps)
  replicate(lm_fit(x=x, y=y), n=reps)
  replicate(lm.fit(x=x, y=y), n=reps)
  
  invisible()
}
burnin()


reps <- 10

m <- 2500
n <- 250

x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)



cat("------------------ RRQR ------------------\n")

benchmark(
#          fastLm(X=x,  y=y, method=0), 
          lm_fit(x=x, y=y, check.rank=TRUE),
          lm.fit(x=x, y=y),
          replications=reps,
          columns=c("test", "replications", "elapsed", "relative")
)


