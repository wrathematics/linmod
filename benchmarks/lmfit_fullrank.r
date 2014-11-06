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

m <- 10000
n <- 250
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)



cat("------------ Assume full rank ------------\n")

benchmark(fastLm(X=x,  y=y, method=1), 
          lm_fit(x=x, y=y, check.rank=FALSE),
          lm.fit(x=x, y=y),
          replications=10,
          columns=c("test", "replications", "elapsed", "relative")
)

