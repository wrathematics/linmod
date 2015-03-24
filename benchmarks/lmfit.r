library(memuse)
library(linmod)
#library(RcppArmadillo)
library(RcppEigen)
library(rbenchmark)
library(microbenchmark)


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


reps <- 15

m <- 12000
n <- 250

x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)



cat("------------------ RRQR ------------------\n")
cat(paste0("Data size:  ", object.size(x)+object.size(y), "\n"))
cat(paste0("L3 cache size:  ", Sys.cachesize()$L3, "\n"))


#bench <- benchmark(
#          lm.fit(x=x, y=y),
#          lm_fit(x=x, y=y, check.rank=TRUE),
##          RcppEigen::fastLm(X=x,  y=y, method=0), 
#          replications=reps,
#          columns=c("test", "replications", "elapsed", "relative")
#)


control <- list(warmup=20, order="random")

bench <- microbenchmark(
          lm.fit(x=x, y=y),
          lm_fit(x=x, y=y, check.rank=TRUE),
#          RcppEigen::fastLm(X=x,  y=y, method=0), 
          times=reps,
          unit="s",
          control=control
)

print(bench)

