library(linmod)
library(RcppArmadillo)
library(rbenchmark)

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
###t1 <- system.time(mdl1 <- fastLm(X=x, y=y, method=1))[3]
###cat(paste("fastLm:", round(t1, 3), "\n"))
###t2 <- system.time(mdl2 <- lm_fit(x, y, check.rank=FALSE))[3]
###cat(paste("linmod:", round(t2, 3), "\n"))
###t3 <- system.time(mdl3 <- lm.fit(x, y))[3]
###cat(paste("R core:", round(t3, 3), "\n"))

###all.equal(mdl2, mdl3)

benchmark(mdl1 <- fastLm(X=x,  y=y, method=0), 
          mdl2 <- lm_fit(x=x, y=y, check.rank=TRUE),
          mdl3 <- lm.fit(x=x, y=y),
          replications=10,
          columns=c("test", "replications", "elapsed", "relative")
)

