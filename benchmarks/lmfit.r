library(linmod)
library(RcppEigen)
library(rbenchmark)


burnin <- function(reps=10)
{
  x <- matrix(rnorm(30), 10)
  y <- rnorm(10)
  
  replicate(fastLm(X=x, y=y), n=reps)
  replicate(lm_fit(x=x, y=y), n=reps)
  
  invisible()
}
burnin()

m <- 4000
n <- 250
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)


#cat("----------------- lm.fit -----------------\n")
#t1 <- system.time(lm.fit(x=x, y=y))[3]
#cat(paste("lm.fit:", round(t1, 3), "\n"))
#cat("\n")


cat("------------------ RRQR ------------------\n")
t1 <- system.time(mdl1 <- fastLm(X=x, y=y, method=0))[3]
cat(paste("fastLm:", round(t1, 3), "\n"))
t2 <- system.time(mdl2 <- lm_fit(x, y))[3]
cat(paste("linmod:", round(t2, 3), "\n"))
t3 <- system.time(mdl3 <- lm.fit(x, y))[3]
cat(paste("R core:", round(t3, 3), "\n"))

all.equal(mdl1$coefficients, mdl2$coefficients, check.attributes=FALSE)


#cat("\n")


#cat("------------ Assume full rank ------------\n")
#t1 <- system.time(fastLm(X=x, y=y, method=1))[3]
#t2 <- system.time(lm_fit(x, y, check.rank=FALSE))[3]
#cat(paste("fastLm:", round(t1, 3), "\n"))
#cat(paste("linmod:", round(t2, 3), "\n"))








#m <- 10000
#n <- 300
#x <- matrix(rnorm(m*n), m, n)
#y <- rnorm(m)



#cols <- c("test", "replications", "elapsed")


#benchmark(fastLm(X=x, y=y, method=0), 
#          lm_fit(x, y, check.rank=TRUE),
#          replications=10, columns=cols)


#benchmark(fastLm(X=x, y=y, method=1), 
#          lm_fit(x, y, check.rank=FALSE),
#          replications=10, columns=cols)
