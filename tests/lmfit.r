library(linmod)

set.seed(1234)
nrhs <- 1

m <- 10
n <- 5

#n <- 10
#m <- 5

x <- matrix(rnorm(m*n), m, n)
y <- matrix(rnorm(m*nrhs), m, nrhs)

#x[, 5] <- x[, 4]
#x[, 3] <- x[, 2]
#x[, 4] <- x[, 3] <- x[, 2]
#x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]
#x[, 5] <- x[, 4] <- x[, 3] <- x[, 2] <- x[, 1]


#stats::lm.fit(x, y)$qr
#linmod::lm_fit(x, y)$qr



f <- function(x, y, check.rank=FALSE)
{
  cat("-------------------------------------------------------\n\n")
  
  mdl1 <- stats::lm.fit(x, y)
  mdl2 <- linmod::lm_fit(x, y, check.rank=check.rank)
  
  a <- mdl1$coefficients
  b <- mdl2$coefficients
  test <- all.equal(a, b, check.names=FALSE)
  
  print(test)
  print(a);print(b)
  
  
  a <- mdl1$fitted
  b <- mdl2$fitted
  test <- all.equal(a, b, check.names=FALSE)
  print(test)
  
  a <- mdl1$residuals
  b <- mdl2$residuals
  test <- all.equal(a, b, check.names=FALSE)
  print(test)
  
  a <- mdl1$effects
  b <- mdl2$effects
  test <- all.equal(a, b, check.names=FALSE)
  print(test)
  
  invisible()
}

f(x, y, TRUE)
#f(x, y, FALSE)


#stats::lm.fit(x, y)
#linmod::lm_fit(x, y)
