library(linmod)

set.seed(1234)
m <- 10
n <- 4
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)

x[, 3] <- x[, 2]


#stats::lm.fit(x, y)$qr
#linmod::lm_fit(x, y)$qr



f <- function(x, y, check.rank=TRUE){
  cat("-------------------------------------------------------\n\n")
  
  a <- stats::lm.fit(x, y)$coefficients
  b <- linmod::lm_fit(x, y)$coefficients
  test <- all.equal(a, b, check.names=FALSE)
  print(test)
  print(a)
  print(b)
  
#  a <- stats::lm.fit(x, y)$fitted
#  b <- linmod::lm_fit(x, y)$fitted
#  test <- all.equal(a, b, check.names=FALSE)
#  print(test)
#  
#  a <- stats::lm.fit(x, y)$residuals
#  b <- linmod::lm_fit(x, y)$residuals
#  test <- all.equal(a, b, check.names=FALSE)
#  print(test)
#  
#  a <- stats::lm.fit(x, y)$effects
#  b <- linmod::lm_fit(x, y)$effects
#  test <- all.equal(a, b, check.names=FALSE)
#  print(test)
#  print(a)
#  print(b)
  
  invisible()
}

f(x, y, TRUE)
#f(x, y, FALSE)


#stats::lm.fit(x, y)
#linmod::lm_fit(x, y)
