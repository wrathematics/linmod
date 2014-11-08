library(linmod)

set.seed(1234)
nrhs <- 2

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



#whichprint <- NULL
whichprint <- c("coefficients")

linmod:::lmfit_test(x, y, offset=NULL, TRUE, whichprint)
#f(x, y, FALSE)

#stats::lm.fit(x, y)
#linmod::lm_fit(x, y)
