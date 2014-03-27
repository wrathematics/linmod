library(linmod)

set.seed(1234)
m <- 10
n <- 3
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)

#x[, 3] <- x[, 2]

#stats::lm.fit(x, y)$coefficients
stats::lm.fit(x, y)$fitted
cat("-------------------------------------------------------\n\n")

#lm_fit(x, y)
#cat("-------------------------------------------------------\n\n")

#linmod::lm_fit(x, y)$coefficients[1:3]

#linmod::lm_fit_R(x, y)$coefficients
linmod::lm_fit_R(x, y)$fitted




#m <- 10000
#n <- 300
#x <- matrix(rnorm(m*n), m, n)
#y <- rnorm(m)

#system.time(lm.fit(x, y))[3]

#system.time(lm_fit_R(x, y))[3]
#system.time(lm_fit_R(x, y, checkrank=FALSE))[3]
