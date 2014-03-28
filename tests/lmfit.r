library(linmod)

set.seed(1234)
m <- 10
n <- 3
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)

#x[, 3] <- x[, 2]

#stats::lm.fit(x, y)$coefficients
#stats::lm.fit(x, y)$fitted
stats::lm.fit(x, y)$residuals
#stats::lm.fit(x, y)$effects
#stats::lm.fit(x, y)$qr

cat("-------------------------------------------------------\n\n")

#lm_fit(x, y)
#cat("-------------------------------------------------------\n\n")

#linmod::lm_fit(x, y)$coefficients[1:3]

#linmod::lm_fit(x, y)$coefficients
#linmod::lm_fit(x, y)$fitted
linmod::lm_fit(x, y)$residuals
#linmod::lm_fit(x, y)$effects
#linmod::lm_fit(x, y)$qr




m <- 10000
n <- 300
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)

system.time(lm.fit(x, y))[3]

system.time(lm_fit(x, y))[3]
system.time(lm_fit(x, y, checkrank=FALSE))[3]
