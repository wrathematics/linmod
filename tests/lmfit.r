library(linmod)

set.seed(1234)
m <- 10
n <- 3
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)

stats::lm.fit(x, y)
cat("-------------------------------------------------------\n\n")

#lm_fit(x, y)
#cat("-------------------------------------------------------\n\n")

glm::lm_fit_R(x, y)




m <- 10000
n <- 30
x <- matrix(rnorm(m*n), m, n)
y <- rnorm(m)

system.time(lm.fit(x, y))[3]

system.time(lm_fit_R(x, y))[3]
