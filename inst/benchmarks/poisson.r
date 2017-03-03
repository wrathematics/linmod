library(linmod, quietly=TRUE)

n <- 50000
p <- 300

family <- poisson(log)

glm_test(n, p, family, verbose=TRUE, timings=TRUE)
