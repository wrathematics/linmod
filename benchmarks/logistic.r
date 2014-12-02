library(linmod, quietly=TRUE)

n <- 5000
p <- 300

family <- binomial(logit)

linmod:::glm_test(n, p, family, verbose=TRUE, timings=TRUE)
