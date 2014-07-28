library(linmod, quiet=TRUE)

### TODO use canonical link to generate starting values for the other bullshits

#library(R330)
#data(chd.df)

#x <- cbind(rep(1, nrow(chd.df)), matrix(chd.df[,1]))
#y <- matrix(chd.df[,2])

set.seed(1234)

n <- 10
p <- 3

offset <- rep(1, n)
control <- list(trace=FALSE)


### Binomial
links <- c("cloglog", "log", "logit", "probit", "cauchit")
links <- "logit"
#links <- "probit"
#links <- "cauchit"
families <- lapply(links, binomial)
invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE, offset=offset, control=control))


### Gamma
#links <- c("identity", "inverse", "log")
##links <- "log"
#families <- lapply(links, Gamma)
#invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE))
#glm_test(n=n, p=p, verbose=TRUE, timings=FALSE, Gamma(log))


### Gaussian
#links <- c("identity", "inverse", "log")
#families <- lapply(links, gaussian)
#invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE, offset=offset))


### Poisson
#links <- c("identity", "log", "sqrt")
#families <- lapply(links, poisson)
#invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE))


### Inverse Gaussian
#links <- c("inverse", "log", "identity", "1/mu^2")
#families <- lapply(links, inverse.gaussian)
#invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE))



