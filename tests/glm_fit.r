library(linmod, quiet=TRUE)


#library(R330)
#data(chd.df)

#x <- cbind(rep(1, nrow(chd.df)), matrix(chd.df[,1]))
#y <- matrix(chd.df[,2])

set.seed(1234)

n <- 10
p <- 3


#### Binomial
links <- c("cloglog", "log", "logit", "probit", "cauchit")
links <- "logit"
links <- "cloglog"
families <- lapply(links, binomial)
invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE))


#### Gamma
#links <- c("identity", "inverse", "log")
#families <- lapply(links, Gamma)
#invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE))


### Gaussian
#links <- c("identity", "inverse", "log")
#families <- lapply(links, gaussian)
#invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE))


### Poisson
#links <- c("identity", "log", "sqrt")
#families <- lapply(links, poisson)
#invisible(sapply(families, glm_test, n=n, p=p, verbose=TRUE, timings=FALSE))


