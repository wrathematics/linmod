library(glm, quiet=T)



set.seed(1684)
#set.seed(12334)  # works



n <- 10
p <- 3

x <- matrix(rnorm(n*p, mean=10, sd=1), n)


#fam <- binomial(logit)
#y <- sample(0:1, size=n, replace=T)

fam <- poisson(log)
y <- rpois(n, 5)

#fam <- gaussian()
#y <- rnorm(n)

#fam <- Gamma()
#y <- rgamma(n, 2)


#library(R330)
#data(chd.df)

#x <- cbind(rep(1, nrow(chd.df)), matrix(chd.df[,1]))
#y <- matrix(chd.df[,2])


intercept <- TRUE
#intercept <- FALSE










system.time({
#  mdl <- glm.fit(x=x, y=y, family=binomial(logit), intercept=intercept)
  mdl <- glm.fit(x=x, y=y, family=fam, intercept=intercept)
})[3]

mdl$coefficients
#mdl$deviance
#mdl$null.deviance
mdl$iter


maxiter <- 20

stoprule <- 3


system.time({
  mdl2 <- glm_fit(x=x, y=y, family=fam, maxiter, intercept=intercept, stoprule=stoprule)
})[3]

mdl2














##n <- 10
##p <- 3

##x <- matrix(rnorm(n*p, mean=10, sd=1), n)


##fam <- binomial(logit)
##y <- sample(0:1, size=n, replace=T)
##fam <- poisson(log)
##y <- rpois(n, 5)


##mdl <- glm.fit(x=x, y=y, family=fam, intercept=T)


