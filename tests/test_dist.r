library(glm, quiet=T)


#library(R330)
#data(chd.df)

#x <- cbind(rep(1, nrow(chd.df)), matrix(chd.df[,1]))
#y <- matrix(chd.df[,2])






#set.seed(1684)
#set.seed(12334)  # works

n <- 500000
p <- 300

n <- 50000
p <- 300


n <- 50000
p <- 70

x <- matrix(rnorm(n*p, mean=10, sd=1), n)

#fam <- binomial(logit)
#y <- sample(0:1, size=n, replace=T)

fam <- poisson(log)
y <- rpois(n, 5)

#fam <- gaussian()
#y <- rnorm(n)

#fam <- Gamma()
#y <- rgamma(n, 2)



intercept <- TRUE
#intercept <- FALSE


test_dist <- function(x, y, fam, intercept, verbose=FALSE)
{
  t1 <- system.time({
    mdl <- glm.fit(x=x, y=y, family=fam, intercept=intercept)
  })[3]
  #mdl$deviance
  #mdl$null.deviance
  
  
  if (verbose)
    cat("done with R\n")
  
  
  maxiter <- 20
  stoprule <- 3
  
  
  t2 <- system.time({
    mdl2 <- glm_fit(x=x, y=y, family=fam, maxiter, intercept=intercept, stoprule=stoprule)
  })[3]
  
  
  
#  cat(paste("Riter=", mdl$iter, "  New iter=", 
  
  test.coef <- all.equal(mdl$coefficients, mdl2)
  
  cat(paste("Coefficients:", test.coef, "\n"))
#  cat(paste("R time:", t1, "   my time:", t2, "\n"))
  
  
  invisible()
}







maxiter <- 20
  stoprule <- 3
t2 <- system.time({
  mdl2 <- glm_fit(x=x, y=y, family=fam, maxiter, intercept=intercept, stoprule=stoprule)
})[3]
t2




n <- 10 
p <- 3

x <- matrix(rnorm(n*p, mean=10, sd=1), n)
y <- rpois(n, 5)

test_dist(x, y, fam, intercept, verbose=TRUE)


