glm_test <- function(n, p, family, intercept=TRUE, verbose=FALSE, timings=FALSE)
{
  x <- rnorm(n*p, mean=10, sd=10000)
  dim(x) <- c(n, p)
  
  if (family$family == "binomial")
    y <- sample(as.double(0:1), size=n, replace=TRUE)
  else if (family$family == "poisson")
    y <- rpois(n, 5)
  else if (family$family == "gaussian")
    y <- rnorm(n)
  else if (family$family == "Gamma")
    y <- rgamma(n, 2)
  
  
  if (verbose)
    cat(paste("##################################\nNow testing: ", family$family, "(", family$link, ")\n", sep=""))
  
  
  t1 <- system.time({
    mdl <- glm.fit(x=x, y=y, family=family, intercept=intercept)
  })[3]
  #mdl$deviance
  #mdl$null.deviance
  
  
  if (verbose && timings)
    cat("Done with R\n")
  
  
  maxiter <- 20
  stoprule <- "deviance"
  
  
  t2 <- system.time({
    mdl2 <- glm_fit(x=x, y=y, family=family, maxiter, intercept=intercept, stoprule=stoprule)
  })[3]
  
  
  
#  cat(paste("Riter=", mdl$iter, "  New iter=", 
  
  test.coef <- all.equal(mdl$coefficients, mdl2)
  
  cat(paste("Coefficients:", test.coef, "\n"))
  
  if (timings)
    cat(paste("R time:", t1, "   my time:", t2, "\n"))
  
  
  invisible()
}

