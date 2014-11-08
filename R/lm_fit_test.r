lmtest <- function()
{
  set.seed(1234)
  x <- matrix(rnorm(30), 10, 3)
  y <- rnorm(10)
  
  lm_fit(x, y)
}



lmfit_test_checker <- function(mdl1, mdl2, which, whichprint)
{
  a <- eval(parse(text=paste("mdl1$", which, sep="")))
  b <- eval(parse(text=paste("mdl2$", which, sep="")))
  test <- all.equal(a, b, check.names=FALSE)
  
  if (which %in% whichprint)
  {
    cat(paste("\n", which, ":\n", sep=""))
    print(a)
    print(b)
  }
  
  return( test )
}



lmfit_test <- function(x, y, offset, check.rank=FALSE, whichprint=NULL, verbose=FALSE)
{
  if (verbose)
    cat("-------------------------------------------------------\n\n")
  
  mdl1 <- stats::lm.fit(x=x, y=y, offset=offset)
  mdl2 <- linmod::lm_fit(x=x, y=y, offset=offset, check.rank=check.rank)
  
  ### QR storage will almost certainly be different because of the
  ### way the algorithms differ
  mdl1$qr$qr <- mdl2$qr$qr <- mdl1$qr$qraux <- mdl2$qr$qraux <- NULL
  
  checks <- c("coefficients", "residuals", "effects", "fitted", "rank", "fitted.values", "assign", "qr")
  res <- list()
  
  for (i in 1:length(checks))
    res[i] <- lmfit_test_checker(mdl1=mdl1, mdl2=mdl2, which=checks[i], whichprint=whichprint)
  
  names(res) <- checks
  
  if (all(sapply(res, is.logical)))
  {
    if (all(as.logical(res)) && verbose)
      cat("\nAll checks passed!\n\n")
    
    return(TRUE)
  }
  else
  {
    if (verbose)
    {
      cat("FAILURE:\n\n")
      print(res)
    }
    
    print(mdl1$qr$pivot)
    print(mdl2$qr$pivot)
    
    return(FALSE)
  }
}

