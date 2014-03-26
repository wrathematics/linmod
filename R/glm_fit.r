string_c2f <- function(str, len)
{
  tmp <- len - length(unlist(strsplit(str, split="")))
  
  if (tmp < 0)
    stop("'str' contains more than 'len' chars")
  else if (tmp == 0L)
    return(str)
  else
    return( paste(str, paste(rep(" ", tmp), collapse="", sep=""), collapse="", sep="") )
}


#fit_logistic <- function(x, y, maxiter=50, tol=10*.Machine$double.eps, intercept=TRUE, stoprule='1')
glm_fit <- function(x, y, family, maxiter=50, tol=1e-8, offset=rep(0.0, nobs), intercept=TRUE, stoprule=3)
{
  n <- dim(x)[1L]
  p <- dim(x)[2L]
  
  nobs <- NROW(y)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  if (!is.double(offset))
    storage.mode(offset) <- "double"
  
  if (intercept)
    incpt <- 'Y'
  else
    incpt <- 'N'
  
  fam <- string_c2f(family$family, 8L)
  link <- string_c2f(family$link, 8L)
  
  fam <- tolower(fam)
  link <- tolower(link)
  
  fit <- .Call("R_GLM_FIT", 
               as.character(fam), as.character(link),
               as.character(incpt), as.integer(stoprule), 
               as.integer(n), as.integer(p),
               x, y, offset, as.integer(maxiter), as.double(tol),
               PACKAGE = "linmod")
  
  return(fit)
}
