# No longer needed since switching to F2003 using iso_c_binding, but it's a nice solution
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
