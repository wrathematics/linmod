#' Linear and Generalized Linear Moddels
#' 
#' The linmod package contains wrappers of the linmod library, a collection of
#' utilities for fitting linear and generalized linear models.  These use
#' modified versions of LAPACK routines (compared to R's use of LINPACK
#' routines).
#' 
#' \tabular{ll}{ 
#'    Package: \tab linmod \cr 
#'    Type: \tab Package \cr 
#'    License: \tab MPL 2.0\cr
#' }
#' 
#' @import RNACI
#' @useDynLib linmod, R_LM_FIT, R_GLM_FIT
#' 
#' @name linmod-package
#' @docType package
#' @author Drew Schmidt \email{wrathematics AT gmail.com}
#' @keywords Package
NULL

