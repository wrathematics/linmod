#include "linmod.h"
#include "linmod/src/c_interface/glm_defines.h"


SEXP R_GLM_FIT(SEXP FAMILY, SEXP LINK, SEXP INTERCEPT, SEXP STOPRULE, SEXP TRACE,
  SEXP X, SEXP Y, SEXP OFFSET, SEXP MAXITER, SEXP TOL)
{
  R_INIT;
  int n = nrows(X), p = ncols(X);
  int info = 0;
  
  SEXP QR, QR_NAMES;
  SEXP CP_X, RANK, QRAUX, PIVOT;
  
  SEXP RET, RET_NAMES;
  SEXP BETA, RESIDS, FITTED, EFFECTS, R, LINPRED, DEV, AIC, NULLDEV, 
       ITER, WT, WT_OLD, DFRESID, DFNULL, CONV, BDD;
  
  newRmat(CP_X, n, p, "dbl");
  newRlist(RANK, 1); // list ?
  
  newRvec(BETA, p, "dbl");
  newRvec(RESIDS, n, "dbl");
  newRvec(FITTED, n, "dbl");
  newRvec(EFFECTS, n, "dbl");
  newRmat(R, p, p, "dbl");
  newRvec(LINPRED, n, "dbl");
  newRvec(DEV, 1, "dbl");
  newRvec(AIC, 1, "dbl");
  newRvec(NULLDEV, 1, "dbl");
  newRvec(ITER, 1, "int");
  newRvec(WT, n, "dbl");
  newRvec(WT_OLD, n, "dbl");
  newRvec(DFRESID, 1, "dbl");
  newRvec(DFNULL, 1, "dbl");
  newRvec(CONV, 1, "dbl");
  newRvec(BDD, 1, "dbl");
  
  
  memcpy(REAL(CP_X), REAL(X), n*p * sizeof(double));
  
  glm_fit(INTP(FAMILY), INTP(LINK), 
    INTP(INTERCEPT), INTP(STOPRULE), &n, &p, REAL(CP_X), 
    REAL(Y), REAL(BETA), REAL(WT), REAL(OFFSET), REAL(RESIDS), 
    INTP(MAXITER), REAL(TOL), INTP(TRACE), &info);
  
  
  // wrangle return
  QR_NAMES = make_list_names(5, "qr", "rank", "qraux", "pivot", "tol");
  QR = make_list(QR_NAMES, 5, CP_X, RANK, QRAUX, PIVOT, TOL);
  
  RET_NAMES = make_list_names(20, 
    "coefficients", "residuals", "fitted.values", "effects", "R", "rank", 
    "qr", "family", "linear.predictors", "deviance", "aic", "null.deviance",
    "iter", "weights", "prior.weights", "df.residual", "df.null", "y", 
    "converged", "boundary");
  
  RET = make_list(RET_NAMES, 20, 
    BETA, RESIDS, FITTED, EFFECTS, R, RANK, 
    QR, LINPRED, DEV, AIC, NULLDEV,
    ITER, WT, WT_OLD, DFRESID, DFNULL, Y,
    CONV, BDD);
  
  
  if (info != 0)
    Rprintf("WARNING : returned info = %d\n", info);
  
  R_END;
  return BETA;
} 

