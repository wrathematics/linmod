/*#include <SEXPtools.h> */

#include <R.h>
#include <Rinternals.h>

// character pointers
#define CHARPT(x,i) ((char*)CHAR(STRING_ELT(x,i)))


#define PTCT(X) PROTECT(X);protect_ct++;


SEXP R_GLM_FIT(SEXP FAMILY, SEXP LINK, SEXP INCPT, SEXP STOPRULE,
  SEXP N, SEXP P, SEXP X, SEXP Y, SEXP OFFSET, SEXP MAXITER, SEXP TOL)
{
  int protect_ct = 0;
  
  const int p = INTEGER(P)[0], n = INTEGER(N)[0];
  int info = 0;
  
  // ----------------------------------------------
  // Wrangling the return and SEXP allocation
  // ----------------------------------------------
  
  // QR
  SEXP QR, QR_NAMES;
  PTCT(QR = allocVector(VECSXP, 5));
  PTCT(QR_NAMES = allocVector(STRSXP, 5));
  
  SEXP CP_X, RANK, QRAUX, PIVOT;
  PTCT(CP_X = allocMatrix(REALSXP, n, p));
  PTCT(RANK = allocVector(VECSXP, 1));
  
  SET_VECTOR_ELT(QR, 0, CP_X);
  SET_VECTOR_ELT(QR, 1, RANK);
  SET_VECTOR_ELT(QR, 2, QRAUX);
  SET_VECTOR_ELT(QR, 3, PIVOT);
  SET_VECTOR_ELT(QR, 4, TOL);
  
  SET_STRING_ELT(QR_NAMES, 0, mkChar("qr")); 
  SET_STRING_ELT(QR_NAMES, 1, mkChar("rank")); 
  SET_STRING_ELT(QR_NAMES, 2, mkChar("qraux")); 
  SET_STRING_ELT(QR_NAMES, 3, mkChar("pivot")); 
  SET_STRING_ELT(QR_NAMES, 4, mkChar("tol")); 
  
  setAttrib(QR, R_NamesSymbol, QR_NAMES);
  
  // Return
  SEXP RET, RET_NAMES;
  PTCT(RET = allocVector(VECSXP, 20));
  PTCT(RET_NAMES = allocVector(STRSXP, 20));
  
  SEXP BETA, RESIDS, FITTED, EFFECTS, R, LINPRED, DEV, AIC, NULLDEV, 
       ITER, WT, WT_OLD, DFRESID, DFNULL, CONV, BDD;
  PTCT(BETA = allocVector(REALSXP, p));
  PTCT(RESIDS = allocVector(REALSXP, n));
  PTCT(FITTED = allocVector(REALSXP, n));
  PTCT(EFFECTS = allocVector(REALSXP, n));
  PTCT(R = allocMatrix(REALSXP, p, p));
  PTCT(LINPRED = allocVector(REALSXP, n));
  PTCT(DEV = allocVector(REALSXP, 1));
  PTCT(AIC = allocVector(REALSXP, 1));
  PTCT(NULLDEV = allocVector(REALSXP, 1));
  PTCT(ITER = allocVector(INTSXP, 1));
  PTCT(WT = allocVector(REALSXP, n));
  PTCT(WT_OLD = allocVector(REALSXP, n));
  PTCT(DFRESID = allocVector(INTSXP, 1));
  PTCT(DFNULL = allocVector(INTSXP, 1));
  PTCT(CONV = allocVector(INTSXP, 1));
  PTCT(BDD = allocVector(INTSXP, 1));
  
  SET_VECTOR_ELT(RET, 0, BETA);
  SET_VECTOR_ELT(RET, 1, RESIDS);
  SET_VECTOR_ELT(RET, 2, FITTED);
  SET_VECTOR_ELT(RET, 3, EFFECTS);
  SET_VECTOR_ELT(RET, 4, R);
  SET_VECTOR_ELT(RET, 5, RANK);
  SET_VECTOR_ELT(RET, 6, QR);
//  SET_VECTOR_ELT(RET, ); 
  SET_VECTOR_ELT(RET, 8, LINPRED);
  SET_VECTOR_ELT(RET, 9, DEV);
  SET_VECTOR_ELT(RET, 10, AIC);
  SET_VECTOR_ELT(RET, 11, NULLDEV);
  SET_VECTOR_ELT(RET, 12, ITER);
  SET_VECTOR_ELT(RET, 13, WT);
  SET_VECTOR_ELT(RET, 14, WT_OLD);
  SET_VECTOR_ELT(RET, 15, DFRESID);
  SET_VECTOR_ELT(RET, 16, DFNULL);
  SET_VECTOR_ELT(RET, 17, Y);
  SET_VECTOR_ELT(RET, 18, CONV);
  SET_VECTOR_ELT(RET, 19, BDD);
  
  SET_STRING_ELT(RET_NAMES, 0, mkChar("coefficients")); 
  SET_STRING_ELT(RET_NAMES, 1, mkChar("residuals")); 
  SET_STRING_ELT(RET_NAMES, 2, mkChar("fitted.values")); 
  SET_STRING_ELT(RET_NAMES, 3, mkChar("effects")); 
  SET_STRING_ELT(RET_NAMES, 4, mkChar("R")); 
  SET_STRING_ELT(RET_NAMES, 5, mkChar("rank")); 
  SET_STRING_ELT(RET_NAMES, 6, mkChar("qr")); 
  SET_STRING_ELT(RET_NAMES, 7, mkChar("family")); 
  SET_STRING_ELT(RET_NAMES, 8, mkChar("linear.predictors")); 
  SET_STRING_ELT(RET_NAMES, 9, mkChar("deviance")); 
  SET_STRING_ELT(RET_NAMES, 10, mkChar("aic")); 
  SET_STRING_ELT(RET_NAMES, 11, mkChar("null.deviance")); 
  SET_STRING_ELT(RET_NAMES, 12, mkChar("iter")); 
  SET_STRING_ELT(RET_NAMES, 13, mkChar("weights")); 
  SET_STRING_ELT(RET_NAMES, 14, mkChar("prior.weights")); 
  SET_STRING_ELT(RET_NAMES, 15, mkChar("df.residual")); 
  SET_STRING_ELT(RET_NAMES, 16, mkChar("df.null")); 
  SET_STRING_ELT(RET_NAMES, 17, mkChar("y")); 
  SET_STRING_ELT(RET_NAMES, 18, mkChar("converged")); 
  SET_STRING_ELT(RET_NAMES, 19, mkChar("boundary")); 
  
  setAttrib(RET, R_NamesSymbol, RET_NAMES);
  
  // ----------------------------------------------
  // Actual work
  // ----------------------------------------------
  
  // copy X
  memcpy(REAL(CP_X), REAL(X), n*p * sizeof(double));
  
  
  // call fortran
  glm_fit_(CHARPT(FAMILY, 0), CHARPT(LINK, 0), 
    CHARPT(INCPT, 0), INTEGER(STOPRULE), &n, &p, REAL(CP_X), 
    REAL(Y), REAL(BETA), REAL(WT), REAL(OFFSET), REAL(RESIDS), 
    INTEGER(MAXITER), REAL(TOL), &info);
  
  if (info != 0)
    Rprintf("WARNING : returned info = %d\n", info);
  
  UNPROTECT(protect_ct);
  return(BETA);
} 

