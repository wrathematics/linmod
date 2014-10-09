#include "linmod.h"


#define setDimNames(X, Y, Z) \
  newRlist(X, 2); \
  SET_VECTOR_ELT(X, 0, Y); \
  SET_VECTOR_ELT(X, 1, RNULL); \
  setAttrib(Z, R_DimNamesSymbol, X);

SEXP R_LM_FIT(SEXP a, SEXP b, SEXP tol, SEXP checkrank)
{
  R_INIT;
  
  int m = nrows(a), n = ncols(a);
  int nrhs = ncols(b);
  
  int mn = (m<n ? m : n);
  
  char trans = 'n';
  int info = 0;
  
  SEXP ret, ret_names;
  SEXP qr, qr_names;
  SEXP rank, df_residual;
  SEXP a_out, b_out;
  SEXP coef, coef_names;
  SEXP eff, eff_names;
  SEXP ft, rsd, tau, jpvt, qraux;
  SEXP dimnames, effdimnames;
  
  newRvec(rank, 1, "int");
  newRvec(df_residual, 1, "int");
  
  newRmat(a_out, m, n, "dbl");
  memcpy(DBLP(a_out), DBLP(a), m*n*sizeof(double));
  
  if (nrhs == 1)
  {
    newRvec(b_out, m, "dbl");
    newRvec(coef, n, "dbl");
    newRvec(eff, m, "dbl");
    newRvec(ft, m, "dbl");
    newRvec(rsd, m, "dbl");
  }
  else
  {
    newRmat(b_out, m, nrhs, "dbl");
    newRmat(coef, n, nrhs, "dbl");
    newRmat(eff, m, nrhs, "dbl");
    newRmat(ft, m, nrhs, "dbl");
    newRmat(rsd, m, nrhs, "dbl");
  }
  
  newRvec(tau, mn, "dbl");
  newRvec(jpvt, n, "int");
  
  memcpy(DBLP(b_out), DBLP(b), m*nrhs*sizeof(double));
  
  
  if (INT(checkrank) == 0)
    INT(rank) = -1;
  else
    INT(rank) = 0;
  
  
  // Fit y~x
  lm_fit_(&m, &n, &nrhs, DBLP(a_out), DBLP(b_out), DBLP(tol), DBLP(coef), DBLP(eff), DBLP(ft), DBLP(rsd), DBLP(tau), INTP(jpvt), INTP(rank), &info);
  
  
  if (info != 0)
    Rprintf("WARNING : returned info = %d\n", info);
  
  INT(df_residual) = m - INT(rank);
  
  // FIXME qraux == work(1:rank)
  
  // Manage return
  qr_names = make_list_names(5, "qr", "qraux", "pivot", "tol", "rank");
  qr = make_list(qr_names, 5, a_out, tau, jpvt, tol, rank);
  
  coef_names = make_lmfit_default_rownames(n);
  eff_names = make_lmfit_default_effectnames(m, n, INTP(jpvt));
  
  if (nrhs == 1)
  {
    setAttrib(coef, R_NamesSymbol, coef_names);
    setAttrib(eff, R_NamesSymbol, eff_names);
  }
  else
  {
    setDimNames(dimnames, coef_names, coef);
    setDimNames(effdimnames, eff_names, eff);
  }
  
  ret_names = make_list_names(8, "coefficients", "residuals", "effects", "rank", "fitted.values", "assign", "qr", "df.residual");
  ret = make_list(ret_names, 8, coef, rsd, eff, rank, ft, RNULL, qr, df_residual);
  
  R_END;
  return ret;
}

