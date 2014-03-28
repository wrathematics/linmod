#include <SEXPtools.h> 


SEXP R_LM_FIT(SEXP a, SEXP b, SEXP tol, SEXP checkrank)
{
  R_INIT;
  
  int m = nrows(a), n = ncols(a);
  int nrhs = ncols(b);
  
  int mn = (m<n ? m : n);
  
  char trans = 'n';
  int info = 0;
  
  double *work, tmpwork;
  int lwork = -1;
  
  SEXP ret, ret_names;
  SEXP qr, qr_names;
  SEXP rank;
  SEXP a_out, b_out;
  SEXP coef, coef_names;
  SEXP eff, ft, rsd, tau, jpvt, qraux;
  
  
  newRvec(rank, 1, "int");
  
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
  
  
  if (INT(checkrank, 0) == 0)
    INT(rank, 0) = -1;
  else
    INT(rank, 0) = 0;
  
  
  // Workspace query
  dgels_(&trans, &m, &n, &nrhs, DBLP(a_out), &m, DBLP(b_out), &m, &tmpwork, &lwork, &info);
  
  lwork = (int) tmpwork;
  work = malloc(lwork * sizeof(*work));
  
  
  // Fit y~x
  rdgels_(&trans, &m, &n, &nrhs, DBLP(a_out), &m, DBLP(b_out), &m, work, &lwork, &info, &tol, DBLP(coef), DBLP(eff), DBLP(ft), DBLP(rsd), DBLP(tau), INTP(jpvt), INTP(rank));
  
  
  if (info != 0)
    Rprintf("WARNING : returned info = %d\n", info);
  
  // FIXME qraux == work(1:rank)
  
  // Manage return
  qr_names = make_list_names(5, "qr", "qraux", "pivot", "tol", "rank");
  qr = make_list(qr_names, 5, a_out, tau, jpvt, tol, rank);
  
  coef_names = make_dataframe_default_colnames(n);
  setAttrib(coef, R_NamesSymbol, coef_names);
  
  ret_names = make_list_names(7, "coefficients", "residuals", "effects", "rank", "fitted.values", "assign", "qr");
  ret = make_list(ret_names, 7, coef, rsd, eff, rank, ft, RNULL, qr);
  
  R_END;
  return ret;
}

