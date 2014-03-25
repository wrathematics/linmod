#include <SEXPtools.h> 

int lmfit_classic(char trans, int m, int n, int nrhs, double *a, double *b);
int lmfit_R(char trans, int m, int n, int nrhs, double *a, double *b,
  double tol, double **coef, double **eff, double **ft, double **rsd, double **tau, 
  int rank);


SEXP R_LM_FIT(SEXP a, SEXP b)
{
  R_INIT;
  
  int m = nrows(a), n = ncols(a);
  int nrhs = ncols(b);
  
  char trans = 'n';
  int info = 0;
  
  double *work, tmpwork;
  int lwork = -1;
  
  SEXP ret, ret_names;
  SEXP a_out, b_out;
  
  
  newRmat(a_out, m, n, "dbl");
  memcpy(DBLP(a_out), DBLP(a), m*n*sizeof(double));
  
  if (nrhs == 1)
  {
    newRvec(b_out, n, "dbl");
  }
  else
  {
    newRmat(b_out, n, nrhs, "dbl");
  }
  
  memcpy(DBLP(b_out), DBLP(b), m*nrhs*sizeof(double));
  
  
  info = lmfit_classic(trans, m, n, nrhs, DBLP(a_out), DBLP(b_out));
  
  
  if (info != 0)
    Rprintf("WARNING : returned info = %d\n", info);
  
  ret_names = make_list_names(2, "qr", "coefficients");
  ret = make_list(ret_names, 2, a_out, b_out);
  
  R_END;
  return ret;
}












SEXP R_LM_FIT_R(SEXP a, SEXP b, SEXP tol)
{
  R_INIT;
  
  int m = nrows(a), n = ncols(a);
  int nrhs = ncols(b);
  
  char trans = 'n';
  int info = 0;
  
  double *work, tmpwork;
  int lwork = -1;
  
  SEXP ret, ret_names;
  SEXP qr, qr_names;
  SEXP rank;
  SEXP a_out, b_out;
  SEXP coef, eff, ft, rsd, tau, qraux;
  
  
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
    newRvec(tau, m, "dbl");
  }
  else
  {
    newRmat(b_out, m, nrhs, "dbl");
    newRmat(coef, n, nrhs, "dbl");
    newRmat(eff, m, nrhs, "dbl");
    newRmat(ft, m, nrhs, "dbl");
    newRmat(rsd, m, nrhs, "dbl");
    newRmat(tau, m, nrhs, "dbl");
  }
  
  memcpy(DBLP(b_out), DBLP(b), m*nrhs*sizeof(double));
  
  
  // Workspace query
  dgels_(&trans, &m, &n, &nrhs, DBLP(a_out), &m, DBLP(b_out), &m, &tmpwork, &lwork, &info);
  
  lwork = (int) tmpwork;
  work = malloc(lwork * sizeof(*work));
  
  
  // Fit y~x
  rdgels_(&trans, &m, &n, &nrhs, DBLP(a_out), &m, DBLP(b_out), &m, work, &lwork, &info, &tol, DBLP(coef), DBLP(eff), DBLP(ft), DBLP(rsd), DBLP(tau), INTP(rank));
  
  
  if (info != 0)
    Rprintf("WARNING : returned info = %d\n", info);
  
  // FIXME qraux == work(1:rank)
  
  // Manage return
  qr_names = make_list_names(5, "qr", "qraux", "pivot", "tol", "rank");
  qr = make_list(qr_names, 5, a_out, tau, RNULL, tol, rank);
  
  ret_names = make_list_names(7, "coefficients", "residuals", "effects", "rank", "fitted.values", "assign", "qr");
  ret = make_list(ret_names, 7, coef, rsd, eff, rank, ft, RNULL, qr);
  
  R_END;
  return ret;
}

