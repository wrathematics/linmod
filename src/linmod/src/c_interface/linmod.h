#ifndef __LINMOD__
#define __LINMOD__


// Fortran
void lm_fit_(int *m, int *n, int *nrhs, double *a, int *lda, 
  double *b, int *ldb, double *tol, double *coef, 
  double *eff, double *ft, double *rsd, double *tau, int *jpvt, 
  int *rank, int *info);

void glm_fit_(int *family, int *link, int *intercept, int *stoprule,
  int *n, int *p, double *x, double *y, double *beta, double *wt,
  double *offset, double *resids, int *maxiter, double *tol, int *trace,
  int *info);


// lm_fit.c
int lmfit_classic(char trans, int m, int n, int nrhs, double *a, double *b);

int lmfit_R(char trans, int m, int n, int nrhs, double *a, double *b,
  double tol, double **eff, double **ft, double **rsd, double **tau,
  int rank);


#endif
