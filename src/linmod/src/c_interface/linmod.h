#ifndef __LINMOD__
#define __LINMOD__


typedef struct qr_t
{
  int p;
  double *qraux;
  int *pivot;
  double tol;
  int rank;
} qr_t;


// Fortran
void lm_fit_(int *m, int *n, int *nrhs, double *a, double *b, 
  double *tol, double *coef, double *eff, double *ft, double *rsd,
  double *tau, int *jpvt, int *rank, int *info);

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
