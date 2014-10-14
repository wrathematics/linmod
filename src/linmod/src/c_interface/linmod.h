#ifndef __LINMOD_H_
#define __LINMOD_H_


typedef struct qr_t
{
  int p;
  double *qraux;
  int *pivot;
  double tol;
  int rank;
} qr_t;


// Fortran
void rdgels(int *m, int *n, int *nrhs, double *a, double *b, double *tol,
  double *coef, double *eff, double *ft, double *rsd, double *tau, 
  int *jpvt, int *rank, double *work, int *lwork, int *info);

void lm_fit(int *m, int *n, int *nrhs, double *a, double *b, 
  double *tol, double *coef, double *eff, double *ft, double *rsd,
  double *tau, int *jpvt, int *rank, int *info);

void glm_fit(int *family, int *link, int *intercept, int *stoprule,
  int *n, int *p, double *x, double *y, double *beta, double *wt,
  double *offset, double *resids, int *maxiter, double *tol, int *trace,
  int *info);


// lm_fit.c
int lm_fit_classic(char trans, int m, int n, int nrhs, double *a, double *b);

int lm_fit_alloc(char trans, int m, int n, int nrhs, double *a, double *b,
  double tol, double **eff, double **ft, double **rsd, double **tau,
  int rank);


#endif
