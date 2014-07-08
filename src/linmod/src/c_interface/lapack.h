#ifndef __LINMOD_LAPACK__
#define __LINMOD_LAPACK__


void dgels_(char* trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);


#endif
