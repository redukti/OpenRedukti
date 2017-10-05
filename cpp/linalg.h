/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#ifndef _REDUKTI_LINALG_H
#define _REDUKTI_LINALG_H

#if !defined(__APPLE__) && !defined(USE_OPENBLAS) && !defined(USE_MKL)

// using reference BLAS and LAPACK

#ifdef __cplusplus
extern "C" {
#endif

extern void ddisna_(char *job, int *m, int *n, double *d, double *sep, int *info);
extern void dgesv_(int *N, int *nrhs, double *A, int *lda, int *ipiv, double *b, int *ldb, int *info);
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b,
		   int *ldb, double *beta, double *c, int *ldc);
extern void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta,
		   double *y, int *incy);
extern void dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *s, double *rcond,
		    int *rank, double *work, int *lwork, int *iwork, int *info);
extern void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *jpvt, double *rcond,
		    int *rank, double *work, int *lwork, int *info);
extern void dgesvx_(char *fact, char *trans, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv,
		    char *equed, double *r__, double *c__, double *b, int *ldb, double *x, int *ldx, double *rcond,
		    double *ferr, double *berr, double *work, int *iwork, int *info);
extern void dgecon_(const char *norm, const int *n, double *a, const int *lda, const double *anorm, double *rcond,
		    double *work, int *iwork, int *info);
extern double dlange_(const char *norm, const int *m, const int *n, const double *a, const int *lda, double *work);
extern void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt,
		    int *ldvt, double *work, int *lwork, int *iwork, int *info);
extern void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
extern void dscal_(int *n, const double *alpha, double *a, int *inca);
extern void daxpy_(int *n, const double *alpha, const double *x, int *incx, double *y, int *incy);
extern void dger_(int *m, int *n, const double *alpha, const double *x, int *incx, const double *y, int *incy,
		  double *A, int *lda);

#ifdef __cplusplus
}
#endif

#elif defined(USE_OPENBLAS)

#include <cblas.h>
#include <f77blas.h>
#include <openblas_config.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *s, double *rcond,
		    int *rank, double *work, int *lwork, int *iwork, int *info);
extern void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *jpvt, double *rcond,
		    int *rank, double *work, int *lwork, int *info);
extern void dgecon_(const char *norm, const int *n, double *a, const int *lda, const double *anorm, double *rcond,
		    double *work, int *iwork, int *info);
extern double dlange_(const char *norm, const int *m, const int *n, const double *a, const int *lda, double *work);
extern void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt,
		    int *ldvt, double *work, int *lwork, int *iwork, int *info);
extern void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

#ifdef __cplusplus
}
#endif

#elif defined(USE_MKL)

#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>

#elif defined(__APPLE__)

// Apple provides cblas and clapack interfaces
#include <Accelerate/Accelerate.h>
#define cblas_daxpby catlas_daxpby

#endif /* __APPLE__ */

#endif