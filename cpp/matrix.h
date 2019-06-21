/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017-2019 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#ifndef _REDUKTI_MATRIX_H_
#define _REDUKTI_MATRIX_H_

#include <linalg.h>
#include <logger.h>

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

#ifdef _MSC_VER
#include <malloc.h>
#define alloca _alloca
#else
#include <alloca.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Note that the matrix data is assumed to be
// column major
typedef struct redukti_matrix_t redukti_matrix_t;
struct redukti_matrix_t {
	int32_t m; /* rows */
	int32_t n; /* columns */
	double *data;
};

static inline double redukti_matrix_get(const struct redukti_matrix_t *A, int32_t row, int32_t col)
{
	assert(row < A->m && col < A->n);
	int32_t pos = col * A->m + row;
	return A->data[pos];
}

static inline void redukti_matrix_set(struct redukti_matrix_t *A, int32_t row, int32_t col, double v)
{
	assert(row < A->m && col < A->n);
	int32_t pos = col * A->m + row;
	A->data[pos] = v;
}

// As the matrix is column major we can
// efficiently iterate over a column
static inline double *redukti_matrix_column_begin(struct redukti_matrix_t *A, int32_t col)
{
	assert(col < A->n);
	int32_t pos = col * A->m;
	return &A->data[pos];
}
static inline double *redukti_matrix_column_end(struct redukti_matrix_t *A, int32_t col)
{
	assert(col < A->n);
	int32_t pos = col * A->m + (A->m - 1) + 1;
	return &A->data[pos];
}

// Dump contents of the matrix
extern void redukti_matrix_dump_to(redukti_matrix_t *A, const char *desc, FILE *file);

enum redukti_matrix_norm_type { RAVI_MATRIX_NORM_ONE, RAVI_MATRIX_NORM_INFINITY, RAVI_MATRIX_NORM_FROBENIUS };
typedef enum redukti_matrix_norm_type redukti_matrix_norm_type;

// workhorse for matrix multiplication
// C=alpha*A*B + beta*C
// Simple wrapper around dgemm
// To perform C=A*B set alpha to 1.0 and beta to 0.0
static inline void redukti_matrix_multiply(redukti_matrix_t *A, redukti_matrix_t *B, redukti_matrix_t *C,
					   bool transposeA, bool transposeB, double alpha, double beta)
{
	int m = transposeA ? A->n : A->m; /* If transposing then m = columns(A) else rows(A) */
	int n = transposeB ? B->m : B->n; /* If transposing then n = rows(B) else columns(B) */
	int k = transposeA ? A->m : A->n; /* If transposing A then k = rows(A) else columns(A) */
	assert(C->m == m);
	assert(C->n == n);
	int lda = A->m;
	int ldb = B->m;
	int ldc = C->m;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
	cblas_dgemm(CblasColMajor, transposeA ? CblasTrans : CblasNoTrans, transposeB ? CblasTrans : CblasNoTrans, m, n,
		    k, alpha, A->data, lda, B->data, ldb, beta, C->data, ldc);
#else
	char transa = transposeA ? 'T' : 'N';
	char transb = transposeB ? 'T' : 'N';
	dgemm_(&transa, &transb, &m, &n, &k, &alpha, A->data, &lda, B->data, &ldb, &beta, C->data, &ldc);
#endif
}

static inline double redukti_matrix_norm(redukti_matrix_t *A, redukti_matrix_norm_type normType)
{
	char normt = normType == RAVI_MATRIX_NORM_ONE ? '1' : (normType == RAVI_MATRIX_NORM_INFINITY ? 'i' : 'f');
	int lda = A->m;
	double *work = (double *)alloca(sizeof(double) * A->m);
	return dlange_(&normt, &A->m, &A->n, A->data, &lda, work);
}

// LU Factorisation
// The matrix 'A' will be updated
// ipsize must be at least min(m,n)
// ipiv must be an array of size ipsize
// Returns 0 on success
static inline int redukti_matrix_lufactor(redukti_matrix_t *A, int32_t ipsize, int *ipiv)
{
	int lda, info = 0;
	int size = A->m < A->n ? A->m : A->n;
	if (ipsize < size) {
		info = -5;
		goto Lerror;
	}
	for (int i = 0; i < ipsize; i++)
		ipiv[i] = 0;
	lda = (1 < A->m ? A->m : 1);
	dgetrf_(&A->m, &A->n, A->data, &lda, ipiv, &info);
	if (info != 0)
		trace("Error ocurred during LU factorization: %d\n", info);
Lerror:
	return info;
}

// A = A*scalar
static inline void redukti_matrix_scalar_multiply(redukti_matrix_t *A, double scalar)
{
	int n = A->m * A->n;
	int inca = 1;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
	cblas_dscal(n, scalar, A->data, inca);
#else
	dscal_(&n, &scalar, A->data, &inca);
#endif
}

// A = A + alpha*B
static inline void redukti_matrix_add(redukti_matrix_t *A, redukti_matrix_t *B, double alpha)
{
	int n = A->m * A->n;
	int incx = 1;
	int incy = 1;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
	cblas_daxpy(n, alpha, B->data, incx, A->data, incy);
#elif 0
	const double *e = &A->data[A->m * A->n];
	double *p = A->data;
	const double *q = B->data;
	for (; p != e; p++, q++)
		*p += alpha * (*q);
#else
	daxpy_(&n, &alpha, B->data, &incx, A->data, &incy);
#endif
}

// Copy B into A
// A = copy(B)
static inline void redukti_matrix_copy(redukti_matrix_t *A, const redukti_matrix_t *B)
{
	assert(A->m == B->m && A->n == B->n);
	const double *e = &A->data[A->m * A->n];
	double *p = &A->data[0];
	const double *q = &B->data[0];
	for (; p != e; p++, q++)
		*p = *q;
}

extern bool redukti_matrix_estimate_rcond(const redukti_matrix_t *A, double *rcond);

// SVD
// S must be matrix min(m,n) x 1 (vector)
// U must be matrix m x m
// V must be matrix n x n
extern bool redukti_matrix_svd(const redukti_matrix_t *A, redukti_matrix_t *S, redukti_matrix_t *U,
			       redukti_matrix_t *V);

// M must be a square matrix
// V a column vector with rows same as M
extern bool redukti_matrix_solve(redukti_matrix_t *M, redukti_matrix_t *V);

// M must rows > columns
// V a column vector with rows same as M
bool redukti_matrix_lsq_solve(redukti_matrix_t *M, redukti_matrix_t *V, double rcond, bool use_svd);

// A must be square matrix
bool redukti_matrix_inverse(redukti_matrix_t *A);

// transposed must be size nxm where original is sized mxn
// B = transpose(A)
static inline void redukti_matrix_transpose(redukti_matrix_t *B, const redukti_matrix_t *A)
{
#if defined(USE_OPENBLAS)
	double scale = 1.0;
	int lda = A->m > 1 ? A->m : 1; // max(1, rows);
	int ldb = A->n > 1 ? A->n : 1; // max(1, cols);
	cblas_domatcopy(CblasColMajor, CblasTrans, A->m, A->n, scale, A->data, lda, B->data, ldb);
#elif defined(USE_MKL)
	double scale = 1.0;
	int lda = A->m > 1 ? A->m : 1; // max(1, rows);
	int ldb = A->n > 1 ? A->n : 1; // max(1, cols);
	mkl_domatcopy(CblasColMajor, CblasTrans, A->m, A->n, scale, A->data, lda, B->data, ldb);
#else
	/*
	result is column order so we can use
	more optimised iterator
	*/
	double *rp = B->data;
	for (int32_t i = 0; i < A->m; i++) {
		for (int32_t j = 0; j < A->n; j++) {
			*rp++ = A->data[j * A->m + i];
		}
	}
#endif
	B->m = A->n;
	B->n = A->m;
}

static inline void redukti_vector_outer_product(int32_t m, const double *x, int32_t n, const double *y, double *a,
						double alpha)
{
	int incx = 1, incy = 1;
#if defined(__APPLE__) || defined(USE_OPENBLAS)
	cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, a, m);
#else
	dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &m);
#endif
}

#ifdef __cplusplus
}
#endif

#endif
