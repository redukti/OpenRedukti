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

#include <logger.h>
#include <matrix.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static redukti_matrix_t *matrix_dup(const redukti_matrix_t *src)
{
	size_t msize = sizeof(double) * src->m * src->n;
	redukti_matrix_t *A = (redukti_matrix_t *)calloc(1, sizeof(redukti_matrix_t));
	double *data = (double *)calloc(1, msize);
	A->m = src->m;
	A->n = src->n;
	memcpy(data, src->data, msize);
	A->data = data;
	return A;
}

static inline void free_matrix(redukti_matrix_t *A)
{
	free(A->data);
	free(A);
}

// SVD
// S must be matrix min(m,n) x 1 (vector)
// U must be matrix m x m
// V must be matrix n x n
bool redukti_matrix_svd(const redukti_matrix_t *input, redukti_matrix_t *S, redukti_matrix_t *U, redukti_matrix_t *V)
{
	double workSize, *work = NULL, *a = NULL;
	int ldu, ldvt, lwork, info = -1, *iwork = NULL;
	char job = 'A';
	redukti_matrix_t *A = NULL;
	int m = input->m, n = input->n, lda = m;
	int ss = m < n ? m : n;
	int s_size = 1 < ss ? ss : 1;
	if (S->m != s_size || S->n != 1) {
		error("The vector S must be a column vector of size %d\n", s_size);
		goto done;
	}
	if (U->m != m || U->n != m) {
		error("The matrix U must be of size %dx%d\n", m, m);
		goto done;
	}
	ldu = U->m;
	if (V->m != n || V->n != n) {
		error("The matrix V must be of size %dx%d\n", n, n);
		goto done;
	}
	A = matrix_dup(input);
	if (!A) {
		error("Failed to allocate memory\n");
		goto done;
	}
	a = &A->data[0];
	ldvt = n;
	workSize = 0;
	work = &workSize;
	lwork = -1; // we want to estimate workSize first
	info = 0;
	int min_sz = 8 * (m < n ? m : n);    // min(m,n)*8
	int works = 1 < min_sz ? min_sz : 1; // max
	iwork = (int *)calloc(works, sizeof(int));
	if (iwork == NULL) {
		error("Failed to allocate memory\n");
		goto done;
	}
	dgesdd_(&job, &m, &n, a, &lda, &S->data[0], &U->data[0], &ldu, &V->data[0], &ldvt, work, &lwork, iwork, &info);
	if (info == 0) {
		lwork = (int)workSize;
		work = (double *)calloc(lwork, sizeof(double));
		if (work) {
			dgesdd_(&job, &m, &n, a, &lda, &S->data[0], &U->data[0], &ldu, &V->data[0], &ldvt, work, &lwork,
				iwork, &info);
		} else {
			error("Failed to allocate memory\n");
			info = -1;
		}
	} else {
		error("Failed to estimate work size for SVD: info=%d\n", info);
	}
	if (info != 0) {
		error("Failed to compute SVD: info=%d\n", info);
	}
done:
	if (A)
		free_matrix(A);
	if (iwork)
		free(iwork);
	if (work && work != &workSize)
		free(work);
	return info == 0;
}

bool redukti_matrix_estimate_rcond(const redukti_matrix_t *A, double *rcond)
{
	bool ok = false;
	redukti_matrix_t *copy_of_A = matrix_dup(A);
	double anorm = redukti_matrix_norm(copy_of_A, RAVI_MATRIX_NORM_ONE);
	int ipsize = A->m < A->n ? A->m : A->n;
	int *ipiv = (int *)alloca(sizeof(int) * ipsize);
	int info = redukti_matrix_lufactor(copy_of_A, ipsize, ipiv);
	if (info != 0) {
		error("failed to estimate rcond (LU factor failed)\n");
		goto done;
	}
	char norm[] = {'1'};
	int n = A->n;
	double *a = &copy_of_A->data[0];
	int lda = A->m;
	if (lda < n) {
		error("failed to estimate rcond (LDA < n)\n");
		goto done;
	}
	double *work = (double *)alloca(sizeof(double) * n * 4);
	int *iwork = (int *)alloca(sizeof(int) * n);
	dgecon_(&norm[0], &n, a, &lda, &anorm, rcond, work, iwork, &info);
	if (info != 0) {
		error("failed to estimate rcond (DGECON failed)\n");
		goto done;
	}
	ok = true;
done:
	if (copy_of_A)
		free_matrix(copy_of_A);
	return ok;
}

bool redukti_matrix_solve(redukti_matrix_t *m, redukti_matrix_t *v)
{
	if (m->m == 0 || m->n == 0 || v->m == 0 || v->n != 1) {
		error("The matrix A and vector y must have rows > 0\n");
		assert(false);
		return false;
	}
	if (m->m != m->n || m->m != v->m) {
		error("The default solver only accepts n x n matrix\n");
		assert(false);
		return false;
	}

	int N = m->n;				    // order of the matrix, number of rows of vector v
	int nrhs = 1;				    // number of columns of vector v
	int lda = m->m;				    // leading dimension of the array for A (rows)
	int *ipiv = (int *)alloca(N * sizeof(int)); // return value: containing pivot indices
	memset(ipiv, 0, N * sizeof(int));

	int ldb = v->m; // leading dimension of vector v
	int info = 0;   // return value: if 0 successful
	double *A = &m->data[0];
	double *b = &v->data[0];

	dgesv_(&N, &nrhs, A, &lda, ipiv, b, &ldb, &info);

	return info == 0;
}

bool redukti_matrix_lsq_solve(redukti_matrix_t *A, redukti_matrix_t *y, double rcond, bool svd)
{
	if (A->m == 0 || A->n == 0 || y->m == 0 || y->n != 1) {
		error("The matrix A and vector y must have rows > 0\n");
		assert(false);
		return false;
	}
	if (A->m < A->n) {
		error("The matrix A must have rows >= cols\n");
		assert(false);
		return false;
	}
	if (y->m != A->m) {
		error("The vector y must have rows = A.rows\n");
		assert(false);
		return false;
	}
	if (rcond <= 0) {
		if (!redukti_matrix_estimate_rcond(A, &rcond))
			return false;
	}
	int m = A->m;
	int n = A->n;
	int nrhs = 1; // number of columns in b
	double *a = &A->data[0];
	int lda = m;
	double *b = &y->data[0];
	int ldb = m;
	int *jpvt = (int *)alloca(sizeof(int) * n);
	memset(jpvt, 0, sizeof(int) * n);
	double *s = (double *)alloca(sizeof(double) * (m > n ? m : n));
	int lwork = -1;
	int liwork = -1;
	int rank = 0;
	double temp = 0;
	int info = 0;
	// First estimate the work space required
	if (svd) {
		dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, &temp, &lwork, &liwork, &info);
	} else {
		dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, &temp, &lwork, &info);
	}
	if (info != 0) {
		error("failed to estimate work space requirement\n");
		return false;
	}
	// allocate work space
	lwork = (int)temp;
	assert(lwork > 0);
	double *ptr = (double *)calloc(lwork, sizeof(double));
	assert(ptr);
	// solve
	if (svd) {
		if (liwork == -1) {
			// Looks like on MAC OSX the dgelsd call above does not
			// set liwork.
			double minmn = m < n ? m : n;  // min
			minmn = 1 < minmn ? minmn : 1; // max
			double smlsiz_plus_1 = 26.;
			double nlvl = (int)(log(minmn / smlsiz_plus_1) / log(2.0)) + 1;
			nlvl = 0.0 < nlvl ? nlvl : 0.0; // max
			liwork = (int)(3.0 * minmn * nlvl + 11.0 * minmn);
			// Double it for safety
			liwork *= 2;
		}
		int *iwork_ptr = (int *)calloc(liwork, sizeof(int));
		assert(iwork_ptr);
		dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, ptr, &lwork, iwork_ptr, &info);
		free(iwork_ptr);
	} else {
		dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, ptr, &lwork, &info);
	}
	free(ptr);
	if (info != 0) {
		error("failed to solve linear least squares problem\n");
		return false;
	}
	return true;
}

bool redukti_matrix_inverse(redukti_matrix_t *a)
{
	int m = a->m;
	int n = a->n;
	if (m != n) {
		error("matrix is not square");
		return false;
	}
	int info = 0;
	int size = m < n ? m : n; // min
	int *ipiv = (int *)alloca(sizeof(int) * size);
	memset(ipiv, 0, sizeof(int) * size);
	int lda = (1 < m ? m : 1);
	dgetrf_(&m, &n, &a->data[0], &lda, ipiv, &info);
	if (info != 0) {
		error("failed LU factorization of input matrix");
		return false;
	}
	n = a->n;
	lda = n;
	int lwork = n * n;
	double *work = (double *)alloca(sizeof(double) * lwork);
	info = 0;
	dgetri_(&n, &a->data[0], &lda, ipiv, work, &lwork, &info);
	if (info != 0) {
		error("failed to compute inverse of matrix");
		return false;
	}
	return true;
}

void redukti_matrix_dump_to(redukti_matrix_t *A, const char *desc, FILE *file)
{
	fprintf(file, "%s\n", desc);
	for (int32_t i = 0; i < A->m; i++) {
		for (int32_t j = 0; j < A->n; j++) {
			fprintf(file, "%.12f\t", redukti_matrix_get(A, i, j));
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}
