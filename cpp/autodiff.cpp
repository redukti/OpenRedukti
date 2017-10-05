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

/**
 * Parts of the implementation of adouble type is based on:
 * Implementation of a vector-mode version of hyper-dual numbers
 *
 * Written by: Jeffrey A. Fike
 * Stanford University, Department of Aeronautics and Astronautics
 *
 * Copyright (c) 2006 Jeffrey A. Fike
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <autodiff.h>

#include <allocators.h>
#include <linalg.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Compute the memory size of autodiff variable
size_t redukti_adouble_alloc_size(int vars, int order)
{
	static_assert(sizeof(redukti_adouble_t) == sizeof(double) * 2,
		      "redukti_adouble_t size is not equal to 2*sizeof(double)");
	return sizeof(redukti_adouble_t) + (order ? vars * sizeof(double) : 0) +
	       (order == 2 ? vars * vars * sizeof(double) : 0);
}

// Initialize - it is assumed that the autodiff structure has been
// sized correcly
void redukti_adouble_init(redukti_adouble_t *x, int n_vars, int order, int var, double v)
{
	assert(var < n_vars);
	assert(order >= 0);
	assert(order <= 2);

	size_t mysize = redukti_adouble_alloc_size(n_vars, order);
	memset(x, 0, mysize);
	x->order_ = order; // order
	x->vars_ = n_vars; // number of variables
	x->data_[0] = v;   // value
	if (order && var >= 0 && n_vars > 0)
		x->data_[1 + var] = 1.0;
}

void redukti_adouble_assign(redukti_adouble_t *A, const redukti_adouble_t *B)
{
	size_t mysize = redukti_adouble_alloc_size(B->vars_, B->order_);
	memcpy(A, B, mysize);
}

// Compute the number of array elements in autodiff variable
static inline int redukti_adouble_arraysize(int vars, int order)
{
	return 1 + (order ? vars : 0) + (order == 2 ? vars * vars : 0);
}

// A = A + alpha*B
void redukti_adouble_add(redukti_adouble_t *A, redukti_adouble_t *B, double alpha)
{
	assert(A->order_ == B->order_ && A->vars_ == B->vars_);
	if (A->vars_ == 0 || A->order_ == 0) {
		A->data_[0] += alpha * B->data_[0];
	} else {
		int n = redukti_adouble_arraysize(A->vars_, A->order_);
		int incx = 1;
		int incy = 1;
		cblas_daxpy(n, alpha, B->data_, incx, A->data_, incy);
	}
}

// A = A*alpha
void redukti_adouble_scalar_multiply(redukti_adouble_t *A, double alpha)
{
	if (A->vars_ == 0 || A->order_ == 0) {
		A->data_[0] *= alpha;
	} else {
		int n = redukti_adouble_arraysize(A->vars_, A->order_);
		int inca = 1;
		cblas_dscal(n, alpha, A->data_, inca);
	}
}

// A = abs(A)
void redukti_adouble_abs(redukti_adouble_t *A)
{
	if (A->data_[0] < 0.0) {
		redukti_adouble_scalar_multiply(A, -1.0);
	}
}

// A = A*B
void redukti_adouble_multiply(redukti_adouble_t *A, redukti_adouble_t *B, redukti_adouble_t *buffer)
{
	assert(A->order_ == B->order_ && A->vars_ == B->vars_);
	int m = A->vars_, n = A->vars_;
	int mn = m * n;
	double A_scalar = A->data_[0];
	double B_scalar = B->data_[0];
	switch (A->order_) {
	case 2: {
		double *buf = &buffer->data_[0];
		double *A_grad = &A->data_[1];
		double *B_grad = &B->data_[1];
		double *A_hessian = &A->data_[1 + m];
		double *B_hessian = &B->data_[1 + m];
#if 0
		// 
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				auto pos = i * m + j;
				buf[pos] = A_scalar * B_hessian[pos] + B_scalar * A_hessian[pos] +
					   A_grad[i] * B_grad[j] + A_grad[j] * B_grad[i];
			}
		}
		memcpy(A_hessian, buf, mn * sizeof(double));
#else
		size_t size = mn * sizeof(double);
		memset(buf, 0, size);
		cblas_daxpby(mn, A_scalar, B_hessian, 1, B_scalar, A_hessian, 1);
		cblas_dger(CblasColMajor, m, n, 1.0, A_grad, 1, B_grad, 1, buf, m);
		cblas_daxpy(mn, 1.0, buf, 1, A_hessian, 1);
		memset(buf, 0, size);
		cblas_dger(CblasColMajor, m, n, 1.0, B_grad, 1, A_grad, 1, buf, m);
		cblas_daxpy(mn, 1.0, buf, 1, A_hessian, 1);
#endif
	}
	case 1: {
		double *A_grad = &A->data_[1];
		double *B_grad = &B->data_[1];
#if 0
		for (int i = 0; i < m; i++) {
			A_grad[i] = A_scalar * B_grad[i] + B_scalar * A_grad[i];
		}
#else
		cblas_daxpby(m, A_scalar, B_grad, 1, B_scalar, A_grad, 1);
#endif
	}
	case 0:
		A->data_[0] *= B_scalar;
	}
}

// A = A / B
void redukti_adouble_divide(redukti_adouble_t *A, redukti_adouble_t *B, redukti_adouble_t *temp1,
			    redukti_adouble_t *temp2)
{
	assert(A->order_ == B->order_ && A->vars_ == B->vars_);
	if (B->data_[0] == 0.0)
		// TODO division by zero
		return;
	// temp1 = B
	redukti_adouble_assign(temp1, B);
	// temp1 = pow(temp1, -1), i.e. inverse of B
	redukti_adouble_power(temp1, -1, temp2);
	// A = A * temp1;
	redukti_adouble_multiply(A, temp1, temp2);
}

// A = A^p
void redukti_adouble_power(redukti_adouble_t *A, double p, redukti_adouble_t *temp)
{
	int m = A->vars_, n = A->vars_;
	int mn = m * n;
	double tol = 1e-15;
	double value = A->data_[0];
	if (fabs(value) < tol) {
		if (value >= 0)
			value = tol;
		else if (value < 0)
			value = -tol;
	}
	double deriv = p * pow(value, (p - 1.0));
	switch (A->order_) {
	case 2: {
		double *A_grad = &A->data_[1];
		double *A_hessian = &A->data_[1 + m];
		double *temp_hessian = &temp->data_[0];
		double factor = p * (p - 1) * pow(value, p - 2.0);
#if 0
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				auto pos = i * m + j;
				temp_hessian[pos] = deriv * A_hessian[pos] + factor * A_grad[i] * A_grad[j];
			}
		}
		memcpy(A_hessian, temp_hessian, mn * sizeof(double));
#else
		size_t size = mn * sizeof(double);
		memset(temp_hessian, 0, size);
		cblas_dger(CblasColMajor, m, n, factor, A_grad, 1, A_grad, 1, temp_hessian, m);
		cblas_daxpby(mn, 1.0, temp_hessian, 1, deriv, A_hessian, 1);
#endif
	}
	case 1: {
		double *A_grad = &A->data_[1];
#if 0
		for (int i = 0; i < m; i++) {
			A_grad[i] = A_grad[i] * deriv;
		}
#else
		cblas_dscal(m, deriv, A_grad, 1);
#endif
	}
	case 0:
		A->data_[0] = pow(A->data_[0], p);
	}
}

// A = exp(A)
void redukti_adouble_exp(redukti_adouble_t *A, redukti_adouble_t *temp)
{
	double deriv = exp(A->data_[0]);
	int m = A->vars_, n = A->vars_;
	int mn = m * n;

	switch (A->order_) {
	case 2: {
		double *A_grad = &A->data_[1];
		double *A_hessian = &A->data_[1 + m];
		double *temp_hessian = &temp->data_[0];
#if 0
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				auto pos = i * m + j;
				temp_hessian[pos] = deriv * (A_hessian[pos] + A_grad[i] * A_grad[j]);
			}
		}
		memcpy(A_hessian, temp_hessian, mn * sizeof(double));
#else
		size_t size = mn * sizeof(double);
		memset(temp_hessian, 0, size);
		cblas_dger(CblasColMajor, m, n, deriv, A_grad, 1, A_grad, 1, temp_hessian, m);
		cblas_daxpby(mn, 1.0, temp_hessian, 1, deriv, A_hessian, 1);
#endif
	}
	case 1: {
		double *A_grad = &A->data_[1];
#if 0
		for (int i = 0; i < m; i++) {
			A_grad[i] = A_grad[i] * deriv;
		}
#else
		cblas_dscal(m, deriv, A_grad, 1);
#endif
	}
	case 0:
		A->data_[0] = deriv;
	}
}

// A = log(A)
void redukti_adouble_log(redukti_adouble_t *A, redukti_adouble_t *temp)
{
	double deriv = log(A->data_[0]);
	int m = A->vars_, n = A->vars_;
	int mn = m * n;
	switch (A->order_) {
	case 2: {
		double *A_grad = &A->data_[1];
		double *A_hessian = &A->data_[1 + m];
		double *temp_hessian = &temp->data_[0];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				auto pos = i * m + j;
				temp_hessian[pos] = A_hessian[pos] / A->data_[0] -
						    (A_grad[i] * A_grad[j]) / (A->data_[0] * A->data_[0]);
			}
		}
		memcpy(A_hessian, temp_hessian, mn * sizeof(double));
	}
	case 1: {
		double *A_grad = &A->data_[1];
		for (int i = 0; i < m; i++) {
			A_grad[i] = A_grad[i] / A->data_[0];
		}
	}
	case 0:
		A->data_[0] = deriv;
	}
}

// A = sin(A)
void redukti_adouble_sin(redukti_adouble_t *A, redukti_adouble_t *temp)
{
	double funval = sin(A->data_[0]);
	double deriv = cos(A->data_[0]);
	int m = A->vars_, n = A->vars_;
	int mn = m * n;
	switch (A->order_) {
	case 2: {
		double *A_grad = &A->data_[1];
		double *A_hessian = &A->data_[1 + m];
		double *temp_hessian = &temp->data_[0];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				auto pos = i * m + j;
				temp_hessian[pos] = deriv * A_hessian[pos] - funval * A_grad[i] * A_grad[j];
			}
		}
		memcpy(A_hessian, temp_hessian, mn * sizeof(double));
	}
	case 1: {
		double *A_grad = &A->data_[1];
		for (int i = 0; i < m; i++) {
			A_grad[i] = deriv * A_grad[i];
		}
	}
	case 0:
		A->data_[0] = funval;
	}
}

// A = cos(A)
void redukti_adouble_cos(redukti_adouble_t *A, redukti_adouble_t *temp)
{
	double funval = cos(A->data_[0]);
	double deriv = -sin(A->data_[0]);
	int m = A->vars_, n = A->vars_;
	int mn = m * n;
	switch (A->order_) {
	case 2: {
		double *A_grad = &A->data_[1];
		double *A_hessian = &A->data_[1 + m];
		double *temp_hessian = &temp->data_[0];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				auto pos = i * m + j;
				temp_hessian[pos] = deriv * A_hessian[pos] - funval * A_grad[i] * A_grad[j];
			}
		}
		memcpy(A_hessian, temp_hessian, mn * sizeof(double));
	}
	case 1: {
		double *A_grad = &A->data_[1];
		for (int i = 0; i < m; i++) {
			A_grad[i] = deriv * A_grad[i];
		}
	}
	case 0:
		A->data_[0] = funval;
	}
}

// A = tan(A)
void redukti_adouble_tan(redukti_adouble_t *A, redukti_adouble_t *temp)
{
	double funval = tan(A->data_[0]);
	double deriv = funval * funval + 1.0;
	int m = A->vars_, n = A->vars_;
	int mn = m * n;
	switch (A->order_) {
	case 2: {
		double *A_grad = &A->data_[1];
		double *A_hessian = &A->data_[1 + m];
		double *temp_hessian = &temp->data_[0];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				auto pos = i * m + j;
				temp_hessian[pos] =
				    deriv * A_hessian[pos] - (2 * funval * deriv) * A_grad[i] * A_grad[j];
			}
		}
		memcpy(A_hessian, temp_hessian, mn * sizeof(double));
	}
	case 1: {
		double *A_grad = &A->data_[1];
		for (int i = 0; i < m; i++) {
			A_grad[i] = deriv * A_grad[i];
		}
	}
	case 0:
		A->data_[0] = funval;
	}
}

// Dump contentx of variable
void redukti_adouble_dump(redukti_adouble_t *v, FILE *out, const char *desc)
{
	fprintf(out, "%s\n", desc);
	fprintf(out, "scalar = %f\n", v->data_[0]);
	fprintf(out, "vector = {\n");
	for (size_t i = 0; i < v->vars_; i++) {
		if (v->data_[1 + i] != 0.0)
			fprintf(out, "[%u] = %f\n", (unsigned int)i, v->data_[1 + i]);
	}
	fprintf(out, "}\n");
	fprintf(out, "matrix = {\n");
	for (size_t i = 0; i < v->vars_; i++) {
		for (size_t j = 0; j < v->vars_; j++) {
			size_t pos = j * v->vars_ + i;
			if (v->data_[1 + v->vars_ + pos] != 0.0)
				fprintf(out, "[%u=%u,%u] = %f\n", (unsigned int)pos, (unsigned int)i, (unsigned int)j,
					v->data_[1 + v->vars_ + pos]);
		}
	}
	fprintf(out, "}\n");
}

namespace redukti
{

////////////////////////////////////// TESTS

static int unequal(double got, double expected)
{
	if (fabs(got - expected) > 1e-6) {
		fprintf(stderr, "mismatch of values: expected %f, got %f\n", expected, got);
		return 1;
	}
	return 0;
}

static int unequal(double *got, double *expected, int n)
{
	int mismatched = 0;
	for (int i = 0; i < n; i++) {
		if (fabs(got[i] - expected[i]) > 1e-6) {
			fprintf(stderr, "mismatch of values at %d: expected %f, got %f\n", i, expected[i], got[i]);
			mismatched++;
		}
	}
	return mismatched;
}

static int test_basics()
{
	// No variables or derivatives so minimum size
	if (redukti_adouble_alloc_size(0, 0) != 16)
		return 1;
	// single variable and 1st order only
	if (redukti_adouble_alloc_size(1, 1) != 24)
		return 1;
	// single variable, 1st and 2nd order
	if (redukti_adouble_alloc_size(1, 2) != 32)
		return 1;

	size_t n = redukti_adouble_alloc_size(5, 1);
	if (n != 56)
		return 1;
	redukti_adouble_t *v1 = (redukti_adouble_t *)alloca(n);
	redukti_adouble_init(v1, 5, 1, 3, 4.2);
	if (redukti_adouble_get_value(v1) != 4.2)
		return 1;
	double expected1[] = {4.2, 0, 0, 0, 1.0, 0};
	if (unequal(v1->data_, expected1, v1->vars_ + 1))
		return 1;

	redukti_adouble_t *v2 = (redukti_adouble_t *)alloca(n);
	redukti_adouble_init(v2, 5, 1, 2, 3.1);
	if (redukti_adouble_get_value(v2) != 3.1)
		return 1;
	double expected2[] = {3.1, 0, 0, 1.0, 0, 0};
	if (unequal(v2->data_, expected2, v2->vars_ + 1))
		return 1;

	redukti_adouble_add(v1, v2, -1.0);
	double expected3[] = {4.2 - 3.1, 0, 0, -1.0, 1.0, 0};
	if (unequal(v1->data_, expected3, v1->vars_ + 1))
		return 1;

	redukti_adouble_scalar_multiply(v1, 3.0);
	double expected4[] = {(4.2 - 3.1) * 3, 0, 0, -3.0, 3.0, 0};
	if (unequal(v1->data_, expected4, v1->vars_ + 1))
		return 1;

	return 0;
}

// Compare two adoubles and return 0 if equal
// non zero if not
static int comp(redukti_adouble_t *x, redukti_adouble_t *y)
{
	if (x->vars_ != y->vars_)
		return 1;
	if (unequal(redukti_adouble_get_value(x), redukti_adouble_get_value(x)))
		return 1;
	for (int i = 0; i < x->vars_; i++) {
		if (unequal(redukti_adouble_get_derivative1(x, i), redukti_adouble_get_derivative1(x, i)))
			return 1;
	}
	for (int i = 0; i < x->vars_; i++) {
		for (int j = 0; j < x->vars_; j++) {
			if (unequal(redukti_adouble_get_derivative2(x, i, j), redukti_adouble_get_derivative2(x, i, j)))
				return 1;
		}
	}
	return 0;
}

static int test_mult()
{
	size_t n = redukti_adouble_alloc_size(3, 2);
	redukti_adouble_t *x = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *y = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *z = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *temp = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *expect = (redukti_adouble_t *)alloca(n);

	redukti_adouble_init(x, 3, 2, 0, 5.0);
	redukti_adouble_init(y, 3, 2, 1, 3.0);
	redukti_adouble_init(z, 3, 2, 2, 7.0);
	redukti_adouble_init(expect, 3, 2, 2, 0.0);

	// redukti_adouble_dump(x, stderr, "dump of x");

	redukti_adouble_multiply(x, x, temp);
	// redukti_adouble_dump(x, stderr, "dump of x");
	redukti_adouble_multiply(x, y, temp);
	redukti_adouble_multiply(x, z, temp);
	// redukti_adouble_dump(x, stderr, "dump of x*x*y*z");

	double expected[] = {525.0, 210.0, 175.0, 75.0, 42.0, 0, 0, 140.0, 0, 0, 60.0, 50.0, 0};
	for (int i = 0; i < (expect->vars_ + 1 + expect->vars_ * expect->vars_); i++) {
		expect->data_[i] = expected[i];
	}
	if (comp(x, expect))
		return 1;

	redukti_adouble_init(x, 3, 2, 0, 5.0);
	redukti_adouble_init(y, 3, 2, 1, 3.0);
	redukti_adouble_init(z, 3, 2, 2, 7.0);

	redukti_adouble_power(x, 2.0, temp);
	// redukti_adouble_dump(x, stderr, "dump of x");
	redukti_adouble_multiply(x, y, temp);
	redukti_adouble_multiply(x, z, temp);

	if (comp(x, expect))
		return 1;

	return 0;
}

static int test_div()
{
	size_t n = redukti_adouble_alloc_size(2, 2);
	redukti_adouble_t *x = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *y = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *tmp = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *tmpbuf = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *tmpbuf2 = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *expect = (redukti_adouble_t *)alloca(n);

	const double X = 5.0;
	const double Y = 3.0;

	redukti_adouble_init(x, 2, 2, 0, X);
	redukti_adouble_init(y, 2, 2, 1, Y);
	redukti_adouble_init(expect, 2, 2, -1, 0.0);

	redukti_adouble_multiply(x, x, tmpbuf);
	redukti_adouble_assign(tmp, y);

	redukti_adouble_multiply(y, y, tmpbuf);
	redukti_adouble_multiply(y, tmp, tmpbuf);

	redukti_adouble_divide(x, y, tmpbuf, tmpbuf2);

	// double expected_value = (X*X)/(Y*Y*Y);
	double expected1 = (2.0 * X) / (Y * Y * Y);
	double expected2 = -(3.0 * X * X) / (Y * Y * Y * Y);
	double expected11 = 2.0 / (Y * Y * Y);
	double expected12 = -(6 * X) / (Y * Y * Y * Y);
	double expected21 = expected12;
	double expected22 = (12.0 * X * X) / (Y * Y * Y * Y * Y);

	// redukti_adouble_dump(x, stderr, "dump of (x*x) / (y*y*y)");
	double expected[] = {0.925926, 0.370370, -0.925926, 0.074074, 0, -0.740741, 1.234568};
	for (int i = 0; i < (expect->vars_ + 1 + expect->vars_ * expect->vars_); i++) {
		expect->data_[i] = expected[i];
	}
	if (comp(x, expect))
		return 1;

	if (unequal(redukti_adouble_get_derivative1(x, 0), expected1)) {
		printf("Expected 1 %f got %f\n", expected1, redukti_adouble_get_derivative1(x, 0));
		return 1;
	}
	if (unequal(redukti_adouble_get_derivative1(x, 1), expected2)) {
		printf("Expected 2 %f got %f\n", expected2, redukti_adouble_get_derivative1(x, 1));
		return 1;
	}

	if (unequal(redukti_adouble_get_derivative2(x, 0, 0), expected11)) {
		printf("Expected 11 %f got %f\n", expected11, redukti_adouble_get_derivative2(x, 0, 0));
		return 1;
	}
	if (unequal(redukti_adouble_get_derivative2(x, 0, 1), expected12)) {
		printf("Expected 12 %f got %f\n", expected12, redukti_adouble_get_derivative2(x, 0, 1));
		return 1;
	}
	if (unequal(redukti_adouble_get_derivative2(x, 1, 0), expected21)) {
		printf("Expected 21 %f got %f\n", expected21, redukti_adouble_get_derivative2(x, 1, 0));
		return 1;
	}
	if (unequal(redukti_adouble_get_derivative2(x, 1, 1), expected22)) {
		printf("Expected 22 %f got %f\n", expected22, redukti_adouble_get_derivative2(x, 1, 1));
		return 1;
	}
	return 0;
}

static int test_hessian()
{
	size_t n = redukti_adouble_alloc_size(2, 2);
	redukti_adouble_t *x = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *y = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *tmp1 = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *tmp2 = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *tmp3 = (redukti_adouble_t *)alloca(n);
	redukti_adouble_t *tmpbuf = (redukti_adouble_t *)alloca(n);

	const double X = 5.0;
	const double Y = 3.0;

	// x=5, y = 3
	// 3*x^2*y^2 - 4*x*y^3

	redukti_adouble_init(x, 2, 2, 0, X);
	redukti_adouble_init(y, 2, 2, 1, Y);

	// tmp1 = x^2
	redukti_adouble_assign(tmp1, x);
	redukti_adouble_power(tmp1, 2.0, tmpbuf);

	// tmp2 = y^2
	redukti_adouble_assign(tmp2, y);
	redukti_adouble_power(tmp2, 2.0, tmpbuf);

	// tmp1 = tmp1 * tmp2
	redukti_adouble_multiply(tmp1, tmp2, tmpbuf);
	// tmp1 = 3*tmp1
	redukti_adouble_scalar_multiply(tmp1, 3.0);

	// tmp3 = y^3
	redukti_adouble_assign(tmp3, y);
	redukti_adouble_power(tmp3, 3.0, tmpbuf);

	// tmp3 = tmp3 * x * 4
	redukti_adouble_multiply(tmp3, x, tmpbuf);
	redukti_adouble_scalar_multiply(tmp3, 4.0);

	// tmp1 = tmp1 - tmp3
	redukti_adouble_add(tmp1, tmp3, -1.0);

	// d/dx = 6*x*y^2 - 4*y^3
	double expected1 = 6.0 * X * Y * Y - 4.0 * Y * Y * Y;
	// d/dy = 6*x^2*y - 12*x*y^2
	double expected2 = 6.0 * X * X * Y - 12.0 * X * Y * Y;
	if (unequal(redukti_adouble_get_derivative1(tmp1, 0), expected1))
		return 1;
	if (unequal(redukti_adouble_get_derivative1(tmp1, 1), expected2))
		return 1;

	// d^2/dx^2 = 6*y^2
	double expected11 = 6.0 * Y * Y;
	// d^2/dxy = 12*x*y - 12*y^2
	double expected12 = 12.0 * X * Y - 12.0 * Y * Y;
	double expected21 = expected12;
	// d^2/dy^2 = 6*x^2 - 24*x*y
	double expected22 = 6.0 * X * X - 24.0 * X * Y;

	// redukti_adouble_dump(tmp1, stdout, "Dump of tmp1");
	if (unequal(redukti_adouble_get_derivative2(tmp1, 0, 0), expected11)) {
		printf("Expected 11 %f got %f\n", expected11, redukti_adouble_get_derivative2(tmp1, 0, 0));
		return 1;
	}
	if (unequal(redukti_adouble_get_derivative2(tmp1, 0, 1), expected12)) {
		printf("Expected 12 %f got %f\n", expected12, redukti_adouble_get_derivative2(tmp1, 0, 1));
		return 1;
	}
	if (unequal(redukti_adouble_get_derivative2(tmp1, 1, 0), expected21)) {
		printf("Expected 21 %f got %f\n", expected21, redukti_adouble_get_derivative2(tmp1, 1, 0));
		return 1;
	}
	if (unequal(redukti_adouble_get_derivative2(tmp1, 1, 1), expected22)) {
		printf("Expected 22 %f got %f\n", expected22, redukti_adouble_get_derivative2(tmp1, 1, 1));
		return 1;
	}
	return 0;
}

int test_autodiff()
{
	int rc = test_basics();
	rc += test_mult();
	rc += test_div();
	rc += test_hessian();
	if (rc == 0)
		printf("Autodiff Tests OK\n");
	else
		printf("Autodiff Tests FAILED\n");
	return rc == 0 ? 0 : 1;
}
} // namespace redukti
