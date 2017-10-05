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
#ifndef _REDUKTI_AUTODIFF_H_
#define _REDUKTI_AUTODIFF_H_

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>

/*
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

#ifdef __cplusplus
extern "C" {
#endif

// WARNING
//
// This is a low level module that must be used with care.
// In general this module requires the caller to allocate memory
// correctly - as it assumes that all supplied arguments are
// properly sized and allocated.

/* autodiff variable */
struct redukti_adouble_t {
	// derivative order
	uint32_t order_ : 2;
	// number of variables
	uint32_t vars_;
	// data
	double data_[1];

	redukti_adouble_t(const redukti_adouble_t &) = delete;
	redukti_adouble_t &operator=(const redukti_adouble_t &) = delete;
};

// Compute memory requirement for given number of variables and order
// Supported orders are 0,1,2.
size_t redukti_adouble_alloc_size(int vars, int order);

// Initialize A; caller must have allocated memory of correct
// size.
void redukti_adouble_init(redukti_adouble_t *A, int n_vars, int order, int var, double v);

// A = B
// must be same size
void redukti_adouble_assign(redukti_adouble_t *A, const redukti_adouble_t *B);

// A = A + alpha*B
void redukti_adouble_add(redukti_adouble_t *A, redukti_adouble_t *B, double alpha);

// A = A*scalar
void redukti_adouble_scalar_multiply(redukti_adouble_t *A, double alpha);

// A = A*B
// A = A*A also works
// temp must be same size as A
void redukti_adouble_multiply(redukti_adouble_t *A, redukti_adouble_t *B, redukti_adouble_t *temp);

// A = A/B
// temp1, temp2 must be same size as A
void redukti_adouble_divide(redukti_adouble_t *A, redukti_adouble_t *B, redukti_adouble_t *temp1,
			    redukti_adouble_t *temp2);

// A = exp(A)
// temp must be same size as A
void redukti_adouble_exp(redukti_adouble_t *A, redukti_adouble_t *temp);

// A = log(A)
// temp must be same size as A
void redukti_adouble_log(redukti_adouble_t *A, redukti_adouble_t *temp);

// A = A^p
// temp must be same size as A
void redukti_adouble_power(redukti_adouble_t *A, double p, redukti_adouble_t *temp);

// A = abs(A)
void redukti_adouble_abs(redukti_adouble_t *A);

// A = sin(A)
// temp must be same size as A
void redukti_adouble_sin(redukti_adouble_t *A, redukti_adouble_t *temp);

// A = cos(A)
// temp must be same size as A
void redukti_adouble_cos(redukti_adouble_t *A, redukti_adouble_t *temp);

// A = tan(A)
// temp must be same size as A
void redukti_adouble_tan(redukti_adouble_t *A, redukti_adouble_t *temp);

// Dumps contents of A
void redukti_adouble_dump(redukti_adouble_t *v, FILE *out, const char *desc);

// A = A + alpha
static inline void redukti_adouble_scalar_add(redukti_adouble_t *A, double alpha) { A->data_[0] += alpha; }

// Get A's value
static inline double redukti_adouble_get_value(redukti_adouble_t *x) { return x->data_[0]; }

// Get first derivative
static inline double redukti_adouble_get_derivative1(redukti_adouble_t *x, int parameter)
{
	if (x->order_ == 0)
		return 0.0;
	assert(parameter >= 0 && parameter < x->vars_);
	return x->data_[parameter + 1];
}

// Get second derivative
static inline double redukti_adouble_get_derivative2(redukti_adouble_t *x, int parameter1, int parameter2)
{
	if (x->order_ < 2)
		return 0.0;
	int m = x->vars_;
	assert(parameter1 >= 0 && parameter1 < m && parameter2 >= 0 && parameter2 < m);
	double *matrix = &x->data_[1 + m];
	int pos1 = parameter2 * m + parameter1;
	return matrix[pos1];
}

static inline void redukti_adouble_set_value(redukti_adouble_t *x, double v) { x->data_[0] = v; }

static inline void redukti_adouble_set_derivative1(redukti_adouble_t *x, int parameter, double v)
{
	if (x->order_ == 0)
		return;
	assert(parameter >= 0 && parameter < x->vars_);
	x->data_[parameter + 1] = v;
}

static inline void redukti_adouble_set_derivative2(redukti_adouble_t *x, int parameter1, int parameter2, double v)
{
	if (x->order_ < 2)
		return;
	int m = x->vars_;
	assert(parameter1 >= 0 && parameter1 < m && parameter2 >= 0 && parameter2 < m);
	double *matrix = &x->data_[1 + m];
	int pos = parameter1 * m + parameter2;
	matrix[pos] = v;
}

#ifdef __cplusplus
}
#endif

namespace redukti
{

extern int test_autodiff();
}

#endif
