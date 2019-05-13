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
#ifndef _REDUKTI_AUTODIFF_HELPERS_H
#define _REDUKTI_AUTODIFF_HELPERS_H

#include <allocators.h>
#include <autodiff.h>

#include <cmath>

namespace redukti
{

// As the size of adouble is dynamically determined the standard
// vector type cannot be used. We define a special array of adoubles
// that makes it easy to handle an array of adouble values
class ADoubleArray
{
	public:
	ADoubleArray(int vars, int order, int n, Allocator *A) : A_(A), n_(n)
	{
		assert(n > 0);
		assert(vars >= 1 && order >= 0 && order <= 2);
		// TODO do we need to check that the size is correctly aligned?
		element_size_ = redukti_adouble_alloc_size(vars, order);
		auto allocated_size = n * element_size_;
		array_ = (char *)A_->allocate(allocated_size);
		memset(array_, 0, allocated_size);
	}

	~ADoubleArray() { A_->deallocate(array_); }

	redukti_adouble_t &operator[](int i)
	{
		assert(i >= 0 && i < n_);
		char *p = array_ + i * element_size_;
		redukti_adouble_t *element = reinterpret_cast<redukti_adouble_t *>(p);
		return *element;
	}

	const redukti_adouble_t &operator[](int i) const
	{
		assert(i >= 0 && i < n_);
		char *p = array_ + i * element_size_;
		redukti_adouble_t *element = reinterpret_cast<redukti_adouble_t *>(p);
		return *element;
	}

	private:
	Allocator *A_;
	// Memory allocated for the array
	char *array_;
	// Size of each array element
	int element_size_;
	// Number of array elements
	int n_;
};

// This mirrors the ADoubleArray type but deals with
// doubles.
class DoubleArray
{
	public:
	DoubleArray(int vars, int order, int n, Allocator *A) : A_(A), n_(n)
	{
		//    assert(vars >= 1 && order == 0);
		auto element_size_ = sizeof(double);
		auto allocated_size = n * element_size_;
		array_ = (double *)A_->allocate(allocated_size);
		std::fill_n(array_, n_, 0.0);
	}

	~DoubleArray() { A_->deallocate(array_); }

	double &operator[](int i)
	{
		assert(i >= 0 && i < n_);
		return array_[i];
	}

	const double &operator[](int i) const
	{
		assert(i >= 0 && i < n_);
		return array_[i];
	}

	private:
	Allocator *A_;
	double *array_;
	int n_;
};

// General case not implemented
// We use the AutoDiffTraits template to handle
// differences between plain double and adouble implementations
template <typename ad_type> struct AutoDiffTraits;

template <> struct AutoDiffTraits<redukti_adouble_t> {
	typedef redukti_adouble_t ad_type;
	typedef double value_type;
	typedef ADoubleArray array_type;
	typedef Deleter<ad_type> deleter_type;
	typedef std::unique_ptr<ad_type, deleter_type> ptr_type;

	static void init(ad_type &x, int n_vars, int order, int var, double v)
	{
		redukti_adouble_init(&x, n_vars, order, var, v);
	}
	static void assign(ad_type &A, const ad_type &B) { redukti_adouble_assign(&A, &B); }
	static void add(ad_type &A, ad_type &B, double alpha) { redukti_adouble_add(&A, &B, alpha); }
	static void scalar_multiply(ad_type &A, double alpha) { redukti_adouble_scalar_multiply(&A, alpha); }
	static void scalar_add(ad_type &A, double alpha) { redukti_adouble_scalar_add(&A, alpha); }
	static void multiply(ad_type &A, ad_type &B, ad_type &temp) { redukti_adouble_multiply(&A, &B, &temp); }
	static void divide(ad_type &A, ad_type &B, ad_type &temp1, ad_type &temp2)
	{
		redukti_adouble_divide(&A, &B, &temp1, &temp2);
	}
	static void exp(ad_type &A, ad_type &temp) { redukti_adouble_exp(&A, &temp); }
	static void pow(ad_type &A, double p, ad_type &temp) { redukti_adouble_power(&A, p, &temp); }
	static void abs(ad_type &A) { redukti_adouble_get_value(&A); }
	// min function
	static ad_type &min(ad_type &x, ad_type &y)
	{
		if (redukti_adouble_get_value(&x) <= redukti_adouble_get_value(&y))
			return x;
		return y;
	}
	// max function
	static ad_type &max(ad_type &x, ad_type &y)
	{
		if (redukti_adouble_get_value(&x) < redukti_adouble_get_value(&y))
			return y;
		return x;
	}
	static value_type value(ad_type &A) { return redukti_adouble_get_value(&A); }
	static value_type value(const ptr_type &A) { return redukti_adouble_get_value(A.get()); };
	static ad_type &reference(const ptr_type &A) { return *A; }
	static ptr_type allocate(Allocator *A, int n_vars, int order, int var, double v)
	{
		auto element_size = alloc_size(n_vars, order);
		ad_type *p = reinterpret_cast<ad_type *>(A->allocate(element_size));
		assert(p != nullptr);
		redukti_adouble_init(p, n_vars, order, var, v);
		return ptr_type(p, deleter_type(A));
	}
	static size_t alloc_size(int n_vars, int order) { return redukti_adouble_alloc_size(n_vars, order); }
};

template <> struct AutoDiffTraits<double> {
	typedef double ad_type;
	typedef double value_type;
	typedef DoubleArray array_type;
	typedef double ptr_type;

	static void init(ad_type &x, int, int, int, double v) { x = v; }
	static void assign(ad_type &A, const ad_type &B) { A = B; }
	static void add(ad_type &A, ad_type &B, double alpha) { A += alpha * B; }
	static void scalar_multiply(ad_type &A, double alpha) { A *= alpha; }
	static void scalar_add(ad_type &A, double alpha) { A += alpha; }
	static void multiply(ad_type &A, ad_type &B, ad_type &temp) { A *= B; }
	static void divide(ad_type &A, ad_type &B, ad_type &temp1, ad_type &temp2) { A /= B; }
	static void exp(ad_type &A, ad_type &temp1) { A = std::exp(A); }
	static void pow(ad_type &A, double p, ad_type &temp1) { A = std::pow(A, p); }
	static void abs(ad_type &A) { A = std::abs(A); }
	// following needs to use references
	static double &min(double &x, double &y) { return x <= y ? x : y; }
	static double &max(double &x, double &y) { return x >= y ? x : y; }
	static value_type value(ad_type &A) { return A; }
	static value_type &reference(ptr_type &A) { return A; }
	static ptr_type allocate(Allocator *, int, int, int, double v) { return v; }
	static size_t alloc_size(int n_vars, int order) { return sizeof(double); }
};

} // namespace redukti

#endif
