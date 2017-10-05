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
 * Portions derived from Quantlib.
 * License: http://quantlib.org/license.shtml
 */
#include <interpolators.h>

#include <allocators.h>
#include <autodiff.h>

#include <internal/autodiff_helpers.h>

#include <assert.h>
#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>

namespace redukti
{

static void adouble_set_derivatives1(redukti_adouble_t &ad, int index1, double value1, double c1, int index2,
				     double value2, double c2, int index3, double value3, double c3)
{
	double t1 = c1 * value1;
	double t2 = c2 * value2;
	double t3 = c3 * value3;
	redukti_adouble_set_value(&ad, t1 + t2 + t3);
	if (ad.order_ > 0 && index1 >= 0) {
		redukti_adouble_set_derivative1(&ad, index1, c1);
	}
	if (ad.order_ > 0 && index2 >= 0) {
		redukti_adouble_set_derivative1(&ad, index2, redukti_adouble_get_derivative1(&ad, index2) + c2);
	}
	if (ad.order_ > 0 && index3 >= 0) {
		redukti_adouble_set_derivative1(&ad, index3, redukti_adouble_get_derivative1(&ad, index3) + c3);
	}
}

static inline size_t locate_x(const double *xBegin_, const double *xEnd_, double x)
{
	if (x < *xBegin_)
		return 0;
	else if (x > *(xEnd_ - 1))
		return xEnd_ - xBegin_ - 2;
	else
		return std::upper_bound(xBegin_, xEnd_ - 1, x) - xBegin_ - 1;
}

enum BoundaryCondition {
	//! Make second(-last) point an inactive knot
	NotAKnot,
	//! Match value of end-slope
	FirstDerivative,
	//! Match value of second derivative at end
	SecondDerivative
};

inline bool close(double x, double y) { return std::fabs(x - y) < 1e-10; }
static int equal(double got, double expected)
{
	if (std::fabs(expected) > 10) {
		if (std::fabs(got - expected) > 1e-2) {
			fprintf(stderr, "mismatch of values: expected %f, got %f\n", expected, got);
			return 0;
		}
	} else if (std::fabs(got - expected) > 1e-6) {
		fprintf(stderr, "mismatch of values: expected %f, got %f\n", expected, got);
		return 0;
	}
	return 1;
}

class AbstractInterpolator : public Interpolator
{
	typedef AutoDiffTraits<redukti_adouble_t> AD;
	typedef typename AD::ptr_type ad_type_ptr;

      public:
	virtual ~AbstractInterpolator() {}
	virtual double interpolate(double x) = 0;

	virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
	interpolate_with_sensitivities(double x, FixedRegionAllocator *A) = 0;

	/**
	 * If underlying values have changed, this
	 * method can be called to reinitialise the
	 * interpolator.
	 */
	virtual void update() = 0;
	/**
	 * Only available on Monotone Convex interpolator as it is an
	 * interest rate aware interpolator - for everything else
	 * an exception will be thrown.
	 */
	virtual double forward(double x)
	{
		assert(false);
		return 0.0;
	}
	virtual bool coefficients(std::vector<double> &cvec) { return false; }

	virtual double primitive(double) { return 0.0; }
	virtual double derivative(double) { return 0.0; }
	virtual double derivative2(double) { return 0.0; }
	/**
	 * Return the interpolator type
	 */
	virtual InterpolatorType type() const = 0;

	virtual size_t size() const = 0;
	virtual double *xvalues() const = 0;
	virtual double *yvalues() const = 0;
	virtual int order() const = 0;
	virtual void get_options(InterpolationOptions &optons) const = 0;
	virtual std::unique_ptr<Interpolator, Deleter<Interpolator>>
	get_interpolator(double *xa, double *ya, unsigned int size, Allocator *A,
			 const InterpolationOptions &options) const = 0;

	virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
	interpolate_with_numeric_sensitivities(double x, FixedRegionAllocator *A)
	{
		auto v = AD::allocate(A, size(), order(), -1, 0.0);
		redukti_adouble_set_value(v.get(), interpolate(x));

		InterpolationOptions options;
		get_options(options);
		options.differentiation_order = 0;

		if (order() >= 1) {
			double *tmp_y = reinterpret_cast<double *>(alloca(sizeof(double) * size()));
			const double h = 0.0001;
			std::copy(yvalues(), yvalues() + size(), tmp_y);
			double *tmp_z = reinterpret_cast<double *>(alloca(sizeof(double) * size()));
			std::copy(yvalues(), yvalues() + size(), tmp_z);

			DynamicRegionAllocator alloc;
			auto I = get_interpolator(xvalues(), tmp_y, size(), &alloc, options);
			auto J = get_interpolator(xvalues(), tmp_z, size(), &alloc, options);
			for (int i = 0; i < size(); i++) {
				tmp_y[i] = tmp_y[i] + h;
				tmp_z[i] = tmp_z[i] - h;
				I->update();
				J->update();
				double d = (I->interpolate(x) - J->interpolate(x)) / (2.0 * h);
				redukti_adouble_set_derivative1(v.get(), i, d);
				tmp_y[i] = tmp_y[i] - h;
				tmp_z[i] = tmp_z[i] + h;
			}
		}
		if (order() == 2) {
			// see http://en.wikipedia.org/wiki/Finite_difference,
			// for how the calculations can be optimized
			double *tmp_y_ip_jp = reinterpret_cast<double *>(alloca(sizeof(double) * size()));
			double *tmp_y_ip_jn = reinterpret_cast<double *>(alloca(sizeof(double) * size()));
			double *tmp_y_in_jp = reinterpret_cast<double *>(alloca(sizeof(double) * size()));
			double *tmp_y_in_jn = reinterpret_cast<double *>(alloca(sizeof(double) * size()));

			const double h = 0.00001;

			std::copy(yvalues(), yvalues() + size(), tmp_y_ip_jp);
			std::copy(yvalues(), yvalues() + size(), tmp_y_ip_jn);
			std::copy(yvalues(), yvalues() + size(), tmp_y_in_jp);
			std::copy(yvalues(), yvalues() + size(), tmp_y_in_jn);

			DynamicRegionAllocator alloc;
			auto IP_JP = get_interpolator(xvalues(), tmp_y_ip_jp, size(), &alloc, options);
			auto IP_JN = get_interpolator(xvalues(), tmp_y_ip_jn, size(), &alloc, options);
			auto IN_JP = get_interpolator(xvalues(), tmp_y_in_jp, size(), &alloc, options);
			auto IN_JN = get_interpolator(xvalues(), tmp_y_in_jn, size(), &alloc, options);

			for (int i = 0; i < size(); i++) {
				for (int j = 0; j < size(); j++) {
					*(tmp_y_ip_jp + i) = *(tmp_y_ip_jp + i) + h;
					*(tmp_y_ip_jp + j) = *(tmp_y_ip_jp + j) + h;

					*(tmp_y_ip_jn + i) = *(tmp_y_ip_jn + i) + h;
					*(tmp_y_ip_jn + j) = *(tmp_y_ip_jn + j) - h;

					*(tmp_y_in_jp + i) = *(tmp_y_in_jp + i) - h;
					*(tmp_y_in_jp + j) = *(tmp_y_in_jp + j) + h;

					*(tmp_y_in_jn + i) = *(tmp_y_in_jn + i) - h;
					*(tmp_y_in_jn + j) = *(tmp_y_in_jn + j) - h;

					IP_JP->update();
					IP_JN->update();
					IN_JP->update();
					IN_JN->update();

					double result = (IP_JP->interpolate(x) - IP_JN->interpolate(x) -
							 IN_JP->interpolate(x) + IN_JN->interpolate(x)) /
							(2 * h * h);

					redukti_adouble_set_derivative2(v.get(), i, j, result);

					*(tmp_y_ip_jp + i) = *(tmp_y_ip_jp + i) - h;
					*(tmp_y_ip_jp + j) = *(tmp_y_ip_jp + j) - h;

					*(tmp_y_ip_jn + i) = *(tmp_y_ip_jn + i) - h;
					*(tmp_y_ip_jn + j) = *(tmp_y_ip_jn + j) + h;

					*(tmp_y_in_jp + i) = *(tmp_y_in_jp + i) + h;
					*(tmp_y_in_jp + j) = *(tmp_y_in_jp + j) - h;

					*(tmp_y_in_jn + i) = *(tmp_y_in_jn + i) + h;
					*(tmp_y_in_jn + j) = *(tmp_y_in_jn + j) + h;
				}
			}
		}
		return v;
	}
};

/**
 * Implementation of Linear Interpolator
 * Dev: Allocator usage verified.
 */
class LinearInterpolator : public AbstractInterpolator
{
	typedef AutoDiffTraits<redukti_adouble_t> AD;
	typedef typename AD::ptr_type ad_type_ptr;

	double *x_, *y_, *xend_;
	unsigned int size_;
	Allocator *A_;
	double *primitiveConst_, *s_;
	int order_;

      public:
	LinearInterpolator(double *x, double *y, unsigned int size, int order, Allocator *alloc)
	    : x_(x), y_(y), xend_(x + size), size_(size), A_(alloc), primitiveConst_(nullptr), s_(nullptr),
	      order_(order)
	{
		primitiveConst_ = new (*A_) double[size];
		s_ = new (*A_) double[size];
		update_internal();
	}
	int locate(double x) const { return locate_x(x_, xend_, x); }
	virtual ~LinearInterpolator()
	{
		if (primitiveConst_)
			A_->deallocate(primitiveConst_);
		if (s_)
			A_->deallocate(s_);
	}
	virtual double primitive(double x) override final
	{
		int i = locate(x);
		double dx = x - x_[i];
		return primitiveConst_[i] + dx * (y_[i] + 0.5 * dx * s_[i]);
	}
	virtual double interpolate(double x) override final
	{
		int i = locate(x);
		return y_[i] + (x - x_[i]) * s_[i];
	}
	virtual ad_type_ptr interpolate_with_sensitivities(double x, FixedRegionAllocator *A) override final
	{
		int i = locate(x);
		double tmp = (x - x_[i]) / (x_[i + 1] - x_[i]);
		auto ptr = AD::allocate(A, size_, order_, -1, 0.0);
		adouble_set_derivatives1(*ptr, i, y_[i], 1.0, i + 1, y_[i + 1], tmp, i, y_[i], -tmp);
		return ptr;
	}
	virtual bool has_sensitivities() const { return true; }
	// could implement in a matrix format, all dates vs. all knotpoints
	virtual double derivative(double x) override final
	{
		int i = locate(x);
		return s_[i];
	}
	virtual double derivative2(double x) override final { return 0; }
	virtual void update() override final { update_internal(); }
	void update_internal()
	{
		primitiveConst_[0] = 0.0;
		for (unsigned int i = 1; i < size_; ++i) {
			auto dx = x_[i] - x_[i - 1];
			s_[i - 1] = (y_[i] - y_[i - 1]) / dx;
			primitiveConst_[i] = primitiveConst_[i - 1] + dx * (y_[i - 1] + 0.5 * dx * s_[i - 1]);
		}
	}
	virtual InterpolatorType type() const override final { return InterpolatorType::LINEAR; }
	virtual size_t size() const override final { return size_; };
	virtual double *xvalues() const override final { return x_; }
	virtual double *yvalues() const override final { return y_; }
	virtual int order() const override final { return order_; }
	virtual void get_options(InterpolationOptions &options) const override final
	{
		options.differentiation_order = order_;
	}
	virtual std::unique_ptr<Interpolator, Deleter<Interpolator>>
	get_interpolator(double *xa, double *ya, unsigned int size, Allocator *A,
			 const InterpolationOptions &options) const override final
	{
		return make_interpolator(InterpolatorType::LINEAR, xa, ya, size, A, options);
	}

      private:
	LinearInterpolator(const LinearInterpolator &) = delete;
	LinearInterpolator &operator=(const LinearInterpolator &) = delete;
};

// Dev: allocator usage verified
class LogInterpolator : public AbstractInterpolator
{
      private:
	typedef AutoDiffTraits<redukti_adouble_t> AD;
	typedef typename AD::ptr_type ad_type_ptr;
	std::unique_ptr<Interpolator, Deleter<Interpolator>> I_;
	double *x_;
	double *y_;
	double *logy_;
	size_t size_;
	Allocator *A_;
	InterpolatorType type_;

      public:
	LogInterpolator(double *x, double *y, unsigned int size, InterpolatorType subtype, InterpolatorType type,
			Allocator *A, const InterpolationOptions &data = InterpolationOptions())
	    : x_(x), y_(y), logy_(nullptr), size_(size), A_(A), type_(type)
	{
		logy_ = new (*A_) double[size];
		for (unsigned int i = 0; i < size_; ++i) {
			logy_[i] = std::log(y_[i]);
		}
		I_ = make_interpolator(subtype, x_, logy_, size_, A_, data);
	}
	virtual ~LogInterpolator() { A_->deallocate(logy_); }
	void update_internal()
	{
		for (unsigned int i = 0; i < size_; ++i) {
			logy_[i] = std::log(y_[i]);
		}
		I_->update();
	}
	virtual void update() override final { update_internal(); }
	virtual double interpolate(double x) override final { return std::exp(I_->interpolate(x)); }
	virtual ad_type_ptr interpolate_with_sensitivities(double x, FixedRegionAllocator *A) override final
	{
		auto v = I_->interpolate_with_sensitivities(x, A);
		// It turns out that when using log interpolators
		// we need the sensitivities generated by the
		// underlying interpolator rather than
		// So although we convert back the interpolated
		// value, we need to leave the sensitivities as they
		// are. If we were to apply exp() here then
		// the caller will need to undo using log() - so we
		// avoid unnecessary operations. However this approach
		// seems "wrong" somehow, ... there must be an incorrect
		// abstraction somewhere FIXME
#if 0
		FixedRegionAllocatorGuard guard(A);
		redukti_adouble_t *temp = (redukti_adouble_t *)A->allocate(
		    redukti_adouble_alloc_size(v->vars_, v->order_));
		redukti_adouble_exp(v.get(), temp);
#else
		redukti_adouble_set_value(v.get(), std::exp(redukti_adouble_get_value(v.get())));
#endif
		return v;
	}
	virtual double derivative(double x) override final { return interpolate(x) * I_->derivative(x); }
	virtual double derivative2(double x) override final
	{
		return derivative(x) * I_->derivative(x) + interpolate(x) * I_->derivative2(x);
	}
	virtual double primitive(double x) override final
	{
		assert(false);
		return 0;
	}
	virtual InterpolatorType type() const override final { return type_; }
	virtual size_t size() const override final { return size_; };
	virtual double *xvalues() const override final { return x_; }
	virtual double *yvalues() const override final { return y_; }
	virtual int order() const override final { return I_->order(); }
	virtual void get_options(InterpolationOptions &options) const override final { I_->get_options(options); }
	virtual std::unique_ptr<Interpolator, Deleter<Interpolator>>
	get_interpolator(double *xa, double *ya, unsigned int size, Allocator *A,
			 const InterpolationOptions &options) const override final
	{
		return make_interpolator(I_->type(), xa, ya, size, A, options);
	}

      private:
	LogInterpolator(const LogInterpolator &) = delete;
	LogInterpolator &operator=(const LogInterpolator &) = delete;
};

// Dev: FIXME Allocator required
class ForwardFlatInterpolation : public AbstractInterpolator
{
	typedef AutoDiffTraits<redukti_adouble_t> AD;
	typedef typename AD::ptr_type ad_type_ptr;
	double *x_, *xend_, *y_;
	unsigned int size_;
	InterpolatorType type_;
	double *primitive_;
	int order_;

      public:
	ForwardFlatInterpolation(double *x, double *y, size_t size, int order)
	    : x_(x), xend_(x + size), y_(y), size_(size), type_(InterpolatorType::FLAT_RIGHT), primitive_(nullptr),
	      order_(order)
	{
		primitive_ = new double[size_];
		std::fill_n(primitive_, size_, 0.0);
		update_internal();
	}
	virtual ~ForwardFlatInterpolation() { delete[] primitive_; }
	int locate(double x) const { return locate_x(x_, xend_, x); }
	virtual void update() override final { update_internal(); }
	void update_internal()
	{
		primitive_[0] = 0.0;
		for (size_t i = 1; i < size_; i++) {
			double dx = this->x_[i] - this->x_[i - 1];
			primitive_[i] = primitive_[i - 1] + dx * this->y_[i - 1];
		}
	}
	virtual double interpolate(double x) override final
	{
		if (x >= this->x_[size_ - 1])
			return this->y_[size_ - 1];

		auto i = this->locate(x);
		return this->y_[i];
	}
	virtual ad_type_ptr interpolate_with_sensitivities(double x, FixedRegionAllocator *A) override final
	{
		if (x >= this->x_[size_ - 1])
			return AD::allocate(A, size_, order_, size_ - 1, this->y_[size_ - 1]);
		auto i = this->locate(x);
		return AD::allocate(A, size_, order_, i, this->y_[i]);
	}
	virtual double primitive(double x) override final
	{
		auto i = this->locate(x);
		double dx = x - this->x_[i];
		return primitive_[i] + dx * this->y_[i];
	}
	virtual double derivative(double) override final { return 0.0; }
	virtual double derivative2(double) override final { return 0.0; }
	virtual InterpolatorType type() const override final { return type_; }
	virtual size_t size() const override final { return size_; };
	virtual double *xvalues() const override final { return x_; }
	virtual double *yvalues() const override final { return y_; }
	virtual int order() const override final { return order_; }
	virtual void get_options(InterpolationOptions &options) const override final
	{
		options.differentiation_order = order_;
	}
	virtual std::unique_ptr<Interpolator, Deleter<Interpolator>>
	get_interpolator(double *xa, double *ya, unsigned int size, Allocator *A,
			 const InterpolationOptions &options) const override final
	{
		return make_interpolator(type(), xa, ya, size, A, options);
	}
};

// Dev: FIXME Allocator required
class BackwardFlatInterpolation : public AbstractInterpolator
{
	typedef AutoDiffTraits<redukti_adouble_t> AD;
	typedef typename AD::ptr_type ad_type_ptr;
	double *x_, *xend_, *y_;
	unsigned int size_;
	InterpolatorType type_;
	double *primitive_;
	int order_;

      public:
	BackwardFlatInterpolation(double *x, double *y, size_t size, int order)
	    : x_(x), xend_(x + size), y_(y), size_(size), type_(InterpolatorType::FLAT_LEFT), primitive_(nullptr),
	      order_(order)
	{
		primitive_ = new double[size_];
		std::fill_n(primitive_, size_, 0.0);
		update_internal();
	}
	virtual ~BackwardFlatInterpolation() { delete[] primitive_; }
	virtual void update() override final { update_internal(); }
	void update_internal()
	{
		size_t n = size_;
		primitive_[0] = 0.0;
		for (size_t i = 1; i < n; i++) {
			double dx = this->x_[i] - this->x_[i - 1];
			primitive_[i] = primitive_[i - 1] + dx * this->y_[i];
		}
	}
	int locate(double x) const { return locate_x(x_, xend_, x); }
	virtual double interpolate(double x) override final
	{
		if (x <= this->x_[0])
			return this->y_[0];
		int i = this->locate(x);
		if (x == this->x_[i])
			return this->y_[i];
		else
			return this->y_[i + 1];
	}
	virtual ad_type_ptr interpolate_with_sensitivities(double x, FixedRegionAllocator *A) override final
	{
		if (x <= this->x_[0])
			return AD::allocate(A, size_, order_, 0, this->y_[0]);
		int i = this->locate(x);
		if (x == this->x_[i])
			return AD::allocate(A, size_, order_, i, this->y_[i]);
		else
			return AD::allocate(A, size_, order_, i + 1, this->y_[i + 1]);
	}
	virtual double primitive(double x) override final
	{
		int i = this->locate(x);
		double dx = x - this->x_[i];
		return primitive_[i] + dx * this->y_[i + 1];
	}
	virtual double derivative(double) override final { return 0.0; }
	virtual double derivative2(double) override final { return 0.0; }
	virtual InterpolatorType type() const override final { return type_; }
	virtual size_t size() const override final { return size_; };
	virtual double *xvalues() const override final { return x_; }
	virtual double *yvalues() const override final { return y_; }
	virtual int order() const override final { return order_; }
	virtual void get_options(InterpolationOptions &options) const override final
	{
		options.differentiation_order = order_;
	}
	virtual std::unique_ptr<Interpolator, Deleter<Interpolator>>
	get_interpolator(double *xa, double *ya, unsigned int size, Allocator *A,
			 const InterpolationOptions &options) const override final
	{
		return make_interpolator(type(), xa, ya, size, A, options);
	}
};

template <typename T> struct CubicInterpolator {
	typedef AutoDiffTraits<T> AD;
	typedef typename AD::ptr_type ad_type_ptr;

	Allocator *A_;
	size_t n_;
	int order_;
	double *primitiveConst_;
	typename AD::array_type a_, b_, c_;
	double *xBegin_;
	double *xEnd_;
	double *yBegin_;
	BoundaryCondition leftType_, rightType_;
	double leftValue_, rightValue_;
	InterpolatorType type_;
	DoubleArray dxScalar_;
	std::unique_ptr<FixedRegionAllocator> internal_allocator_;

	CubicInterpolator(double *xBegin, double *xEnd, double *yBegin, BoundaryCondition leftCondition,
			  double leftConditionValue, BoundaryCondition rightCondition, double rightConditionValue,
			  InterpolatorType it, int diff_order, Allocator *A)
	    : A_(A), n_(xEnd - xBegin), order_(diff_order), primitiveConst_(nullptr), a_(n_, order_, n_ - 1, A_),
	      b_(n_, order_, n_ - 1, A_), c_(n_, order_, n_ - 1, A_), xBegin_(xBegin), xEnd_(xEnd), yBegin_(yBegin),
	      leftType_(leftCondition), rightType_(rightCondition), leftValue_(leftConditionValue),
	      rightValue_(rightConditionValue), type_(it), dxScalar_(n_, order_, n_ - 1, A_)
	{
		primitiveConst_ = new (*A_) double[n_ - 1];
		std::fill_n(primitiveConst_, n_ - 1, 0.0);
		auto n_temps = n_ * 3 + (n_ - 1) * 3 + 2 + 1 + 2;
		internal_allocator_ = std::unique_ptr<FixedRegionAllocator>(new FixedRegionAllocator(
		    AD::alloc_size(n_, order_) * n_temps + 7 * sizeof(typename AD::array_type)));
		update_internal();
	}

	~CubicInterpolator()
	{
		// printf("CubicInterpolator destroyed\n");
		A_->deallocate(primitiveConst_);
	}

	size_t locate(double x) const { return locate_x(xBegin_, xEnd_, x); }
	double interpolate(double x)
	{
		size_t j = this->locate(x);
		double dx = x - this->xBegin_[j];
		double v = this->yBegin_[j] + dx * (AD::value(a_[j]) + dx * (AD::value(b_[j]) + dx * AD::value(c_[j])));
		return v;
	}
	void update() { update_internal(); }
	void update_internal();
	ad_type_ptr interpolate_with_sensitivities(double x, FixedRegionAllocator *A);
	bool coefficients(std::vector<double> &cvec)
	{
		cvec.clear();
		for (size_t i = 0; i < n_ - 1; i++) {
			cvec.push_back(AD::value(c_[i]));
		}
		return true;
	}
	double primitive(double);
	double derivative(double);
	double derivative2(double);
	InterpolatorType type() const { return type_; }
	size_t size() const { return n_; };
	double *xvalues() const { return xBegin_; }
	double *yvalues() const { return yBegin_; }
	int order() const { return order_; }
	void get_options(InterpolationOptions &options) const
	{
		options.differentiation_order = order_;
		options.cubic_left_condition_value = leftValue_;
		options.cubic_right_condition_value = rightValue_;
	}
};

template <typename T> void CubicInterpolator<T>::update_internal()
{
	FixedRegionAllocator &temp_allocator = *internal_allocator_;
	FixedRegionAllocatorGuard guard(&temp_allocator);

	typename AD::array_type tempvalues(n_, order_, 4, &temp_allocator);
	typename AD::array_type tmp_(n_, order_, n_, &temp_allocator);
	typename AD::array_type S_(n_, order_, n_ - 1, &temp_allocator);
	typename AD::array_type diagonal_(n_, order_, n_, &temp_allocator),
	    lowerDiagonal_(n_, order_, n_ - 1, &temp_allocator), upperDiagonal_(n_, order_, n_ - 1, &temp_allocator),
	    temp_(n_, order_, n_, &temp_allocator);

	auto &tempbuf1 = tempvalues[2]; // multiply
	auto &tempbuf2 = tempvalues[3]; // divide

	for (int i = 0; i < n_ - 1; ++i) {
		// dx_[i] = this->xBegin_[i+1] - this->xBegin_[i];
		dxScalar_[i] = this->xBegin_[i + 1] - this->xBegin_[i];
		assert(dxScalar_[i] > 0.0 && !close(dxScalar_[i], 0.0));
		// S_[i] = (this->yBegin_[i+1] - this->yBegin_[i])/dx_[i];
		AD::init(tempvalues[0], n_, order_, i, this->yBegin_[i]);
		AD::init(S_[i], n_, order_, i + 1, this->yBegin_[i + 1]);
		AD::add(S_[i], tempvalues[0], -1.0);
		AD::scalar_multiply(S_[i], 1.0 / dxScalar_[i]);
		assert(equal(AD::value(S_[i]), (this->yBegin_[i + 1] - this->yBegin_[i]) / dxScalar_[i]));
	}

	for (int i = 1; i < n_ - 1; ++i) {
		//  lowerDiagonal_[i - 1] = dx_[i];
		AD::init(lowerDiagonal_[i - 1], n_, order_, -1, dxScalar_[i]);
		//  diagonal_[i] = 2.0 * (dx_[i] + dx_[i - 1]);
		AD::init(diagonal_[i], n_, order_, -1, 2.0 * (dxScalar_[i] + dxScalar_[i - 1]));
		//  upperDiagonal_[i] = dx_[i - 1];
		AD::init(upperDiagonal_[i], n_, order_, -1, dxScalar_[i - 1]);
		//  tmp_[i] = 3.0 * (dx_[i] * S_[i - 1] + dx_[i - 1] * S_[i]);
		AD::assign(tempvalues[0], S_[i - 1]);
		AD::scalar_multiply(tempvalues[0], dxScalar_[i]);
		AD::assign(tmp_[i], S_[i]);
		AD::scalar_multiply(tmp_[i], dxScalar_[i - 1]);
		AD::add(tmp_[i], tempvalues[0], 1.0);
		AD::scalar_multiply(tmp_[i], 3.0);
		assert(equal(AD::value(tmp_[i]),
			     3.0 * (dxScalar_[i] * AD::value(S_[i - 1]) + dxScalar_[i - 1] * AD::value(S_[i]))));
	}

	// left boundary condition
	switch (leftType_) {
	case NotAKnot:
		// ignoring end condition value
		// diagonal_[0] = dx_[1] * (dx_[1] + dx_[0]);
		AD::init(diagonal_[0], n_, order_, -1, dxScalar_[1] * (dxScalar_[1] + dxScalar_[0]));
		// upperDiagonal_[0] = (dx_[0] + dx_[1]) * (dx_[0] +
		// dx_[1]);
		AD::init(upperDiagonal_[0], n_, order_, -1,
			 (dxScalar_[0] + dxScalar_[1]) * (dxScalar_[0] + dxScalar_[1]));
		// tmp_[0] = S_[0] * dx_[1] * (2.0 * dx_[1] + 3.0 *
		// dx_[0]) +
		//           S_[1] * dx_[0] * dx_[0];
		AD::assign(tempvalues[0], S_[1]);
		AD::scalar_multiply(tempvalues[0], dxScalar_[0] * dxScalar_[0]);
		AD::assign(tmp_[0], S_[0]);
		AD::scalar_multiply(tmp_[0], dxScalar_[1] * (2.0 * dxScalar_[1] + 3.0 * dxScalar_[0]));
		AD::add(tmp_[0], tempvalues[0], 1.0);
		assert(equal(AD::value(tmp_[0]),
			     AD::value(S_[0]) * dxScalar_[1] * (2.0 * dxScalar_[1] + 3.0 * dxScalar_[0]) +
				 AD::value(S_[1]) * dxScalar_[0] * dxScalar_[0]));
		break;
	case FirstDerivative:
		// diagonal_[0] = 1.0;
		AD::init(diagonal_[0], n_, order_, -1, 1.0);
		// upperDiagonal_[0] = 0.0;
		AD::init(upperDiagonal_[0], n_, order_, -1, 0.0);
		// tmp_[0] = AD::create(n_, order_, -1, leftValue_);
		AD::init(tmp_[0], n_, order_, -1, leftValue_);
		break;
	case SecondDerivative:
		// diagonal_[0] = 2.0;
		AD::init(diagonal_[0], n_, order_, -1, 2.0);
		// upperDiagonal_[0] = 1.0;
		AD::init(upperDiagonal_[0], n_, order_, -1, 1.0);
		// tmp_[0] =
		//  3.0 * S_[0] -
		//  AD::create(n_, order_, -1, leftValue_) * dx_[0]
		//  / 2.0;
		AD::assign(tmp_[0], S_[0]);
		AD::scalar_multiply(tmp_[0], 3.0);
		AD::init(tempvalues[0], n_, order_, -1, leftValue_);
		AD::add(tmp_[0], tempvalues[0], -dxScalar_[0] / 2.0);
		assert(equal(AD::value(tmp_[0]), 3.0 * AD::value(S_[0]) - leftValue_ * dxScalar_[0] / 2.0));
		break;
	default:
		assert(false);
	}

	// right boundary condition
	switch (rightType_) {
	case NotAKnot:
		// ignoring end condition value
		// lowerDiagonal_[n_ - 2] =
		//  -(dx_[n_ - 2] + dx_[n_ - 3]) * (dx_[n_ - 2] + dx_[n_
		//  - 3]);
		AD::init(lowerDiagonal_[n_ - 2], n_, order_, -1,
			 -(dxScalar_[n_ - 2] + dxScalar_[n_ - 3]) * (dxScalar_[n_ - 2] + dxScalar_[n_ - 3]));
		// diagonal_[n_ - 1] = -dx_[n_ - 3] * (dx_[n_ - 3] +
		// dx_[n_ - 2]);
		AD::init(diagonal_[n_ - 1], n_, order_, -1,
			 -dxScalar_[n_ - 3] * (dxScalar_[n_ - 3] + dxScalar_[n_ - 2]));
		// tmp_[n_ - 1] = MINUS * S_[n_ - 3] * dx_[n_ - 2] *
		// dx_[n_ - 2] -
		//  S_[n_ - 2] * dx_[n_ - 3] *
		//  (THREE * dx_[n_ - 2] + TWO * dx_[n_ - 3]);
		AD::assign(tempvalues[0], S_[n_ - 2]);
		AD::scalar_multiply(tempvalues[0],
				    dxScalar_[n_ - 3] * (3.0 * dxScalar_[n_ - 2] + 2.0 * dxScalar_[n_ - 3]));
		AD::assign(tmp_[n_ - 1], S_[n_ - 3]);
		AD::scalar_multiply(tmp_[n_ - 1], -(dxScalar_[n_ - 2] * dxScalar_[n_ - 2]));
		AD::add(tmp_[n_ - 1], tempvalues[0], -1.0);
		assert(equal(AD::value(tmp_[n_ - 1]), -AD::value(S_[n_ - 3]) * dxScalar_[n_ - 2] * dxScalar_[n_ - 2] -
							  AD::value(S_[n_ - 2]) * dxScalar_[n_ - 3] *
							      (3.0 * dxScalar_[n_ - 2] + 2.0 * dxScalar_[n_ - 3])));
		break;
	case FirstDerivative:
		// lowerDiagonal_[n_ - 2] = ZERO;
		AD::init(lowerDiagonal_[n_ - 2], n_, order_, -1, 0.0);
		// diagonal_[n_ - 1] = ONE;
		AD::init(diagonal_[n_ - 1], n_, order_, -1, 1.0);
		// tmp_[n_ - 1] = AD::create(n_, order_, -1,
		// rightValue_);
		AD::init(tmp_[n_ - 1], n_, order_, -1, rightValue_);
		break;
	case SecondDerivative:
		// lowerDiagonal_[n_ - 2] = ONE;
		AD::init(lowerDiagonal_[n_ - 2], n_, order_, -1, 1.0);
		// diagonal_[n_ - 1] = TWO;
		AD::init(diagonal_[n_ - 1], n_, order_, -1, 2.0);
		// tmp_[n_ - 1] = THREE * S_[n_ - 2] +
		//  AD::create(n_, order_, -1, rightValue_) *
		//  dx_[n_ - 2] / TWO;
		AD::assign(tmp_[n_ - 1], S_[n_ - 2]);
		AD::scalar_multiply(tmp_[n_ - 1], 3.0);
		AD::init(tempvalues[0], n_, order_, -1, rightValue_);
		AD::add(tmp_[n_ - 1], tempvalues[0], dxScalar_[n_ - 2] / 2.0);
		assert(equal(AD::value(tmp_[n_ - 1]),
			     3.0 * AD::value(S_[n_ - 2]) + rightValue_ * dxScalar_[n_ - 2] / 2.0));
		break;
	default:
		assert(false);
	}

	// solve the system
	// auto bet = diagonal_[0];
	typename AD::ad_type *bet = &diagonal_[0];
	// assert(!close(AD::value(bet), 0.0));
	assert(!close(AD::value(*bet), 0.0));
	// tmp_[0] = tmp_[0] / bet;
	AD::divide(tmp_[0], *bet, tempbuf1, tempbuf2);
	for (size_t j = 1; j <= n_ - 1; ++j) {
		// temp_[j] = upperDiagonal_[j - 1] / bet;
		AD::assign(temp_[j], upperDiagonal_[j - 1]);
		AD::divide(temp_[j], *bet, tempbuf1, tempbuf2);
		// bet = diagonal_[j] - lowerDiagonal_[j - 1] * temp_[j];
		AD::assign(tempvalues[0], lowerDiagonal_[j - 1]);
		AD::multiply(tempvalues[0], temp_[j], tempbuf1);
		AD::assign(tempvalues[1], diagonal_[j]);
		AD::add(tempvalues[1], tempvalues[0], -1.0);
		bet = &tempvalues[1];
		// assert(!close(AD::value(bet), 0.0));
		assert(!close(AD::value(*bet), 0.0));
		// tmp_[j] = (tmp_[j] - lowerDiagonal_[j - 1] * tmp_[j - 1]) /
		// bet;
		AD::assign(tempvalues[0], lowerDiagonal_[j - 1]);
		AD::multiply(tempvalues[0], tmp_[j - 1], tempbuf1);
		AD::add(tmp_[j], tempvalues[0], -1.0);
		AD::divide(tmp_[j], *bet, tempbuf1, tempbuf2);
	}
	// cannot be j>=0 with size_t j
	for (size_t j = n_ - 2; j > 0; --j) {
		// tmp_[j] -= temp_[j + 1] * tmp_[j + 1];
		AD::assign(tempvalues[0], temp_[j + 1]);
		AD::multiply(tempvalues[0], tmp_[j + 1], tempbuf1);
		AD::add(tmp_[j], tempvalues[0], -1.0);
	}
	// tmp_[0] -= temp_[1] * tmp_[1];
	AD::assign(tempvalues[0], temp_[1]);
	AD::multiply(tempvalues[0], tmp_[1], tempbuf1);
	AD::add(tmp_[0], tempvalues[0], -1.0);

	// cubic coefficients
	for (size_t i = 0; i < n_ - 1; ++i) {
		// a_[i] = tmp_[i];
		AD::assign(a_[i], tmp_[i]);
		// b_[i] = (3.0 * S_[i] - tmp_[i + 1] - 2.0 * tmp_[i]) / dx_[i];
		AD::assign(b_[i], S_[i]);
		AD::scalar_multiply(b_[i], 3.0);
		AD::add(b_[i], tmp_[i + 1], -1.0);
		AD::add(b_[i], tmp_[i], -2.0);
		AD::scalar_multiply(b_[i], 1.0 / dxScalar_[i]);
		assert(
		    equal(AD::value(b_[i]),
			  (3.0 * AD::value(S_[i]) - AD::value(tmp_[i + 1]) - 2.0 * AD::value(tmp_[i])) / dxScalar_[i]));
		// c_[i] = (tmp_[i + 1] + tmp_[i] - 2.0 * S_[i]) / (dx_[i] *
		// dx_[i]);
		AD::assign(c_[i], tmp_[i + 1]);
		AD::add(c_[i], tmp_[i], 1.0);
		AD::add(c_[i], S_[i], -2.0);
		AD::scalar_multiply(c_[i], 1.0 / (dxScalar_[i] * dxScalar_[i]));
		assert(equal(AD::value(c_[i]), (AD::value(tmp_[i + 1]) + AD::value(tmp_[i]) - 2.0 * AD::value(S_[i])) /
						   (dxScalar_[i] * dxScalar_[i])));
	}

	primitiveConst_[0] = 0.0;
	for (size_t i = 1; i < n_ - 1; ++i) {
		//    primitiveConst_[i] =
		//      primitiveConst_[i - 1] +
		//      AD::value(dx_[i - 1]) *
		//      (this->yBegin_[i - 1] +
		//        AD::value(dx_[i - 1]) *
		//        (AD::value(a_[i - 1]) / 2.0 +
		//          AD::value(dx_[i - 1]) *
		//          (AD::value(b_[i - 1]) / 3.0 +
		//            AD::value(dx_[i - 1]) *
		//            AD::value(c_[i - 1]) / 4.0)));
		primitiveConst_[i] =
		    primitiveConst_[i - 1] +
		    dxScalar_[i - 1] *
			(this->yBegin_[i - 1] +
			 dxScalar_[i - 1] * (AD::value(a_[i - 1]) / 2.0 +
					     dxScalar_[i - 1] * (AD::value(b_[i - 1]) / 3.0 +
								 dxScalar_[i - 1] * AD::value(c_[i - 1]) / 4.0)));
	}
}

template <typename T>
typename CubicInterpolator<T>::ad_type_ptr CubicInterpolator<T>::interpolate_with_sensitivities(double x,
												FixedRegionAllocator *A)
{
	size_t j = this->locate(x);
	//    adouble dx_ = adouble(n_, 2, -1, x) - adouble(n_, 2, -1,
	//    this->xBegin_[j]);
	//    return adouble(n_, 2, j, this->yBegin_[j]) + dx_ * (a_[j] + dx_ *
	//    (b_[j]
	//    + dx_ * c_[j]));
	double dx = x - this->xBegin_[j];
	// auto v = AD::create(n_, order_, j, this->yBegin_[j]) +
	//  dx_ * (a_[j] + dx_ * (b_[j] + dx_ * c_[j]));
	auto vptr = AD::allocate(A, n_, order_, j, this->yBegin_[j]);

	FixedRegionAllocatorGuard guard(A);
	typename AD::array_type tempvalues(n_, order_, 1, A);

	AD::assign(tempvalues[0], b_[j]);
	AD::add(tempvalues[0], c_[j], dx);
	AD::scalar_multiply(tempvalues[0], dx);
	AD::add(tempvalues[0], a_[j], 1.0);
	AD::scalar_multiply(tempvalues[0], dx);
	AD::add(*vptr, tempvalues[0], 1.0);
	return vptr;
}

template <typename T> double CubicInterpolator<T>::primitive(double x)
{
	size_t j = this->locate(x);
	double dx_ = x - this->xBegin_[j];
	return primitiveConst_[j] +
	       dx_ * (this->yBegin_[j] +
		      dx_ * (AD::value(a_[j]) / 2.0 + dx_ * (AD::value(b_[j]) / 3.0 + dx_ * AD::value(c_[j]) / 4.0)));
}

template <typename T> double CubicInterpolator<T>::derivative(double x)
{
	size_t j = this->locate(x);
	double dx_ = x - this->xBegin_[j];
	return AD::value(a_[j]) + (2.0 * AD::value(b_[j]) + 3.0 * AD::value(c_[j]) * dx_) * dx_;
}

template <typename T> double CubicInterpolator<T>::derivative2(double x)
{
	size_t j = this->locate(x);
	double dx_ = x - this->xBegin_[j];
	return 2.0 * AD::value(b_[j]) + 6.0 * AD::value(c_[j]) * dx_;
}

// Support functions for Monotone Convex interpolator

static inline int bound(int minimum, int v, int maximum)
{
	if (v < minimum)
		return minimum;
	else if (v > maximum)
		return maximum;
	else
		return v;
}

static inline double &bound(double &minimum, double &v, double &maximum)
{
	if (v < minimum)
		return minimum;
	else if (v > maximum)
		return maximum;
	else
		return v;
}

static int last_index(double *terms, int N, double term, int cache_i)
{
	int i = int(bound(0, cache_i, N));
	for (;;) {
		if (term >= terms[i]) {
			if (i == N) {
				if (term == terms[i]) {
					i = N - 1;
					break;
				} else {
					i = N;
					break;
				}
			} else {
				if (term >= terms[i + 1]) {
					i = i + 1;
				} else {
					break;
				}
			}
		} else {
			if (i == 0) {
				i = 0;
				break;
			} else {
				if (term >= terms[i - 1]) {
					i = i - 1;
					break;
				} else {
					i = i - 1;
				}
			}
		}
	}
	return i;
}

static inline bool is_zero(double v)
{
	if (v == 0.0)
		return true;
	// if (::fabs(v) < 1e-13) return true;
	return false;
}

static inline double square(double v) { return v * v; }

static inline double cube(double v) { return v * v * v; }

static inline redukti_adouble_t &bound(redukti_adouble_t &minimum, redukti_adouble_t &v, redukti_adouble_t &maximum)
{
	if (redukti_adouble_get_value(&v) < redukti_adouble_get_value(&minimum))
		return minimum;
	else if (redukti_adouble_get_value(&v) > redukti_adouble_get_value(&maximum))
		return maximum;
	else
		return v;
}

/**
 * Monotone Convex implementation using code from Hagan & West.
 * Inputs can be either zero rates or discrete forwards.
 * Interpolates zero rates.
 * Unlike other interpolators this one is designed specially for
 * zero rate curves, and is not usable in non-financial
 * context.
 */
template <typename T> class MonotoneConvexInterpolator
{
      public:
	typedef AutoDiffTraits<T> AD;
	typedef typename AD::ad_type ad_type;
	typedef typename AD::ptr_type ad_type_ptr;
	typedef typename AD::array_type ad_array_type;

	/**
	 * Inputs must be discrete forward rates or zero rates
	 */
	MonotoneConvexInterpolator(double *x, double *y, unsigned int size, bool inputsAreForwards,
				   bool allowNegativeForwards, InterpolatorType type, int order,
				   bool extrapolate = false, Allocator *alloc = &GlobalAllocator);
	~MonotoneConvexInterpolator() {}

	void update();
	double interpolate(double x)
	{
		FixedRegionAllocator *alloc = get_threadspecific_allocators()->tempspace_allocator;
		auto ptr = interpolate_with_sensitivities(x, alloc);
		return AD::value(ptr);
	}
	ad_type_ptr interpolate_with_sensitivities(double x, FixedRegionAllocator *A);
	ad_type_ptr forward_with_sensitivities(double t, FixedRegionAllocator *A);
	double forward(double t)
	{
		FixedRegionAllocator *alloc = get_threadspecific_allocators()->tempspace_allocator;
		auto ptr = forward_with_sensitivities(t, alloc);
		return AD::value(ptr);
	}
	InterpolatorType type() const { return type_; }
	size_t size() const { return size_; };
	double *xvalues() const { return terms_; }
	double *yvalues() const { return values_; }
	int order() const { return order_; }
	void get_options(InterpolationOptions &options) const
	{
		options.differentiation_order = order_;
		options.monotoneconvex_inputs_are_forwards = inputsAreForwards_;
		options.extrapolate = extrapolate_;
	}

      private:
	void estimateFI();

	bool inputsAreForwards_;
	/**
	 * Input values ( can be zero rates or discrete forwards )
	 */
	double *values_;
	/**
	 * Input knot points (first point is expected to be 0)
	 */
	double *terms_;
	/**
	 * If true negative forward rates will be allowed
	 */
	bool negativeForwardsAllowed_;
	/**
	 * Size of the arrays above
	 */
	unsigned int size_;
	/**
	 * Allocator
	 */
	Allocator *allocator_;
	/**
	 * N = last element in values/terms
	 */
	int N_;
	/**
	 * Should extrapolation be allowed?
	 */
	bool extrapolate_;
	// The type of monotone convex implementation
	InterpolatorType type_;
	// derivative order
	int order_;
	/**
	 * Intermediate value for interpolation
	 */
	ad_array_type interpolantAtNode_;
	/**
	 * Discrete forward rates
	 */
	ad_array_type fdiscrete_;
	/**
	 * instantaneous forward at knot point
	 */
	ad_array_type f_;

	std::unique_ptr<FixedRegionAllocator> internal_allocator_;

      private:
	MonotoneConvexInterpolator(const MonotoneConvexInterpolator &) = delete;
	MonotoneConvexInterpolator &operator=(const MonotoneConvexInterpolator &) = delete;
};

template <typename T>
MonotoneConvexInterpolator<T>::MonotoneConvexInterpolator(double *x, double *y, unsigned int size,
							  bool inputsAreForwards, bool allowNegativeForwards,
							  InterpolatorType type, int order, bool extrapolate,
							  Allocator *alloc)
    : inputsAreForwards_(inputsAreForwards), values_(y), terms_(x), negativeForwardsAllowed_(allowNegativeForwards),
      size_(size), allocator_(alloc), N_(size - 1), extrapolate_(extrapolate), type_(type), order_(order),
      interpolantAtNode_(size_, order_, size_, alloc), fdiscrete_(size_, order_, size_, alloc),
      f_(size_, order_, size_, alloc)
{
	assert(std::is_sorted(x, x + size));
	auto n_temps = 29;
	internal_allocator_ = std::unique_ptr<FixedRegionAllocator>(
	    new FixedRegionAllocator(AD::alloc_size(size_, order_) * n_temps + 7 * sizeof(typename AD::array_type)));
	estimateFI();
}

template <typename T> void MonotoneConvexInterpolator<T>::update() { estimateFI(); }

template <typename T> void MonotoneConvexInterpolator<T>::estimateFI()
{
	// Note that the use of internal allocator here is not thread safe
	// But this is okay as any thread modifying the interpolator should
	// endure that there is no race
	FixedRegionAllocator &temp_allocator = *internal_allocator_;
	FixedRegionAllocatorGuard guard(&temp_allocator);

	typename AD::array_type tempvalues(size_, order_, 2, &temp_allocator);

	terms_[0] = 0.0;
	values_[0] = values_[1];
	AD::init(fdiscrete_[0], size_, order_, 0, 0.0);
	AD::init(f_[0], size_, order_, 0, 0.0);
	AD::init(interpolantAtNode_[0], size_, order_, 0, 0.0);

	if (!inputsAreForwards_) {
		for (int j = 1; j <= N_; j++) {
			// fdiscrete_[j] =
			//  (terms_[j] * values_[j] - terms_[j - 1] * values_[j
			//  - 1])) /
			//    (terms_[j] - terms_[j - 1]);
			AD::init(tempvalues[0], size_, order_, j - 1, values_[j - 1]);
			AD::init(fdiscrete_[j], size_, order_, j, values_[j]);
			AD::scalar_multiply(fdiscrete_[j], terms_[j]);
			AD::add(fdiscrete_[j], tempvalues[0], -terms_[j - 1]);
			AD::scalar_multiply(fdiscrete_[j], 1.0 / (terms_[j] - terms_[j - 1]));
			assert(
			    equal(AD::value(fdiscrete_[j]), (terms_[j] * values_[j] - terms_[j - 1] * values_[j - 1]) /
								(terms_[j] - terms_[j - 1])));
			// interpolantAtNode_[j] =
			//  AD::create(size_, order_, j, values_[j]);
			AD::init(interpolantAtNode_[j], size_, order_, j, values_[j]);
		}
	} else {
		// T termrate = AD::create(size_, order_, -1, value_type(0.0));
		AD::init(tempvalues[0], size_, order_, -1, 0.0);
		for (int j = 1; j <= N_; j++) {
			// fdiscrete_[j] = AD::create(size_, order_, j,
			// values_[j]);
			AD::init(fdiscrete_[j], size_, order_, j, values_[j]);
			// termrate = termrate + fdiscrete_[j] * (terms_[j] -
			// terms_[j - 1]);
			AD::add(tempvalues[0], fdiscrete_[j], terms_[j] - terms_[j - 1]);
			// interpolantAtNode_[j] = termrate / terms_[j];
			AD::assign(interpolantAtNode_[j], tempvalues[0]);
			AD::scalar_multiply(interpolantAtNode_[j], 1.0 / terms_[j]);
		}
	}
	for (int j = 1; j <= N_ - 1; j++) {
		// f_[j] = ((terms_[j] - terms_[j - 1]) / (terms_[j + 1] -
		// terms_[j - 1])) *
		//  fdiscrete_[j + 1] +
		//  ((terms_[j + 1] - terms_[j]) / (terms_[j + 1] - terms_[j -
		//  1])) * fdiscrete_[j];
		AD::assign(f_[j], fdiscrete_[j + 1]);
		AD::scalar_multiply(f_[j], (terms_[j] - terms_[j - 1]) / (terms_[j + 1] - terms_[j - 1]));
		AD::add(f_[j], fdiscrete_[j], (terms_[j + 1] - terms_[j]) / (terms_[j + 1] - terms_[j - 1]));
		assert(equal(
		    AD::value(f_[j]),
		    ((terms_[j] - terms_[j - 1]) / (terms_[j + 1] - terms_[j - 1])) * AD::value(fdiscrete_[j + 1]) +
			((terms_[j + 1] - terms_[j]) / (terms_[j + 1] - terms_[j - 1])) * AD::value(fdiscrete_[j])));
	}
	// f_[0] = fdiscrete_[1] - 0.5 * (f_[1] - fdiscrete_[1]);
	AD::assign(tempvalues[0], f_[1]);
	AD::add(tempvalues[0], fdiscrete_[1], -1.0);
	AD::assign(f_[0], fdiscrete_[1]);
	AD::add(f_[0], tempvalues[0], -0.5);

	assert(equal(AD::value(f_[0]), AD::value(fdiscrete_[1]) - 0.5 * (AD::value(f_[1]) - AD::value(fdiscrete_[1]))));

	// f_[N_] = fdiscrete_[N_] - 0.5 * (f_[N_ - 1] - fdiscrete_[N_]);
	AD::assign(tempvalues[0], f_[N_ - 1]);
	AD::add(tempvalues[0], fdiscrete_[N_], -1.0);
	AD::assign(f_[N_], fdiscrete_[N_]);
	AD::add(f_[N_], tempvalues[0], -0.5);

	assert(equal(AD::value(f_[N_]),
		     AD::value(fdiscrete_[N_]) - 0.5 * (AD::value(f_[N_ - 1]) - AD::value(fdiscrete_[N_]))));

	if (!negativeForwardsAllowed_) {
		// f_[0] = bound(AD::create(size_, order_, 0, 0.0),
		//  f_[0], 2.0 * fdiscrete_[1]);
		AD::init(tempvalues[0], size_, order_, 0, 0.0);
		AD::assign(tempvalues[1], fdiscrete_[1]);
		AD::scalar_multiply(tempvalues[1], 2.0);
		auto &bounded = bound(tempvalues[0], f_[0], tempvalues[1]);
		if (&bounded != &f_[0]) {
			AD::assign(f_[0], bounded);
		}
		for (int j = 1; j <= N_ - 1; j++) {
			// f_[j] =
			//  bound(AD::create(size_, order_, j, 0.0), f_[j],
			//    2.0 *
			//    AD::min(fdiscrete_[j], fdiscrete_[j + 1]));
			auto &min_fdiscrete = AD::min(fdiscrete_[j], fdiscrete_[j + 1]);
			AD::assign(tempvalues[0], min_fdiscrete);
			AD::scalar_multiply(tempvalues[0], 2.0);
			AD::init(tempvalues[1], size_, order_, j, 0.0);
			auto &bounded_j = bound(tempvalues[1], f_[j], tempvalues[0]);
			if (&bounded_j != &f_[j]) {
				AD::assign(f_[j], bounded_j);
			}
		}
		// f_[N_] = bound(AD::create(size_, order_, N_, 0.0),
		//  f_[N_], 2.0 * fdiscrete_[N_]);
		AD::assign(tempvalues[0], fdiscrete_[N_]);
		AD::scalar_multiply(tempvalues[0], 2.0);
		AD::init(tempvalues[1], size_, order_, N_, 0.0);
		auto &bounded_N = bound(tempvalues[1], f_[N_], tempvalues[0]);
		if (&bounded_N != &f_[N_]) {
			AD::assign(f_[N_], bounded_N);
		}
	}
}

template <typename T>
typename MonotoneConvexInterpolator<T>::ad_type_ptr
MonotoneConvexInterpolator<T>::forward_with_sensitivities(double term, FixedRegionAllocator *A)
{
	if (term > terms_[N_])
		term = terms_[N_];

	// The result - this will be retained
	auto vptr = AD::allocate(A, size_, order_, -1, 0.0);

	// memory allocated below will be freed
	FixedRegionAllocatorGuard guard(A);
	StackRegionWithFallbackAllocator<128> membuf(A); // Optimize for the non derivative case
	typename AD::array_type tempvalues(size_, order_, 9, &membuf);

	auto &G = AD::reference(vptr);

	auto &g0 = tempvalues[0];
	auto &g1 = tempvalues[1];
	auto &eta = tempvalues[3];
	auto &tempbuf1 = tempvalues[7]; // for multiply
	auto &tempbuf2 = tempvalues[8]; // for divide

	if (term <= 0) {
		// return f_[0];
		AD::assign(G, f_[0]);
		return vptr;
	} else {
		int i = last_index(terms_, N_, term, 0);
		double x = (term - terms_[i]) / (terms_[i + 1] - terms_[i]);
		// T g0 = f_[i] - fdiscrete_[i + 1];
		AD::assign(g0, f_[i]);
		AD::add(g0, fdiscrete_[i + 1], -1.0);

		assert(equal(AD::value(g0), AD::value(f_[i]) - AD::value(fdiscrete_[i + 1])));

		// T g1 = f_[i + 1] - fdiscrete_[i + 1];
		AD::assign(g1, f_[i + 1]);
		AD::add(g1, fdiscrete_[i + 1], -1.0);

		assert(equal(AD::value(g1), AD::value(f_[i + 1]) - AD::value(fdiscrete_[i + 1])));

		// T G = AD::create(size_, order_, -1, value_type(0.0));
		// AD::init(G, size_, order_, -1, 0.0); // We don't need this
		if (is_zero(x)) {
			// G = g0;
			AD::assign(G, g0);
		} else if (is_zero(x - 1.0)) {
			// G = g1;
			AD::assign(G, g1);
		} else if ((AD::value(g0) < 0 && -0.5 * AD::value(g0) <= AD::value(g1) &&
			    AD::value(g1) <= -2.0 * AD::value(g0)) ||
			   (AD::value(g0) > 0 && -0.5 * AD::value(g0) >= AD::value(g1) &&
			    AD::value(g1) >= -2.0 * AD::value(g0))) {
			// G = g0 * (1.0 - 4.0 * x + 3.0 * square(x)) +
			//    g1 * (-2.0 * x + 3.0 * square(x));
			AD::assign(G, g0);
			AD::scalar_multiply(G, (1.0 - 4.0 * x + 3.0 * square(x)));
			AD::add(G, g1, (-2.0 * x + 3.0 * square(x)));

			assert(equal(AD::value(G), AD::value(g0) * (1.0 - 4.0 * x + 3.0 * square(x)) +
						       AD::value(g1) * (-2.0 * x + 3.0 * square(x))));

		} else if ((AD::value(g0) < 0 && AD::value(g1) > -2.0 * AD::value(g0)) ||
			   (AD::value(g0) > 0 && AD::value(g1) < -2.0 * AD::value(g0))) {
			// T eta = (g1 + 2.0 * g0) / (g1 - g0);
			AD::assign(eta, g1);
			AD::add(eta, g0, 2.0); // eta = (g1 + 2.0 * g0)
			AD::assign(tempvalues[4], g1);
			AD::add(tempvalues[4], g0, -1.0); // (g1 - g0)
			AD::divide(eta, tempvalues[4], tempbuf1,
				   tempbuf2); // eta = (g1 + 2.0 * g0) / (g1 - g0)

			assert(equal(AD::value(eta),
				     (AD::value(g1) + 2.0 * AD::value(g0)) / (AD::value(g1) - AD::value(g0))));

			if (x <= AD::value(eta)) {
				// G = g0;
				AD::assign(G, g0);
			} else {
				// G = g0 + (g1 - g0) * square((x - eta) / (1.0
				// - eta));
				AD::init(tempvalues[5], size_, order_, -1, x);
				AD::add(tempvalues[5], eta, -1.0); // x - eta
				AD::init(tempvalues[6], size_, order_, -1, 1.0);
				AD::add(tempvalues[6], eta, -1.0); // 1.0 - eta
				AD::divide(tempvalues[5], tempvalues[6], tempbuf1,
					   tempbuf2); // (x - eta) / (1.0 - eta)
				AD::multiply(tempvalues[5], tempvalues[5],
					     tempbuf1); // square((x - eta) /
				// (1.0 - eta))
				AD::multiply(tempvalues[4], tempvalues[5],
					     tempbuf1); // (g1 - g0) *
				// square((x - eta) /
				// (1.0 - eta))
				AD::assign(G, g0);
				AD::add(G, tempvalues[4], 1.0);

				assert(equal(AD::value(G), AD::value(g0) + (AD::value(g1) - AD::value(g0)) *
									       square((x - AD::value(eta)) /
										      (1.0 - AD::value(eta)))));
			}
		} else if ((AD::value(g0) > 0 && 0.0 > AD::value(g1) && AD::value(g1) > -0.5 * AD::value(g0)) ||
			   (AD::value(g0) < 0 && 0.0 < AD::value(g1) && AD::value(g1) < -0.5 * AD::value(g0))) {
			// T eta = 3.0 * g1 / (g1 - g0);
			AD::assign(eta, g1);
			AD::scalar_multiply(eta, 3.0); // eta = 3.0 * g1
			AD::assign(tempvalues[4], g1);
			AD::add(tempvalues[4], g0, -1.0); // (g1 - g0)
			AD::divide(eta, tempvalues[4], tempbuf1,
				   tempbuf2); // eta = 3.0 * g1 / (g1 - g0)

			assert(equal(AD::value(eta), 3.0 * AD::value(g1) / (AD::value(g1) - AD::value(g0))));

			if (x < AD::value(eta)) {
				// G = g1 + (g0 - g1) * square((eta - x) / eta);
				AD::assign(tempvalues[5], eta);
				AD::scalar_add(tempvalues[5], -x); // (eta - x)
				AD::divide(tempvalues[5], eta, tempbuf1,
					   tempbuf2); // (eta - x) / eta
				AD::multiply(tempvalues[5], tempvalues[5],
					     tempbuf1); // square((eta - x) / eta)
				AD::multiply(tempvalues[4], tempvalues[5],
					     tempbuf1); // (g0 - g1) *
				// square((eta - x) /
				// eta)
				AD::assign(G, g1);
				AD::add(G, tempvalues[4], 1.0);

				assert(equal(AD::value(G),
					     AD::value(g1) + (AD::value(g0) - AD::value(g1)) *
								 square((AD::value(eta) - x) / AD::value(eta))));

			} else {
				// G = g1;
				AD::assign(G, g1);
			}
		} else if (is_zero(AD::value(g0)) || is_zero(AD::value(g1))) {
			AD::init(G, size_, order_, -1, 0.0); // should -1 be i?
		} else {
			// T eta = g1 / (g1 + g0);
			AD::assign(tempvalues[4], g1);
			AD::add(tempvalues[4], g0, 1.0); // (g1 + g0)
			AD::assign(eta, g1);
			AD::divide(eta, tempvalues[4], tempbuf1,
				   tempbuf2); // eta = g1 / (g1 + g0)

			assert(equal(AD::value(eta), AD::value(g1) / (AD::value(g1) + AD::value(g0))));

			// T A = -1.0 * g0 * g1 / (g0 + g1);
			auto &A = tempvalues[5];
			AD::assign(A, g0);
			AD::scalar_multiply(A, -1.0);
			AD::multiply(A, g1, tempbuf1);
			AD::divide(A, tempvalues[4], tempbuf1, tempbuf2);

			assert(equal(AD::value(A),
				     -1.0 * AD::value(g0) * AD::value(g1) / (AD::value(g0) + AD::value(g1))));

			if (x <= AD::value(eta)) {
				// tempvalues 2 and 4 can be used
				// G = A + (g0 - A) * square((eta - x) / eta);
				AD::assign(tempvalues[2], eta);
				AD::scalar_add(tempvalues[2], -x); // (eta - x)
				AD::divide(tempvalues[2], eta, tempbuf1,
					   tempbuf2); // (eta - x) / eta
				AD::multiply(tempvalues[2], tempvalues[2],
					     tempbuf1); // square((eta - x) / eta)
				AD::assign(tempvalues[4], g0);
				AD::add(tempvalues[4], A, -1.0); // (g0 - A)
				AD::multiply(tempvalues[4], tempvalues[2],
					     tempbuf1); // (g0 - A) *
				// square((eta - x) /
				// eta)
				AD::assign(G, A);
				AD::add(G, tempvalues[4], 1.0);

				assert(equal(AD::value(G),
					     AD::value(A) + (AD::value(g0) - AD::value(A)) *
								square((AD::value(eta) - x) / AD::value(eta))));

			} else {
				// G = A + (g1 - A) * square((eta - x) / (1.0 -
				// eta));
				AD::init(tempvalues[2], size_, order_, -1, 1.0);
				AD::add(tempvalues[2], eta,
					-1.0); // (1.0 - eta)
				AD::assign(tempvalues[4], eta);
				AD::scalar_add(tempvalues[4], -x); // (eta - x)
				AD::divide(tempvalues[4], tempvalues[2], tempbuf1,
					   tempbuf2); // (eta - x) / (1.0 - eta)
				AD::multiply(tempvalues[4], tempvalues[4],
					     tempbuf1); // square((eta - x) /
				// (1.0 - eta))
				AD::assign(tempvalues[2], g1);
				AD::add(tempvalues[2], A, -1.0); // (g1 - A)
				AD::multiply(tempvalues[2], tempvalues[4],
					     tempbuf1); // (g1 - A) *
				// square((eta - x) /
				// (1.0 - eta))
				AD::assign(G, A);
				AD::add(G, tempvalues[2], 1.0);

				assert(equal(AD::value(G),
					     AD::value(A) + (AD::value(g1) - AD::value(A)) *
								square((AD::value(eta) - x) / (1.0 - AD::value(eta)))));
			}
		}
		// return G + fdiscrete_[i + 1];
		AD::add(G, fdiscrete_[i + 1], 1.0);
		return vptr;
	}
}

template <typename T>
typename MonotoneConvexInterpolator<T>::ad_type_ptr
MonotoneConvexInterpolator<T>::interpolate_with_sensitivities(double term, FixedRegionAllocator *A)
{
	if (term > terms_[N_]) {
		if (extrapolate_) {
			// return interpolate_with_sensitivities(terms_[N_]) *
			// terms_[N_] / term +
			//  forward(terms_[N_]) * (1.0 - terms_[N_] / term);
			auto vptr2 = interpolate_with_sensitivities(terms_[N_], A);
			AD::scalar_multiply(AD::reference(vptr2), terms_[N_] / term);
			auto fptr = forward_with_sensitivities(terms_[N_], A);
			AD::add(AD::reference(vptr2), AD::reference(fptr), (1.0 - terms_[N_] / term));
			return vptr2;
		} else
			term = terms_[N_];
	}

	// The interpolation result
	auto vptr = AD::allocate(A, size_, order_, -1, 0.0);

	// Memory allocated below will be freed
	FixedRegionAllocatorGuard guard(A);
	StackRegionWithFallbackAllocator<128> membuf(A); // Optimize for the non derivative case
	typename AD::array_type tempvalues(size_, order_, 10, &membuf);

	auto &retvalue = AD::reference(vptr);

	auto &g0 = tempvalues[0];
	auto &g1 = tempvalues[1];
	auto &G = tempvalues[2];
	auto &eta = tempvalues[3];
	auto &tempbuf1 = tempvalues[8]; // For multiply
	auto &tempbuf2 = tempvalues[9]; // For divide

	if (term <= 0.0) {
		// return AD::convert(f_[0]);
		AD::assign(retvalue, f_[0]);
		return vptr;
	} else {
		int i = last_index(terms_, N_, term, 0);

		// HW 'the x in (25)
		// HW x = (Term - Terms(i)) / (Terms(i + 1) - Terms(i))
		double x = (term - terms_[i]) / (terms_[i + 1] - terms_[i]);

		// HW g0 = f(i) - fdiscrete(i + 1)
		// T g0 = f_[i] - fdiscrete_[i + 1];
		AD::assign(g0, f_[i]);
		AD::add(g0, fdiscrete_[i + 1], -1.0);

		assert(equal(AD::value(g0), AD::value(f_[i]) - AD::value(fdiscrete_[i + 1])));

		// HW g1 = f(i + 1) - fdiscrete(i + 1)
		// T g1 = f_[i + 1] - fdiscrete_[i + 1];
		AD::assign(g1, f_[i + 1]);
		AD::add(g1, fdiscrete_[i + 1], -1.0);

		assert(equal(AD::value(g1), AD::value(f_[i + 1]) - AD::value(fdiscrete_[i + 1])));

		// T G;

		// HW If x = 0 Or x = 1 Then
		if (is_zero(x) || is_zero(x - 1.0)) {
			// HW G = 0
			// G = AD::create(size_, order_, -1, 0.0);
			AD::init(G, size_, order_, -1, 0.0);
		}

		// HW ElseIf g0 = 0 Or g1 = 0 Then
		else if (is_zero(AD::value(g0)) || is_zero(AD::value(g1))) {
			// HW G = 0
			// G = AD::create(size_, order_, -1, 0.0);
			AD::init(G, size_, order_, -1, 0.0);
		}

		// HW ElseIf (g0 < 0 And -0.5 * g0 <= g1 And g1 <= -2 * g0) _
		// HW Or(g0 > 0 And - 0.5 * g0 >= g1 And g1 >= -2 * g0) Then
		else if ((AD::value(g0) < 0.0 && -0.5 * AD::value(g0) <= AD::value(g1) &&
			  AD::value(g1) <= -2.0 * AD::value(g0)) ||
			 (AD::value(g0) > 0.0 && -0.5 * AD::value(g0) >= AD::value(g1) &&
			  AD::value(g1) >= -2.0 * AD::value(g0))) {
			// HW 'zone (i)
			// HW G = g0 * (x - 2 * x ^ 2 + x ^ 3) + g1 * (-x ^ 2 +
			// x ^ 3) G = g0 * (x - 2.0 * square(x) + cube(x)) +
			//    g1 * (-1.0 * square(x) + cube(x));
			AD::assign(G, g0);
			AD::scalar_multiply(G, (x - 2.0 * square(x) + cube(x)));
			AD::add(G, g1, (-1.0 * square(x) + cube(x)));

			assert(equal(AD::value(G), AD::value(g0) * (x - 2.0 * square(x) + cube(x)) +
						       AD::value(g1) * (-1.0 * square(x) + cube(x))));
		}

		// HW ElseIf (g0 < 0 And g1 > -2 * g0) Or (g0 > 0 And g1 < -2 *
		// g0) Then
		else if ((AD::value(g0) < 0.0 && AD::value(g1) > -2.0 * AD::value(g0)) ||
			 (AD::value(g0) > 0.0 && AD::value(g1) < -2.0 * AD::value(g0))) {
			// HW 'zone (ii)
			// HW '(29)
			// HW eta = (g1 + 2 * g0) / (g1 - g0)
			// T eta = (g1 + 2.0 * g0) / (g1 - g0);
			AD::assign(tempvalues[4], g1);
			AD::add(tempvalues[4], g0, -1.0);
			AD::assign(eta, g1);
			AD::add(eta, g0, 2.0);
			AD::divide(eta, tempvalues[4], tempbuf1, tempbuf2);

			assert(equal(AD::value(eta),
				     (AD::value(g1) + 2.0 * AD::value(g0)) / (AD::value(g1) - AD::value(g0))));

			// HW '(28)
			// HW If x <= eta Then
			if (x <= AD::value(eta)) {
				// HW G = g0 * x
				// G = g0 * x;
				AD::assign(G, g0);
				AD::scalar_multiply(G, x);
			} else {
				// HW G = g0 * x + (g1 - g0) * (x - eta) ^ 3 /
				// (1 - eta) ^ 2 / 3 G = g0 * x + (g1 - g0) *
				// cube(x - eta) / square(1.0 - eta) / 3.0;
				AD::init(tempvalues[5], size_, order_, -1, x);
				AD::add(tempvalues[5], eta, -1.0); // x - eta
				AD::pow(tempvalues[5], 3.0,
					tempbuf1); // cube(x - eta)
				AD::init(tempvalues[6], size_, order_, -1, 1.0);
				AD::add(tempvalues[6], eta, -1.0); // 1.0 - eta
				AD::multiply(tempvalues[6], tempvalues[6],
					     tempbuf1); // square(1.0 - eta)
				AD::multiply(tempvalues[4], tempvalues[5],
					     tempbuf1); // (g1 - g0) * cube(x - eta)
				AD::divide(tempvalues[4], tempvalues[6], tempbuf1,
					   tempbuf2); // (g1 - g0) * cube(x - eta) /
				// square(1.0 - eta)
				AD::assign(G, g0);
				AD::scalar_multiply(G, x);
				AD::add(G, tempvalues[4], 1.0 / 3.0);

				assert(equal(AD::value(G), AD::value(g0) * x + (AD::value(g1) - AD::value(g0)) *
										   cube(x - AD::value(eta)) /
										   square(1.0 - AD::value(eta)) / 3.0));
			}
		}

		// HW ElseIf (g0 > 0 And 0 > g1 And g1 > -0.5 * g0) _
		// HW Or(g0 < 0 And 0 < g1 And g1 < -0.5 * g0) Then
		else if ((AD::value(g0) > 0.0 && 0.0 > AD::value(g1) && AD::value(g1) > -0.5 * AD::value(g0)) ||
			 (AD::value(g0) < 0.0 && 0.0 < AD::value(g1) && AD::value(g1) < -0.5 * AD::value(g0))) {
			// HW 'zone (iii)
			// HW '(31)
			// HW eta = 3 * g1 / (g1 - g0)
			// T eta = 3.0 * g1 / (g1 - g0);
			AD::assign(tempvalues[4], g1);
			AD::add(tempvalues[4], g0, -1.0);
			AD::assign(eta, g1);
			AD::scalar_multiply(eta, 3.0);
			AD::divide(eta, tempvalues[4], tempbuf1, tempbuf2);

			assert(equal(AD::value(eta), 3.0 * AD::value(g1) / (AD::value(g1) - AD::value(g0))));

			// HW '(30)
			// HW If x < eta Then
			if (x < AD::value(eta)) {
				// HW G = g1 * x - 1 / 3 * (g0 - g1) * ((eta -
				// x) ^ 3 / eta ^ 2 - eta) G = g1 * x -
				//    1.0 / 3.0 * (g0 - g1) * (cube(eta - x) /
				//    square(eta) - eta);
				AD::assign(tempvalues[4], eta);
				AD::multiply(tempvalues[4], tempvalues[4],
					     tempbuf1); // square(eta)
				AD::assign(tempvalues[5], eta);
				AD::scalar_add(tempvalues[5], -x);
				AD::pow(tempvalues[5], 3.0,
					tempbuf1); // cube(eta - x)
				AD::divide(tempvalues[5], tempvalues[4], tempbuf1,
					   tempbuf2); // cube(eta - x) / square(eta)
				AD::add(tempvalues[5], eta,
					-1.0); // (cube(eta - x) / square(eta)
				// - eta)
				AD::assign(tempvalues[4], g0);
				AD::add(tempvalues[4], g1, -1.0); // (g0 - g1)
				AD::multiply(tempvalues[4], tempvalues[5],
					     tempbuf1); // (g0 - g1) *
				// (cube(eta - x) /
				// square(eta) - eta)
				AD::assign(G, g1);
				AD::scalar_multiply(G, x); // G = g1 * x
				AD::add(G, tempvalues[4],
					-1.0 / 3.0); // G = G - 1.0 / 3.0 * (g0 - g1)
				// * (cube(eta - x) /
				// square(eta) - eta);

				assert(equal(AD::value(G), AD::value(g1) * x -
							       1.0 / 3.0 * (AD::value(g0) - AD::value(g1)) *
								   (cube(AD::value(eta) - x) / square(AD::value(eta)) -
								    AD::value(eta))));

			} else {
				// HW G = (2 / 3 * g1 + 1 / 3 * g0) * eta + g1 *
				// (x - eta) G = (2.0 / 3.0 * g1 + 1.0 / 3.0 *
				// g0) * eta + g1 * (x - eta);
				AD::assign(G, g1);
				AD::scalar_multiply(G, 2.0 / 3.0);
				AD::add(G, g0, 1.0 / 3.0);
				AD::multiply(G, eta, tempbuf1);
				AD::init(tempvalues[4], size_, order_, -1, x);
				AD::add(tempvalues[4], eta, -1.0);
				AD::assign(tempvalues[5], g1);
				AD::multiply(tempvalues[5], tempvalues[4], tempbuf1);
				AD::add(G, tempvalues[5], 1.0);

				assert(equal(AD::value(G),
					     (2.0 / 3.0 * AD::value(g1) + 1.0 / 3.0 * AD::value(g0)) * AD::value(eta) +
						 AD::value(g1) * (x - AD::value(eta))));
			}
		} else {
			// HW 'zone (iv)
			// HW '(33)
			// HW eta = g1 / (g1 + g0)
			// T eta = g1 / (g1 + g0);
			AD::assign(tempvalues[4], g1);
			AD::add(tempvalues[4], g0, 1.0);
			AD::assign(eta, g1);
			AD::divide(eta, tempvalues[4], tempbuf1, tempbuf2);

			assert(equal(AD::value(eta), AD::value(g1) / (AD::value(g1) + AD::value(g0))));

			// HW '(34)
			// HW A = -g0 * g1 / (g0 + g1)
			// T A = -1.0 * g0 * g1 / (g0 + g1);
			auto &A = tempvalues[6];
			AD::assign(A, g0);
			AD::scalar_multiply(A, -1.0);
			AD::multiply(A, eta, tempbuf1);

			assert(equal(AD::value(A), -AD::value(g0) * AD::value(g1) / (AD::value(g0) + AD::value(g1))));

			// HW '(32)
			// HW If x <= eta Then
			if (x <= AD::value(eta)) {
				// HW G = A * x - 1 / 3 * (g0 - A) * ((eta - x)
				// ^ 3 / eta ^ 2 - eta) G = A * x - 1.0 / 3.0 *
				// (g0 - A) * (cube(eta - x) / square(eta) -
				// eta);
				AD::assign(tempvalues[4], eta);
				AD::multiply(tempvalues[4], tempvalues[4],
					     tempbuf1); // square(eta)
				AD::assign(tempvalues[5], eta);
				AD::scalar_add(tempvalues[5], -x);
				AD::pow(tempvalues[5], 3.0,
					tempbuf1); // cube(eta - x)
				AD::divide(tempvalues[5], tempvalues[4], tempbuf1,
					   tempbuf2); // cube(eta - x) / square(eta)
				AD::add(tempvalues[5], eta,
					-1.0); // (cube(eta - x) / square(eta)
				// - eta)
				AD::assign(tempvalues[4], g0);
				AD::add(tempvalues[4], A, -1.0); // (g0 - A)
				AD::multiply(tempvalues[4], tempvalues[5],
					     tempbuf1); // (g0 - A) * (cube(eta
				// - x) / square(eta) -
				// eta)
				AD::assign(G, A);
				AD::scalar_multiply(G, x);
				AD::add(G, tempvalues[4], -1.0 / 3.0);

				assert(equal(AD::value(G),
					     AD::value(A) * x - 1.0 / 3.0 * (AD::value(g0) - AD::value(A)) *
								    (cube(AD::value(eta) - x) / square(AD::value(eta)) -
								     AD::value(eta))));

			} else {
				// HW G = (2 / 3 * A + 1 / 3 * g0) * eta + A *
				// (x - eta) + (g1 - A) / 3
				// * (x - eta) ^ 3 / (1 - eta) ^ 2
				// G = (2.0 / 3.0 * A + 1.0 / 3.0 * g0) * eta +
				// A * (x - eta) +
				//    (g1 - A) / 3.0 * cube(x - eta) /
				//    square(1.0 - eta);
				AD::init(tempvalues[4], size_, order_, -1, x);
				AD::add(tempvalues[4], eta, -1.0); // x - eta
				AD::assign(tempvalues[5], tempvalues[4]);
				AD::pow(tempvalues[5], 3.0,
					tempbuf1); // cube(x - eta)
				AD::init(tempvalues[7], size_, order_, -1, 1.0);
				AD::add(tempvalues[7], eta, -1.0); // 1.0 - eta
				AD::multiply(tempvalues[7], tempvalues[7],
					     tempbuf1); // square(1.0 - eta)
				AD::divide(tempvalues[5], tempvalues[7], tempbuf1,
					   tempbuf2); // cube(x - eta) /
				// square(1.0 - eta)
				AD::assign(tempvalues[7], g1);
				AD::add(tempvalues[7], A, -1.0); // (g1 - A)
				AD::scalar_multiply(tempvalues[7],
						    1.0 / 3.0); // (g1 - A) / 3.0
				AD::multiply(tempvalues[7], tempvalues[5],
					     tempbuf1); // (g1 - A) / 3.0 * cube(x
				// - eta) / square(1.0 -
				// eta)
				AD::assign(tempvalues[5], A);
				AD::multiply(tempvalues[5], tempvalues[4],
					     tempbuf1); // A * (x - eta)
				AD::assign(G, A);
				AD::scalar_multiply(G, 2.0 / 3.0);
				AD::add(G, g0, 1.0 / 3.0);
				AD::multiply(G, eta,
					     tempbuf1); // (2.0 / 3.0 * A + 1.0
				// / 3.0 * g0) * eta
				AD::add(G, tempvalues[5],
					1.0); // + A * (x - eta)
				AD::add(G, tempvalues[7],
					1.0); // + (g1 - A) / 3.0 * cube(x -
				// eta) / square(1.0 - eta)

				assert(equal(AD::value(G),
					     (2.0 / 3.0 * AD::value(A) + 1.0 / 3.0 * AD::value(g0)) * AD::value(eta) +
						 AD::value(A) * (x - AD::value(eta)) +
						 (AD::value(g1) - AD::value(A)) / 3.0 * cube(x - AD::value(eta)) /
						     square(1.0 - AD::value(eta))));
			}
		}

		// HW '(12)
		// HW Interpolant = 1 / Term * (Terms(i) * dInterpolantatNode(i)
		// + (Term - Terms(i)) * fdiscrete(i + 1) + (Terms(i + 1) -
		// Terms(i)) * G) auto retvalue = 1.0 / term * (terms_[i] *
		// interpolantAtNode_[i] +
		//                              (term - terms_[i]) *
		//                              fdiscrete_[i + 1] + (terms_[i +
		//                              1] - terms_[i]) * G);
		AD::assign(retvalue, interpolantAtNode_[i]);
		AD::scalar_multiply(retvalue, terms_[i]);
		AD::add(retvalue, fdiscrete_[i + 1], term - terms_[i]);
		AD::add(retvalue, G, terms_[i + 1] - terms_[i]);
		AD::scalar_multiply(retvalue, 1.0 / term);

		assert(equal(AD::value(retvalue), 1.0 / term *
						      (terms_[i] * AD::value(interpolantAtNode_[i]) +
						       (term - terms_[i]) * AD::value(fdiscrete_[i + 1]) +
						       (terms_[i + 1] - terms_[i]) * AD::value(G))));

		return vptr;
	}
}

class CubicInterpolatorFacade : public AbstractInterpolator
{
	Allocator *A_;
	std::unique_ptr<CubicInterpolator<redukti_adouble_t>, Deleter<CubicInterpolator<redukti_adouble_t>>> ad_;
	std::unique_ptr<CubicInterpolator<double>, Deleter<CubicInterpolator<double>>> non_ad_;

      public:
	CubicInterpolatorFacade(double *xBegin, double *xEnd, double *yBegin, BoundaryCondition leftCondition,
				double leftConditionValue, BoundaryCondition rightCondition, double rightConditionValue,
				InterpolatorType it, int order, Allocator *A)
	    : A_(A)
	{
		non_ad_ = std::unique_ptr<CubicInterpolator<double>, Deleter<CubicInterpolator<double>>>(
		    new (*A) CubicInterpolator<double>(xBegin, xEnd, yBegin, leftCondition, leftConditionValue,
						       rightCondition, rightConditionValue, it, order, A_),
		    Deleter<CubicInterpolator<double>>(A));
		if (order > 0)
			ad_ = std::unique_ptr<CubicInterpolator<redukti_adouble_t>,
					      Deleter<CubicInterpolator<redukti_adouble_t>>>(
			    new (*A) CubicInterpolator<redukti_adouble_t>(xBegin, xEnd, yBegin, leftCondition,
									  leftConditionValue, rightCondition,
									  rightConditionValue, it, order, A_),
			    Deleter<CubicInterpolator<redukti_adouble_t>>(A));
	}
	virtual ~CubicInterpolatorFacade()
	{
		// printf("CubicInterpolatorFacade destroyed\n");
	}

	/**
	 * Interpolate at x
	 */
	virtual double interpolate(double x) override final { return non_ad_->interpolate(x); }

	virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
	interpolate_with_sensitivities(double x, FixedRegionAllocator *A) override final
	{
		if (!ad_) {
			return std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>(nullptr);
		}
		return ad_->interpolate_with_sensitivities(x, A);
	}
	/**
	 * Return the first derivative of the interpolated
	 * function at x. Note - not implemented for Monotone
	 * Convex
	 */
	virtual double derivative(double x) override final { return non_ad_->derivative(x); }
	/**
	 * Return the second derivative of the interpolated
	 * function at x. Note - not implemented for Monotone
	 * Convex
	 */
	virtual double derivative2(double x) override final { return non_ad_->derivative2(x); }
	/**
	 * Note - not implemented for Log* interpolators
	 * or for Monotone Convex.
	 */
	virtual double primitive(double x) override final { return non_ad_->primitive(x); }
	/**
	 * If underlying values have changed, this
	 * method can be called to reinitialise the
	 * interpolator.
	 */
	virtual void update() override final
	{
		non_ad_->update();
		if (ad_)
			ad_->update();
	}

	virtual bool coefficients(std::vector<double> &cvec) override final { return non_ad_->coefficients(cvec); }

	/**
	 * Return the interpolator type
	 */
	virtual InterpolatorType type() const override final { return non_ad_->type(); }

	virtual size_t size() const override final { return non_ad_->size(); };
	virtual double *xvalues() const override final { return non_ad_->xvalues(); }
	virtual double *yvalues() const override final { return non_ad_->yvalues(); }
	virtual int order() const override final { return non_ad_->order(); }
	virtual void get_options(InterpolationOptions &options) const override final { non_ad_->get_options(options); }
	virtual std::unique_ptr<Interpolator, Deleter<Interpolator>>
	get_interpolator(double *xa, double *ya, unsigned int size, Allocator *A,
			 const InterpolationOptions &options) const override final
	{
		return make_interpolator(type(), xa, ya, size, A, options);
	}
};

// convenience classes
class CubicSplineNotAKnotImpl : public CubicInterpolatorFacade
{
      public:
	/*! \pre the \f$ x \f$ values must be sorted. */
	CubicSplineNotAKnotImpl(double *xBegin, double *xEnd, double *yBegin, Allocator *A, int order = 0,
				double k1 = 0.0, double k2 = 0.0)
	    : CubicInterpolatorFacade(xBegin, xEnd, yBegin, NotAKnot, k1, NotAKnot, k2,
				      InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT, order, A)
	{
	}
};

// convenience classes
class CubicSplineClampedImpl : public CubicInterpolatorFacade
{
      public:
	/*! \pre the \f$ x \f$ values must be sorted. */
	CubicSplineClampedImpl(double *xBegin, double *xEnd, double *yBegin, Allocator *A, int order = 0,
			       double k1 = 0.0, double k2 = 0.0)
	    : CubicInterpolatorFacade(xBegin, xEnd, yBegin, FirstDerivative, k1, FirstDerivative, k2,
				      InterpolatorType::CUBIC_SPLINE_CLAMPED, order, A)
	{
	}
};

// convenience classes
class CubicNaturalSplineImpl : public CubicInterpolatorFacade
{
      public:
	/*! \pre the \f$ x \f$ values must be sorted. */
	CubicNaturalSplineImpl(double *xBegin, double *xEnd, double *yBegin, Allocator *A, int order = 0,
			       double k1 = 0.0, double k2 = 0.0)
	    : CubicInterpolatorFacade(xBegin, xEnd, yBegin, SecondDerivative, k1, SecondDerivative, k2,
				      InterpolatorType::CUBIC_SPLINE_NATURAL, order, A)
	{
	}
};

// For performance reasons we don't want to use the
// differentiable implementation all the time; the facade
// keeps two implementations - a plain double one and
// an autodiff one, so that it can efficiently handle
// simple interpolation requests
class MonotoneConvexInterpolatorFacade : public AbstractInterpolator
{
	// Differentiable interpolator
	std::unique_ptr<MonotoneConvexInterpolator<redukti_adouble_t>,
			Deleter<MonotoneConvexInterpolator<redukti_adouble_t>>>
	    ad_;

	// Plain interpolator
	std::unique_ptr<MonotoneConvexInterpolator<double>, Deleter<MonotoneConvexInterpolator<double>>> non_ad_;
	/**
	 * Allocator
	 */
	Allocator *A_;

      public:
	MonotoneConvexInterpolatorFacade(double *x, double *y, unsigned int size, bool inputsAreForwards,
					 bool allowNegativeForwards, InterpolatorType type, int order, Allocator *alloc,
					 bool extrapolate = false)
	    : A_(alloc)
	{
		non_ad_ =
		    std::unique_ptr<MonotoneConvexInterpolator<double>, Deleter<MonotoneConvexInterpolator<double>>>(
			new (*A_) MonotoneConvexInterpolator<double>(
			    x, y, size, inputsAreForwards, allowNegativeForwards, type, order, extrapolate, alloc),
			Deleter<MonotoneConvexInterpolator<double>>(alloc));
		if (order > 0)
			ad_ = std::unique_ptr<MonotoneConvexInterpolator<redukti_adouble_t>,
					      Deleter<MonotoneConvexInterpolator<redukti_adouble_t>>>(
			    new (*A_) MonotoneConvexInterpolator<redukti_adouble_t>(
				x, y, size, inputsAreForwards, allowNegativeForwards, type, order, extrapolate, alloc),
			    Deleter<MonotoneConvexInterpolator<redukti_adouble_t>>(alloc));
	}

	virtual ~MonotoneConvexInterpolatorFacade()
	{
		// printf("MonotoneConvexInterpolatorFacade destroyed\n");
	}

	/**
	 * Interpolate at x
	 */
	virtual double interpolate(double x) override final
	{
		// return interpolate_with_sensitivities(x).getValue();
		return non_ad_->interpolate(x);
	}

	virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
	interpolate_with_sensitivities(double x, FixedRegionAllocator *A) override final
	{
		if (!ad_) {
			return std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>(nullptr);
		}
		return ad_->interpolate_with_sensitivities(x, A);
	}
	virtual void update() override final
	{
		non_ad_->update();
		if (ad_)
			ad_->update();
	}
	/**
	 * Return the interpolator type
	 */
	virtual InterpolatorType type() const override final { return non_ad_->type(); }
	virtual size_t size() const override final { return non_ad_->size(); };
	virtual double *xvalues() const override final { return non_ad_->xvalues(); }
	virtual double *yvalues() const override final { return non_ad_->yvalues(); }
	virtual int order() const override final { return non_ad_->order(); }
	virtual void get_options(InterpolationOptions &options) const override final { non_ad_->get_options(options); }
	virtual std::unique_ptr<Interpolator, Deleter<Interpolator>>
	get_interpolator(double *xa, double *ya, unsigned int size, Allocator *A,
			 const InterpolationOptions &options) const override final
	{
		return make_interpolator(type(), xa, ya, size, A, options);
	}

	virtual double forward(double t) override final { return non_ad_->forward(t); }

	// virtual adouble forward_with_sensitivities(double t) {
	//  if (!ad_) {
	//    ad_ = new MonotoneConvexInterpolatorT<adouble>(
	//      terms_, values_, size_, inputsAreForwards_,
	//      negativeForwardsAllowed_, type_, order_, allocator_,
	//      extrapolate_);
	//  }
	//  return ad_->forward_with_sensitivities(t);
	//}
};

std::unique_ptr<Interpolator, Deleter<Interpolator>> make_interpolator(InterpolatorType type, double *xa, double *ya,
								       unsigned int size, Allocator *A,
								       const InterpolationOptions &options)
{
	double k1 = 0.0, k2 = 0.0;
	int order = options.differentiation_order;
	bool extrapolate = options.extrapolate;
	switch (type) {
	case InterpolatorType::LINEAR:
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) LinearInterpolator(xa, ya, size, order, A), Deleter<Interpolator>(A));
	case InterpolatorType::LOG_LINEAR:
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) LogInterpolator(xa, ya, size, InterpolatorType::LINEAR, InterpolatorType::LOG_LINEAR, A,
					     options),
		    Deleter<Interpolator>(A));
	case InterpolatorType::FLAT_LEFT:
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) BackwardFlatInterpolation(xa, ya, size, order), Deleter<Interpolator>(A));
	case InterpolatorType::FLAT_RIGHT:
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) ForwardFlatInterpolation(xa, ya, size, order), Deleter<Interpolator>(A));
	case InterpolatorType::MONOTONE_CONVEX:
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A)
			MonotoneConvexInterpolatorFacade(xa, ya, size, options.monotoneconvex_inputs_are_forwards, true,
							 InterpolatorType::MONOTONE_CONVEX, order, A, extrapolate),
		    Deleter<Interpolator>(A));
	case InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT:
		k1 = options.cubic_left_condition_value;
		k2 = options.cubic_right_condition_value;
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) CubicSplineNotAKnotImpl(xa, xa + size, ya, A, order, k1, k2), Deleter<Interpolator>(A));
	case InterpolatorType::CUBIC_SPLINE_CLAMPED:
		k1 = options.cubic_left_condition_value;
		k2 = options.cubic_right_condition_value;
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) CubicSplineClampedImpl(xa, xa + size, ya, A, order, k1, k2), Deleter<Interpolator>(A));
	case InterpolatorType::CUBIC_SPLINE_NATURAL:
		k1 = options.cubic_left_condition_value;
		k2 = options.cubic_right_condition_value;
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) CubicNaturalSplineImpl(xa, xa + size, ya, A, order, k1, k2), Deleter<Interpolator>(A));
	case InterpolatorType::LOG_CUBIC_SPLINE_NATURAL:
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(
		    new (*A) LogInterpolator(xa, ya, size, InterpolatorType::CUBIC_SPLINE_NATURAL,
					     InterpolatorType::LOG_CUBIC_SPLINE_NATURAL, A, options),
		    Deleter<Interpolator>(A));
	default:
		assert(false);
		return std::unique_ptr<Interpolator, Deleter<Interpolator>>(nullptr);
	}
}

//////////////////////////////// TESTS //////////////////////////////

static int unequal(double got, double expected)
{
	if (std::fabs(got - expected) > 1e-6) {
		fprintf(stderr, "mismatch of values: expected %f, got %f\n", expected, got);
		return 1;
	}
	return 0;
}

namespace
{

std::vector<double> xRange(double start, double finish, size_t points)
{
	std::vector<double> x(points);
	double dx = (finish - start) / (points - 1);
	for (size_t i = 0; i < points - 1; i++)
		x[i] = start + i * dx;
	x[points - 1] = finish;
	return x;
}

std::vector<double> gaussian(const std::vector<double> &x)
{
	std::vector<double> y(x.size());
	for (size_t i = 0; i < x.size(); i++)
		y[i] = std::exp(-x[i] * x[i]);
	return y;
}

std::vector<double> parabolic(const std::vector<double> &x)
{
	std::vector<double> y(x.size());
	for (size_t i = 0; i < x.size(); i++)
		y[i] = -x[i] * x[i];
	return y;
}

template <class I, class J> int checkValues(const char *type, Interpolator &cubic, I xBegin, I xEnd, J yBegin)
{
	int failure_count = 0;
	double tolerance = 2.0e-15;
	while (xBegin != xEnd) {
		double interpolated = cubic.interpolate(*xBegin);
		if (std::fabs(interpolated - *yBegin) > tolerance) {
			std::cerr << " interpolation failed at x = " << *xBegin
				  << "\n    interpolated value: " << interpolated
				  << "\n    expected value:     " << *yBegin
				  << "\n    error:              " << std::fabs(interpolated - *yBegin);
		}
		++xBegin;
		++yBegin;
	}
	return failure_count;
}

int check1stDerivativeValue(const char *type, Interpolator &cubic, double x, double value)
{
	int failure_count = 0;
	double tolerance = 1.0e-14;
	double interpolated = cubic.derivative(x);
	double error = std::fabs(interpolated - value);
	if (error > tolerance) {
		std::cerr << " interpolation first derivative failure\n"
			  << "at x = " << x << "\n    interpolated value: " << interpolated
			  << "\n    expected value:     " << value << "\n    error:              " << error;
		failure_count++;
	}
	return failure_count;
}

int check2ndDerivativeValue(const char *type, Interpolator &cubic, double x, double value)
{
	int failure_count = 0;
	double tolerance = 1.0e-13;
	double interpolated = cubic.derivative2(x);
	double error = std::fabs(interpolated - value);
	if (error > tolerance) {
		std::cerr << " interpolation second derivative failure\n"
			  << "at x = " << x << "\n    interpolated value: " << interpolated
			  << "\n    expected value:     " << value << "\n    error:              " << error;
		failure_count++;
	}
	return failure_count;
}

int checkNotAKnotCondition(const char *type, Interpolator &cubic)
{
	int failure_count = 0;
	double tolerance = 1.0e-14;
	std::vector<double> c;
	cubic.coefficients(c);
	if (std::fabs(c[0] - c[1]) > tolerance) {
		std::cerr << " interpolation failure"
			  << "\n    cubic coefficient of the first"
			  << " polinomial is " << c[0] << "\n    cubic coefficient of the second"
			  << " polinomial is " << c[1];
		failure_count++;
	}
	size_t n = c.size();
	if (std::fabs(c[n - 2] - c[n - 1]) > tolerance) {
		std::cerr << " interpolation failure"
			  << "\n    cubic coefficient of the 2nd to last"
			  << " polinomial is " << c[n - 2] << "\n    cubic coefficient of the last"
			  << " polinomial is " << c[n - 1];
		failure_count++;
	}
	return failure_count;
}

int checkSymmetry(const char *type, Interpolator &cubic, double xMin)
{
	int failure_count = 0;
	double tolerance = 1.0e-15;
	for (double x = xMin; x < 0.0; x += 0.1) {
		double y1 = cubic.interpolate(x), y2 = cubic.interpolate(-x);
		if (std::fabs(y1 - y2) > tolerance) {
			std::cerr << " interpolation not symmetric"
				  << "\n    x = " << x << "\n    g(x)  = " << y1 << "\n    g(-x) = " << y2
				  << "\n    error:  " << std::fabs(y1 - y2);
			failure_count++;
		}
	}
	return failure_count;
}

template <class F> class errorFunction : public std::unary_function<double, double>
{
      public:
	errorFunction(F &f) : f_(f) {}
	double operator()(double x) const
	{
		double temp = f_(x) - std::exp(-x * x);
		return temp * temp;
	}

      private:
	F &f_;
};

template <class F> errorFunction<F> make_error_function(F &f) { return errorFunction<F>(f); }

#if 0
double multif(double s, double t, double u, double v, double w) {
  return std::sqrt(s * std::sinh(std::log(t)) +
                   std::exp(std::sin(u) * std::sin(3 * v)) +
                   std::sinh(std::log(v * w)));
}

double epanechnikovKernel(double u) {
  if (std::fabs(u) <= 1) {
    return (3.0 / 4.0) * (1 - u * u);
  } else {
    return 0.0;
  }
}
#endif
} // namespace

class Integrator
{
      public:
	Integrator(double absoluteAccuracy, size_t maxEvaluations);
	virtual ~Integrator() {}

	double operator()(const std::function<double(double)> &f, double a, double b) const;

	//! \name Modifiers
	//@{
	void setAbsoluteAccuracy(double);
	void setMaxEvaluations(size_t);
	//@}

	//! \name Inspectors
	//@{
	double absoluteAccuracy() const;
	size_t maxEvaluations() const;
	//@}

	double absoluteError() const;

	size_t numberOfEvaluations() const;

	virtual bool integrationSuccess() const;

      protected:
	virtual double integrate(const std::function<double(double)> &f, double a, double b) const = 0;
	void setAbsoluteError(double error) const;
	void setNumberOfEvaluations(size_t evaluations) const;
	void increaseNumberOfEvaluations(size_t increase) const;

      private:
	double absoluteAccuracy_;
	mutable double absoluteError_;
	size_t maxEvaluations_;
	mutable double evaluations_;
};

Integrator::Integrator(double absoluteAccuracy, size_t maxEvaluations)
    : absoluteAccuracy_(absoluteAccuracy), maxEvaluations_(maxEvaluations)
{
	// QL_REQUIRE(absoluteAccuracy > QL_EPSILON,
	//           std::scientific << "required tolerance (" <<
	//           absoluteAccuracy << ") not allowed. It must be > " <<
	//           QL_EPSILON);
}

void Integrator::setAbsoluteAccuracy(double accuracy) { absoluteAccuracy_ = accuracy; }

void Integrator::setMaxEvaluations(size_t maxEvaluations) { maxEvaluations_ = maxEvaluations; }

double Integrator::absoluteAccuracy() const { return absoluteAccuracy_; }

size_t Integrator::maxEvaluations() const { return maxEvaluations_; }

double Integrator::absoluteError() const { return absoluteError_; }

void Integrator::setAbsoluteError(double error) const { absoluteError_ = error; }

size_t Integrator::numberOfEvaluations() const { return evaluations_; }

void Integrator::setNumberOfEvaluations(size_t evaluations) const { evaluations_ = evaluations; }

void Integrator::increaseNumberOfEvaluations(size_t increase) const { evaluations_ += increase; }

bool Integrator::integrationSuccess() const
{
	return evaluations_ <= maxEvaluations_ && absoluteError_ <= absoluteAccuracy_;
}

double Integrator::operator()(const std::function<double(double)> &f, double a, double b) const
{
	evaluations_ = 0;
	if (a == b)
		return 0.0;
	if (b > a)
		return integrate(f, a, b);
	else
		return -integrate(f, b, a);
}

//! Integral of a one-dimensional function
/*! Given a target accuracy \f$ \epsilon \f$, the integral of
a function \f$ f \f$ between \f$ a \f$ and \f$ b \f$ is
calculated by means of the trapezoid formula
\f[
\int_{a}^{b} f \mathrm{d}x =
\frac{1}{2} f(x_{0}) + f(x_{1}) + f(x_{2}) + \dots
+ f(x_{N-1}) + \frac{1}{2} f(x_{N})
\f]
where \f$ x_0 = a \f$, \f$ x_N = b \f$, and
\f$ x_i = a+i \Delta x \f$ with
\f$ \Delta x = (b-a)/N \f$. The number \f$ N \f$ of intervals
is repeatedly increased until the target accuracy is reached.

\test the correctness of the result is tested by checking it
against known good values.
*/
template <class IntegrationPolicy> class TrapezoidIntegral : public Integrator
{
      public:
	TrapezoidIntegral(double accuracy, size_t maxIterations) : Integrator(accuracy, maxIterations) {}

      protected:
	double integrate(const std::function<double(double)> &f, double a, double b) const
	{
		// start from the coarsest trapezoid...
		size_t N = 1;
		double I = (f(a) + f(b)) * (b - a) / 2.0, newI;
		// ...and refine it
		size_t i = 1;
		do {
			newI = IntegrationPolicy::integrate(f, a, b, I, N);
			N *= IntegrationPolicy::nbEvalutions();
			// good enough? Also, don't run away immediately
			if (std::fabs(I - newI) <= absoluteAccuracy() && i > 5)
				// ok, exit
				return newI;
			// oh well. Another step.
			I = newI;
			i++;
		} while (i < maxEvaluations());
		assert(false);
		return 0.0;
	}
};

// Integration policies
struct Default {
	inline static double integrate(const std::function<double(double)> &f, double a, double b, double I, size_t N)
	{
		double sum = 0.0;
		double dx = (b - a) / N;
		double x = a + dx / 2.0;
		for (size_t i = 0; i < N; x += dx, ++i)
			sum += f(x);
		return (I + dx * sum) / 2.0;
	}
	inline static size_t nbEvalutions() { return 2; }
};

struct MidPoint {
	inline static double integrate(const std::function<double(double)> &f, double a, double b, double I, size_t N)
	{
		double sum = 0.0;
		double dx = (b - a) / N;
		double x = a + dx / 6.0;
		double D = 2.0 * dx / 3.0;
		for (size_t i = 0; i < N; x += dx, ++i)
			sum += f(x) + f(x + D);
		return (I + dx * sum) / 3.0;
	}
	inline static size_t nbEvalutions() { return 3; }
};

//! Integral of a one-dimensional function
/*! \test the correctness of the result is tested by checking it
against known good values.
*/
class SimpsonIntegral : public TrapezoidIntegral<Default>
{
      public:
	SimpsonIntegral(double accuracy, size_t maxIterations) : TrapezoidIntegral<Default>(accuracy, maxIterations) {}

      protected:
	double integrate(const std::function<double(double)> &f, double a, double b) const
	{
		// start from the coarsest trapezoid...
		size_t N = 1;
		double I = (f(a) + f(b)) * (b - a) / 2.0, newI;
		double adjI = I, newAdjI;
		// ...and refine it
		size_t i = 1;
		do {
			newI = Default::integrate(f, a, b, I, N);
			N *= 2;
			newAdjI = (4.0 * newI - I) / 3.0;
			// good enough? Also, don't run away immediately
			if (std::fabs(adjI - newAdjI) <= absoluteAccuracy() && i > 5)
				// ok, exit
				return newAdjI;
			// oh well. Another step.
			I = newI;
			adjI = newAdjI;
			i++;
		} while (i < maxEvaluations());
		assert(false);
		return 0.0;
	}
};

/* See J. M. Hyman, "Accurate monotonicity preserving cubic interpolation"
SIAM J. of Scientific and Statistical Computing, v. 4, 1983, pp. 645-654.
http://math.lanl.gov/~mac/papers/numerics/H83.pdf
*/
int testSplineErrorOnGaussianValues()
{
	int failure_count = 0;

	// fprintf(stderr, "Testing spline approximation on Gaussian data
	// sets...\n");

	size_t points[] = {5, 9, 17, 33};

	// complete spline data from the original 1983 Hyman paper
	double tabulatedErrors[] = {3.5e-2, 2.0e-3, 4.0e-5, 1.8e-6};
	double toleranceOnTabErr[] = {0.1e-2, 0.1e-3, 0.1e-5, 0.1e-6};

	// (complete) MC spline data from the original 1983 Hyman paper
	// NB: with the improved Hyman filter from the Dougherty, Edelman, and
	//     Hyman 1989 paper the n=17 nonmonotonicity is not filtered anymore
	//     so the error agrees with the non MC method.
	// double tabulatedMCErrors[] = {1.7e-2, 2.0e-3, 4.0e-5, 1.8e-6};
	// double toleranceOnTabMCErr[] = {0.1e-2, 0.1e-3, 0.1e-5, 0.1e-6};

	SimpsonIntegral integral(1e-12, 10000);
	std::vector<double> x, y;

	// still unexplained scale factor needed to obtain the numerical
	// results from the paper
	double scaleFactor = 1.9;

	size_t len = std::end(points) - std::begin(points);
	for (size_t i = 0; i < len; i++) {
		size_t n = points[i];
		std::vector<double> x = xRange(-1.7, 1.9, n);
		std::vector<double> y = gaussian(x);

		// Not-a-knot
		auto f = make_interpolator(InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT, x.data(), y.data(),
					   x.end() - x.begin(), &GlobalAllocator);
		double result = std::sqrt(integral(make_error_function(*f), -1.7, 1.9));
		result /= scaleFactor;
		if (std::fabs(result - tabulatedErrors[i]) > toleranceOnTabErr[i]) {
			std::cerr << "Not-a-knot spline interpolation "
				  << "\n    sample points:      " << n << "\n    norm of difference: " << result
				  << "\n    it should be:       " << tabulatedErrors[i];
			failure_count++;
		}
// MC not-a-knot
#if 0
    f = make_interpolator(InterpolatorType::MonotonicCubicSplineNotAKnot, x.data(), y.data(), x.end() - x.begin(), &GlobalAllocator);
    result = std::sqrt(integral(make_error_function(*f), -1.7, 1.9));
    result /= scaleFactor;
    if (std::fabs(result - tabulatedMCErrors[i]) > toleranceOnTabMCErr[i]) {
      std::cerr << "MC Not-a-knot spline interpolation "
        << "\n    sample points:      " << n
        << "\n    norm of difference: " << result
        << "\n    it should be:       "
        << tabulatedMCErrors[i];
      failure_count++;
    }
#endif
	}
	return failure_count;
}

/* See J. M. Hyman, "Accurate monotonicity preserving cubic interpolation"
SIAM J. of Scientific and Statistical Computing, v. 4, 1983, pp. 645-654.
http://math.lanl.gov/~mac/papers/numerics/H83.pdf
*/
int testSplineOnRPN15AValues()
{
	int failure_count = 0;
	// fprintf(stderr, "Testing spline interpolation on RPN15A data
	// set...\n");

	double RPN15A_x[] = {7.99, 8.09, 8.19, 8.7, 9.2, 10.0, 12.0, 15.0, 20.0};
	double RPN15A_y[] = {0.0, 2.76429e-5, 4.37498e-5, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994};

	double interpolated;

	// Natural spline
	auto f = make_interpolator(InterpolatorType::CUBIC_SPLINE_NATURAL, RPN15A_x, RPN15A_y,
				   std::end(RPN15A_x) - std::begin(RPN15A_x), &GlobalAllocator);
	failure_count +=
	    checkValues("Natural spline", *f, std::begin(RPN15A_x), std::end(RPN15A_x), std::begin(RPN15A_y));
	failure_count += check2ndDerivativeValue("Natural spline", *f, *std::begin(RPN15A_x), 0.0);
	failure_count += check2ndDerivativeValue("Natural spline", *f, *(std::end(RPN15A_x) - 1), 0.0);
	// poor performance
	double x_bad = 11.0;
	interpolated = f->interpolate(x_bad);
	if (interpolated < 1.0) {
		std::cerr << "Natural spline interpolation "
			  << "poor performance unverified\n"
			  << "at x = " << x_bad << "\ninterpolated value: " << interpolated << "\nexpected value > 1.0";
	}

	// Clamped spline
	f = make_interpolator(InterpolatorType::CUBIC_SPLINE_CLAMPED, RPN15A_x, RPN15A_y,
			      std::end(RPN15A_x) - std::begin(RPN15A_x), &GlobalAllocator);
	failure_count +=
	    checkValues("Clamped spline", *f, std::begin(RPN15A_x), std::end(RPN15A_x), std::begin(RPN15A_y));
	failure_count += check1stDerivativeValue("Clamped spline", *f, *std::begin(RPN15A_x), 0.0);
	failure_count += check1stDerivativeValue("Clamped spline", *f, *(std::end(RPN15A_x) - 1), 0.0);
	// poor performance
	interpolated = f->interpolate(x_bad);
	if (interpolated < 1.0) {
		std::cerr << "Clamped spline interpolation "
			  << "poor performance unverified\n"
			  << "at x = " << x_bad << "\ninterpolated value: " << interpolated << "\nexpected value > 1.0";
		failure_count++;
	}

	// Not-a-knot spline
	f = make_interpolator(InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT, RPN15A_x, RPN15A_y,
			      std::end(RPN15A_x) - std::begin(RPN15A_x), &GlobalAllocator);
	failure_count +=
	    checkValues("Not-a-knot spline", *f, std::begin(RPN15A_x), std::end(RPN15A_x), std::begin(RPN15A_y));
	failure_count += checkNotAKnotCondition("Not-a-knot spline", *f);
	// poor performance
	interpolated = f->interpolate(x_bad);
	if (interpolated < 1.0) {
		std::cerr << "Not-a-knot spline interpolation "
			  << "poor performance unverified\n"
			  << "at x = " << x_bad << "\ninterpolated value: " << interpolated << "\nexpected value > 1.0";
		failure_count++;
	}

#if 0
  // MC natural spline values
  f = make_interpolator(InterpolatorType::MonotonicCubicSplineNatural, RPN15A_x, RPN15A_y, std::end(RPN15A_x) - std::begin(RPN15A_x), &GlobalAllocator);
  failure_count += checkValues("MC natural spline", *f,
    std::begin(RPN15A_x), std::end(RPN15A_x), std::begin(RPN15A_y));
  // good performance
  interpolated = f->interpolate(x_bad);
  if (interpolated>1.0) {
    std::cerr << "MC natural spline interpolation "
      << "good performance unverified\n"
      << "at x = " << x_bad
      << "\ninterpolated value: " << interpolated
      << "\nexpected value < 1.0";
    failure_count++;
  }


  // MC clamped spline values
  f = make_interpolator(InterpolatorType::MonotonicCubicSplineClamped, RPN15A_x, RPN15A_y, std::end(RPN15A_x) - std::begin(RPN15A_x), &GlobalAllocator);
  failure_count += checkValues("MC clamped spline", *f,
    std::begin(RPN15A_x), std::end(RPN15A_x), std::begin(RPN15A_y));
  failure_count += check1stDerivativeValue("MC clamped spline", *f,
    *std::begin(RPN15A_x), 0.0);
  failure_count += check1stDerivativeValue("MC clamped spline", *f,
    *(std::end(RPN15A_x) - 1), 0.0);
  // good performance
  interpolated = f->interpolate(x_bad);
  if (interpolated>1.0) {
    std::cerr << "MC clamped spline interpolation "
      << "good performance unverified\n"
      << "at x = " << x_bad
      << "\ninterpolated value: " << interpolated
      << "\nexpected value < 1.0";
    failure_count++;
  }


  // MC not-a-knot spline values
  f = make_interpolator(InterpolatorType::MonotonicCubicSplineNotAKnot, RPN15A_x, RPN15A_y, std::end(RPN15A_x) - std::begin(RPN15A_x), &GlobalAllocator);
  failure_count += checkValues("MC not-a-knot spline", *f,
    std::begin(RPN15A_x), std::end(RPN15A_x), std::begin(RPN15A_y));
  // good performance
  interpolated = f->interpolate(x_bad);
  if (interpolated>1.0) {
    std::cerr << "MC clamped spline interpolation "
      << "good performance unverified\n"
      << "at x = " << x_bad
      << "\ninterpolated value: " << interpolated
      << "\nexpected value < 1.0";
    failure_count++;
  }
#endif
	return failure_count;
}

/* Blossey, Frigyik, Farnum "A Note On CubicSpline Splines"
Applied Linear Algebra and Numerical Analysis AMATH 352 Lecture Notes
http://www.amath.washington.edu/courses/352-winter-2002/spline_note.pdf
*/
int testSplineOnGenericValues()
{
	int failure_count = 0;
	// fprintf(stderr, "Testing spline interpolation on generic
	// values...\n");

	double generic_x[] = {0.0, 1.0, 3.0, 4.0};
	double generic_y[] = {0.0, 0.0, 2.0, 2.0};
	double generic_natural_y2[] = {0.0, 1.5, -1.5, 0.0};

	double interpolated, error;
	size_t i, n = std::end(generic_x) - std::begin(generic_x);
	std::vector<double> x35(3);

	InterpolationOptions data;
	data.cubic_left_condition_value = generic_natural_y2[0];
	data.cubic_right_condition_value = generic_natural_y2[n - 1];
	// Natural spline
	auto f = make_interpolator(InterpolatorType::CUBIC_SPLINE_NATURAL, generic_x, generic_y,
				   std::end(generic_x) - std::begin(generic_x), &GlobalAllocator, data);
	failure_count +=
	    checkValues("Natural spline", *f, std::begin(generic_x), std::end(generic_x), std::begin(generic_y));
	// cached second derivative
	for (i = 0; i < n; i++) {
		interpolated = f->derivative2(generic_x[i]);
		error = interpolated - generic_natural_y2[i];
		if (std::fabs(error) > 3e-16) {
			std::cerr << "Natural spline interpolation "
				  << "second derivative failed at x=" << generic_x[i]
				  << "\ninterpolated value: " << interpolated
				  << "\nexpected value:     " << generic_natural_y2[i]
				  << "\nerror:              " << error;
		}
	}
	x35[1] = f->interpolate(3.5);

	// Clamped spline
	double y1a = 0.0, y1b = 0.0;
	data.cubic_left_condition_value = y1a;
	data.cubic_right_condition_value = y1b;
	f = make_interpolator(InterpolatorType::CUBIC_SPLINE_CLAMPED, generic_x, generic_y,
			      std::end(generic_x) - std::begin(generic_x), &GlobalAllocator, data);
	failure_count +=
	    checkValues("Clamped spline", *f, std::begin(generic_x), std::end(generic_x), std::begin(generic_y));
	failure_count += check1stDerivativeValue("Clamped spline", *f, *std::begin(generic_x), 0.0);
	failure_count += check1stDerivativeValue("Clamped spline", *f, *(std::end(generic_x) - 1), 0.0);
	x35[0] = f->interpolate(3.5);

	// Not-a-knot spline
	f = make_interpolator(InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT, generic_x, generic_y,
			      std::end(generic_x) - std::begin(generic_x), &GlobalAllocator);
	failure_count +=
	    checkValues("Not-a-knot spline", *f, std::begin(generic_x), std::end(generic_x), std::begin(generic_y));
	failure_count += checkNotAKnotCondition("Not-a-knot spline", *f);

	x35[2] = f->interpolate(3.5);

	if (x35[0] > x35[1] || x35[1] > x35[2]) {
		std::cerr << "Spline interpolation failure"
			  << "\nat x = " << 3.5 << "\nclamped spline    " << x35[0] << "\nnatural spline    " << x35[1]
			  << "\nnot-a-knot spline " << x35[2] << "\nvalues should be in increasing order";
		failure_count++;
	}
	return failure_count;
}

int testSimmetricEndConditions()
{
	int failure_count = 0;

	// fprintf(stderr,
	//        "Testing symmetry of spline interpolation
	//        end-conditions...\n");

	size_t n = 9;

	std::vector<double> x, y;
	x = xRange(-1.8, 1.8, n);
	y = gaussian(x);

	// Not-a-knot spline
	auto f = make_interpolator(InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT, x.data(), y.data(), x.size(),
				   &GlobalAllocator);
	failure_count += checkValues("Not-a-knot spline", *f, x.begin(), x.end(), y.begin());
	failure_count += checkNotAKnotCondition("Not-a-knot spline", *f);
	failure_count += checkSymmetry("Not-a-knot spline", *f, x[0]);

#if 0
  // MC not-a-knot spline
  f = make_interpolator(InterpolatorType::MonotonicCubicSplineNotAKnot, x.data(), y.data(), x.size(), &GlobalAllocator);
  failure_count += checkValues("MC not-a-knot spline", *f,
    x.begin(), x.end(), y.begin());
  failure_count += checkSymmetry("MC not-a-knot spline", *f, x[0]);
#endif
	return failure_count;
}

int testDerivativeEndConditions()
{
	int failure_count = 0;

	// fprintf(stderr,
	//        "Testing derivative end-conditions for spline
	//        interpolation...\n");

	size_t n = 4;

	std::vector<double> x, y;
	x = xRange(-2.0, 2.0, n);
	y = parabolic(x);

	// Not-a-knot spline
	auto f = make_interpolator(InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT, x.data(), y.data(), x.size(),
				   &GlobalAllocator);
	failure_count += checkValues("Not-a-knot spline", *f, x.begin(), x.end(), y.begin());
	failure_count += check1stDerivativeValue("Not-a-knot spline", *f, x[0], 4.0);
	failure_count += check1stDerivativeValue("Not-a-knot spline", *f, x[n - 1], -4.0);
	failure_count += check2ndDerivativeValue("Not-a-knot spline", *f, x[0], -2.0);
	failure_count += check2ndDerivativeValue("Not-a-knot spline", *f, x[n - 1], -2.0);

	// Clamped spline
	InterpolationOptions data;
	data.cubic_left_condition_value = 4.0;
	data.cubic_right_condition_value = -4.0;
	f = make_interpolator(InterpolatorType::CUBIC_SPLINE_CLAMPED, x.data(), y.data(), x.size(), &GlobalAllocator,
			      data);
	failure_count += checkValues("Clamped spline", *f, x.begin(), x.end(), y.begin());
	failure_count += check1stDerivativeValue("Clamped spline", *f, x[0], 4.0);
	failure_count += check1stDerivativeValue("Clamped spline", *f, x[n - 1], -4.0);
	failure_count += check2ndDerivativeValue("Clamped spline", *f, x[0], -2.0);
	failure_count += check2ndDerivativeValue("Clamped spline", *f, x[n - 1], -2.0);

	// SecondDerivative spline
	// f = CubicInterpolation(x.begin(), x.end(), y.begin(),
	//	CubicInterpolation::Spline, false,
	//	CubicInterpolation::SecondDerivative, -2.0,
	//	CubicInterpolation::SecondDerivative, -2.0);
	// f.update();
	// checkValues("SecondDerivative spline", f,
	//	x.begin(), x.end(), y.begin());
	// check1stDerivativeValue("SecondDerivative spline", f,
	//	x[0], 4.0);
	// check1stDerivativeValue("SecondDerivative spline", f,
	//	x[n - 1], -4.0);
	// check2ndDerivativeValue("SecondDerivative spline", f,
	//	x[0], -2.0);
	// check2ndDerivativeValue("SecondDerivative spline", f,
	//	x[n - 1], -2.0);

#if 0
  // MC Not-a-knot spline
  f = make_interpolator(InterpolatorType::MonotonicCubicSplineNotAKnot, x.data(), y.data(), x.size(), &GlobalAllocator);
  failure_count += checkValues("MC Not-a-knot spline", *f,
    x.begin(), x.end(), y.begin());
  failure_count += check1stDerivativeValue("MC Not-a-knot spline", *f,
    x[0], 4.0);
  failure_count += check1stDerivativeValue("MC Not-a-knot spline", *f,
    x[n - 1], -4.0);
  failure_count += check2ndDerivativeValue("MC Not-a-knot spline", *f,
    x[0], -2.0);
  failure_count += check2ndDerivativeValue("MC Not-a-knot spline", *f,
    x[n - 1], -2.0);

  // MC Clamped spline
  data.cubic_left_condition_value = 4.0;
  data.cubic_right_condition_value = -4.0;
  f = make_interpolator(InterpolatorType::MonotonicCubicSplineClamped, x.data(), y.data(), x.size(), &GlobalAllocator, data);
  failure_count += checkValues("MC Clamped spline", *f,
    x.begin(), x.end(), y.begin());
  failure_count += check1stDerivativeValue("MC Clamped spline", *f,
    x[0], 4.0);
  failure_count += check1stDerivativeValue("MC Clamped spline", *f,
    x[n - 1], -4.0);
  failure_count += check2ndDerivativeValue("MC Clamped spline", *f,
    x[0], -2.0);
  failure_count += check2ndDerivativeValue("MC Clamped spline", *f,
    x[n - 1], -2.0);
#endif

	// MC SecondDerivative spline
	// f = CubicInterpolation(x.begin(), x.end(), y.begin(),
	//	CubicInterpolation::Spline, true,
	//	CubicInterpolation::SecondDerivative, -2.0,
	//	CubicInterpolation::SecondDerivative, -2.0);
	// f.update();
	// checkValues("MC SecondDerivative spline", f,
	//	x.begin(), x.end(), y.begin());
	// check1stDerivativeValue("MC SecondDerivative spline", f,
	//	x[0], 4.0);
	// check1stDerivativeValue("MC SecondDerivative spline", f,
	//	x[n - 1], -4.0);
	// check2ndDerivativeValue("SecondDerivative spline", f,
	//	x[0], -2.0);
	// check2ndDerivativeValue("MC SecondDerivative spline", f,
	//	x[n - 1], -2.0);
	return failure_count;
}

int test_simple_cubicspline()
{
	int failure_count = 0;

	double xa[] = {0, 2, 4, 6};
	double ya[] = {5, 8, 9, 11.5};

	double expected[] = {5, 6.7375, 8, 8.5375, 9, 10.05, 11.5};

	CubicInterpolator<redukti_adouble_t> interp(
	    std::begin(xa), std::end(xa), ya, BoundaryCondition::SecondDerivative, 0.0,
	    BoundaryCondition::SecondDerivative, 0.0, InterpolatorType::CUBIC_SPLINE_NATURAL, 2, &GlobalAllocator);

	for (int x = 0; x < 7; x++) {
		double y = interp.interpolate(x);
		if (::fabs(y - expected[x]) > 0.0001) {
			printf("Expected %f got %f\n", expected[x], y);
			failure_count++;
		}

#if 0
    StackRegionWithFallbackAllocator<1024> buf;
    auto ad = interp.interpolate_with_sensitivities(x, &buf);
    
    for (int i = 0; i < 4; i++) {
      double sens = interp.sensitivity_to(x, i);
      double sens2 = redukti_adouble_get_derivative1(ad.get(), i);
      if (::fabs(sens - sens2) > 1e-10) {
        printf("Numerical sens to pillar %d is %f\n", i, sens);
        printf("Auto diff sens to pillar %d is %f\n", i, sens2);
        failure_count++;
      }
      for (int j = 0; j < 4; j++) {
        double sens3 = interp.sensitivity_to_to(x, i, j);
        double sens4 = redukti_adouble_get_derivative2(ad.get(), i, j);
        if (::fabs(sens3 - sens4) > 1e-6) {
          printf("Numerical sens to pillar %d,%d is %f\n", i, j, sens3);
          printf("Auto diff sens to pillar %d,%d is %f\n", i, j, sens4);
          failure_count++;
        }
      }
    }
#endif
	}
	return failure_count;
}

static int test_interp(double *x_data, double *y_data, size_t n, InterpolatorType T, double *test_x, double *test_y,
		       double *test_dy, double *test_di, size_t test_n)
{
	int status = 0;
	size_t i;

	auto interp = make_interpolator(T, x_data, y_data, n, &GlobalAllocator);

	for (i = 0; i < test_n; i++) {
		double x = test_x[i];
		double diff_y, diff_deriv, diff_integ;

		double y = interp->interpolate(x);
		double deriv = interp->derivative(x);
		double integ = interp->primitive(x);

		diff_y = y - test_y[i];
		diff_deriv = deriv - test_dy[i];
		diff_integ = integ - test_di[i];
		if (fabs(diff_y) > 1.e-10 || fabs(diff_deriv) > 1.0e-10 || fabs(diff_integ) > 1.0e-10) {
			fprintf(stderr, "Failed to match: %d %.10f %.10f %.10f\n", (int)i, diff_y, diff_deriv,
				diff_integ);
			status++;
		}
	}

	return status;
}

static int test_linear_sensitivities()
{
#if 0
	double xa[] = {0, 2, 4, 6};
	double ya[] = {5, 8, 9, 11.5};

	double xValues[] = {0, 0, 1, 1, 2, 3, 4, 5, 6};
	int knotPoint[] = {0, 1, 0, 1, 1, 2, 2, 2, 2};
	double expected[] = {1, 0, .5, .5, 1., .5, 1, .5, 0};
	double interpolated[] = {5., 5., 6.5, 6.5, 8., 8.5, 9., 10.25, 11.5};

	InterpolationOptions data;
	data.differentiation_order = 1;
	auto interp = make_interpolator(InterpolatorType::LINEAR, xa, ya, 4,
					&GlobalAllocator, data);
	int failure_count = 0;
	auto n = std::end(xValues) - std::begin(xValues);
	for (int i = 0; i < n; i++) {
		double y = interp->interpolate(xValues[i]);
		double result =
		    interp->sensitivity_to(xValues[i], knotPoint[i]);
		FixedRegionAllocator *allocator =
		    get_threadspecific_allocators()->tempspace_allocator;
		FixedRegionAllocatorGuard guard(allocator);
		auto sens = interp->interpolate_with_sensitivities(xValues[i],
								   allocator);
		if (!equal(y, interpolated[i])) {
			printf("Interpolated value at %f: Expected %f got %f\n",
			       xValues[i], y, interpolated[i]);
			failure_count++;
		}
		if (!equal(y, redukti_adouble_get_value(sens.get()))) {
			printf("Interpolated value at %f: Expected %f got %f\n",
			       xValues[i], y,
			       redukti_adouble_get_value(sens.get()));
			failure_count++;
		}
		if (!equal(result, expected[i])) {
			printf(
			    "Numeric sensitivity to node %d at %f: Expected %f "
			    "got %f\n",
			    knotPoint[i], xValues[i], expected[i], result);
			failure_count++;
		}
		double result2 =
		    redukti_adouble_get_derivative1(sens.get(), knotPoint[i]);
		if (!equal(result, result2)) {
			printf(
			    "Autodiff sensitivity to node %d at %f: Expected "
			    "%f got %f\n",
			    knotPoint[i], xValues[i], result, result2);
			failure_count++;
		}
		// sens.dump(stderr, "node sensitivity");
		// printf("\n");
	}
	return failure_count;
#endif
	return 0;
}

static int test_linear(void)
{
	int s;

	double data_x[4] = {0.0, 1.0, 2.0, 3.0};
	double data_y[4] = {0.0, 1.0, 2.0, 3.0};
	double test_x[6] = {0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
	double test_y[6] = {0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
	double test_dy[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double test_iy[6] = {0.0, 0.125, 0.5, 9.0 / 8.0, 25.0 / 8.0, 9.0 / 2.0};

	s = test_interp(data_x, data_y, std::end(data_x) - std::begin(data_x), InterpolatorType::LINEAR, test_x, test_y,
			test_dy, test_iy, std::end(test_x) - std::begin(test_x));
	return s;
}

static int test_cspline(void)
{
	int s;

	double data_x[3] = {0.0, 1.0, 2.0};
	double data_y[3] = {0.0, 1.0, 2.0};
	double test_x[4] = {0.0, 0.5, 1.0, 2.0};
	double test_y[4] = {0.0, 0.5, 1.0, 2.0};
	double test_dy[4] = {1.0, 1.0, 1.0, 1.0};
	double test_iy[4] = {0.0, 0.125, 0.5, 2.0};

	s = test_interp(data_x, data_y, std::end(data_x) - std::begin(data_x), InterpolatorType::CUBIC_SPLINE_NATURAL,
			test_x, test_y, test_dy, test_iy, std::end(test_x) - std::begin(test_x));
	return s;
}

static int test_cspline2(void)
{
	int s;

	double data_x[6] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

	double data_y[6] = {
	    1.0, 0.961538461538461, 0.862068965517241, 0.735294117647059, 0.609756097560976, 0.500000000000000};

	double test_x[50] = {0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24,
			     0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50,
			     0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76,
			     0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98};

	double test_y[50] = {
	    1.000000000000000, 0.997583282975581, 0.995079933416512, 0.992403318788142, 0.989466806555819,
	    0.986183764184894, 0.982467559140716, 0.978231558888635, 0.973389130893999, 0.967853642622158,
	    0.961538461538461, 0.954382579685350, 0.946427487413627, 0.937740299651188, 0.928388131325928,
	    0.918438097365742, 0.907957312698524, 0.897012892252170, 0.885671950954575, 0.874001603733634,
	    0.862068965517241, 0.849933363488199, 0.837622973848936, 0.825158185056786, 0.812559385569085,
	    0.799846963843167, 0.787041308336369, 0.774162807506023, 0.761231849809467, 0.748268823704033,
	    0.735294117647059, 0.722328486073082, 0.709394147325463, 0.696513685724764, 0.683709685591549,
	    0.671004731246381, 0.658421407009825, 0.645982297202442, 0.633709986144797, 0.621627058157454,
	    0.609756097560976, 0.598112015427308, 0.586679029833925, 0.575433685609685, 0.564352527583445,
	    0.553412100584061, 0.542588949440392, 0.531859618981294, 0.521200654035625, 0.510588599432241};

	double test_dy[50] = {
	    -0.120113913432180, -0.122279726798445, -0.128777166897241, -0.139606233728568, -0.154766927292426,
	    -0.174259247588814, -0.198083194617734, -0.226238768379184, -0.258725968873165, -0.295544796099676,
	    -0.336695250058719, -0.378333644186652, -0.416616291919835, -0.451543193258270, -0.483114348201955,
	    -0.511329756750890, -0.536189418905076, -0.557693334664512, -0.575841504029200, -0.590633926999137,
	    -0.602070603574326, -0.611319695518765, -0.619549364596455, -0.626759610807396, -0.632950434151589,
	    -0.638121834629033, -0.642273812239728, -0.645406366983674, -0.647519498860871, -0.648613207871319,
	    -0.648687494015019, -0.647687460711257, -0.645558211379322, -0.642299746019212, -0.637912064630930,
	    -0.632395167214473, -0.625749053769843, -0.617973724297039, -0.609069178796061, -0.599035417266910,
	    -0.587872439709585, -0.576731233416743, -0.566762785681043, -0.557967096502484, -0.550344165881066,
	    -0.543893993816790, -0.538616580309654, -0.534511925359660, -0.531580028966807, -0.529820891131095};

	double test_iy[50] = {
	    0.000000000000000, 0.019975905023535, 0.039902753768792, 0.059777947259733, 0.079597153869625,
	    0.099354309321042, 0.119041616685866, 0.138649546385285, 0.158166836189794, 0.177580491219196,
	    0.196875783942601, 0.216036382301310, 0.235045759060558, 0.253888601161251, 0.272550937842853,
	    0.291020140643388, 0.309284923399436, 0.327335342246135, 0.345162795617181, 0.362760024244829,
	    0.380121111159890, 0.397241442753010, 0.414117280448683, 0.430745332379281, 0.447122714446318,
	    0.463246950320456, 0.479115971441505, 0.494728117018421, 0.510082134029305, 0.525177177221407,
	    0.540012809111123, 0.554589001813881, 0.568906157172889, 0.582965126887879, 0.596767214344995,
	    0.610314174616794, 0.623608214462242, 0.636651992326715, 0.649448618342004, 0.662001654326309,
	    0.674315113784241, 0.686393423540581, 0.698241001711602, 0.709861835676399, 0.721259443710643,
	    0.732436874986582, 0.743396709573044, 0.754141058435429, 0.764671563435718, 0.774989397332469};

	s = test_interp(data_x, data_y, std::end(data_x) - std::begin(data_x), InterpolatorType::CUBIC_SPLINE_NATURAL,
			test_x, test_y, test_dy, test_iy, std::end(test_x) - std::begin(test_x));
	return s;
}

static int test_cspline3(void)
{
	/* This data has been chosen to be random (uneven spacing in x and y)
	and the exact cubic spine solution computed using Octave */

	int s;

	double data_x[7] = {-1.2139767065644265, -0.792590494453907, -0.250954683125019, 0.665867809951305,
			    0.735655088722706,   0.827622053027153,  1.426592227816582};

	double data_y[7] = {-0.00453877449035645, 0.49763182550668716, 0.17805472016334534, 0.40514493733644485,
			    -0.21595209836959839, 0.47405586764216423, 0.46561462432146072};

	double test_x[19] = {
	    -1.2139767065644265, -1.0735146358609200, -0.9330525651574135, -0.7925904944539071, -0.6120452240109444,
	    -0.4314999535679818, -0.2509546831250191, 0.0546528145670890,  0.3602603122591972,  0.6658678099513053,
	    0.6891302362084388,  0.7123926624655723,  0.7356550887227058,  0.7663107434908548,  0.7969663982590039,
	    0.8276220530271530,  1.0272787779569625,  1.2269355028867721,  1.4265922278165817,
	};

	double test_y[19] = {-0.00453877449035645, 0.25816917628390590,  0.44938881397673230, 0.49763182550668716,
			     0.31389980410075147,  0.09948951681196887,  0.17805472016334534, 1.27633142487980233,
			     2.04936553432792001,  0.40514493733644485,  0.13322324792901385, -0.09656315924697809,
			     -0.21595209836959839, -0.13551147728045118, 0.13466779030061801, 0.47405586764216423,
			     1.68064089899304370,  1.43594739539458649,  0.46561462432146072};

	double test_dy[19] = {
	    1.955137555965937,   1.700662049790549,  0.937235531264386,  -0.335141999612553, -1.401385073563169,
	    -0.674982149482761,  1.844066772628670,  4.202528085784793,  -0.284432022227558, -11.616813551408383,
	    -11.272731243226174, -7.994209291156876, -1.781247695200491, 6.373970868827501,  10.597456848997197,
	    10.889210245308570,  1.803124267866902,  -3.648527318598099, -5.465744514086432,
	};

	double test_iy[19] = {
	    0.000000000000000, 0.018231117234914, 0.069178822023139, 0.137781019634897, 0.213936442847744,
	    0.249280997744777, 0.267492946016120, 0.471372708120518, 1.014473660088477, 1.477731933018837,
	    1.483978291717981, 1.484256847945450, 1.480341742628893, 1.474315901028282, 1.473972210647307,
	    1.483279773396950, 1.728562698140330, 2.057796448999396, 2.253662886537457,
	};

	s = test_interp(data_x, data_y, std::end(data_x) - std::begin(data_x), InterpolatorType::CUBIC_SPLINE_NATURAL,
			test_x, test_y, test_dy, test_iy, std::end(test_x) - std::begin(test_x));
	return s;
}

// Following is data in triplets, x, rate, forward
static double testdata_monotone_convex[] = {
    0.000000000000, 0.025000000000, 0.025000000000, 0.010000000000, 0.025000500000, 0.025001500000, 0.020000000000,
    0.025002000000, 0.025006000000, 0.030000000000, 0.025004500000, 0.025013500000, 0.040000000000, 0.025008000000,
    0.025024000000, 0.050000000000, 0.025012500000, 0.025037500000, 0.060000000000, 0.025018000000, 0.025054000000,
    0.070000000000, 0.025024500000, 0.025073500000, 0.080000000000, 0.025032000000, 0.025096000000, 0.090000000000,
    0.025040500000, 0.025121500000, 0.100000000000, 0.025050000000, 0.025150000000, 0.110000000000, 0.025060500000,
    0.025181500000, 0.120000000000, 0.025072000000, 0.025216000000, 0.130000000000, 0.025084500000, 0.025253500000,
    0.140000000000, 0.025098000000, 0.025294000000, 0.150000000000, 0.025112500000, 0.025337500000, 0.160000000000,
    0.025128000000, 0.025384000000, 0.170000000000, 0.025144500000, 0.025433500000, 0.180000000000, 0.025162000000,
    0.025486000000, 0.190000000000, 0.025180500000, 0.025541500000, 0.200000000000, 0.025200000000, 0.025600000000,
    0.210000000000, 0.025220500000, 0.025661500000, 0.220000000000, 0.025242000000, 0.025726000000, 0.230000000000,
    0.025264500000, 0.025793500000, 0.240000000000, 0.025288000000, 0.025864000000, 0.250000000000, 0.025312500000,
    0.025937500000, 0.260000000000, 0.025338000000, 0.026014000000, 0.270000000000, 0.025364500000, 0.026093500000,
    0.280000000000, 0.025392000000, 0.026176000000, 0.290000000000, 0.025420500000, 0.026261500000, 0.300000000000,
    0.025450000000, 0.026350000000, 0.310000000000, 0.025480500000, 0.026441500000, 0.320000000000, 0.025512000000,
    0.026536000000, 0.330000000000, 0.025544500000, 0.026633500000, 0.340000000000, 0.025578000000, 0.026734000000,
    0.350000000000, 0.025612500000, 0.026837500000, 0.360000000000, 0.025648000000, 0.026944000000, 0.370000000000,
    0.025684500000, 0.027053500000, 0.380000000000, 0.025722000000, 0.027166000000, 0.390000000000, 0.025760500000,
    0.027281500000, 0.400000000000, 0.025800000000, 0.027400000000, 0.410000000000, 0.025840500000, 0.027521500000,
    0.420000000000, 0.025882000000, 0.027646000000, 0.430000000000, 0.025924500000, 0.027773500000, 0.440000000000,
    0.025968000000, 0.027904000000, 0.450000000000, 0.026012500000, 0.028037500000, 0.460000000000, 0.026058000000,
    0.028174000000, 0.470000000000, 0.026104500000, 0.028313500000, 0.480000000000, 0.026152000000, 0.028456000000,
    0.490000000000, 0.026200500000, 0.028601500000, 0.500000000000, 0.026250000000, 0.028750000000, 0.510000000000,
    0.026300500000, 0.028901500000, 0.520000000000, 0.026352000000, 0.029056000000, 0.530000000000, 0.026404500000,
    0.029213500000, 0.540000000000, 0.026458000000, 0.029374000000, 0.550000000000, 0.026512500000, 0.029537500000,
    0.560000000000, 0.026568000000, 0.029704000000, 0.570000000000, 0.026624500000, 0.029873500000, 0.580000000000,
    0.026682000000, 0.030046000000, 0.590000000000, 0.026740500000, 0.030221500000, 0.600000000000, 0.026800000000,
    0.030400000000, 0.610000000000, 0.026860500000, 0.030581500000, 0.620000000000, 0.026922000000, 0.030766000000,
    0.630000000000, 0.026984500000, 0.030953500000, 0.640000000000, 0.027048000000, 0.031144000000, 0.650000000000,
    0.027112500000, 0.031337500000, 0.660000000000, 0.027178000000, 0.031534000000, 0.670000000000, 0.027244500000,
    0.031733500000, 0.680000000000, 0.027312000000, 0.031936000000, 0.690000000000, 0.027380500000, 0.032141500000,
    0.700000000000, 0.027450000000, 0.032350000000, 0.710000000000, 0.027520500000, 0.032561500000, 0.720000000000,
    0.027592000000, 0.032776000000, 0.730000000000, 0.027664500000, 0.032993500000, 0.740000000000, 0.027738000000,
    0.033214000000, 0.750000000000, 0.027812500000, 0.033437500000, 0.760000000000, 0.027888000000, 0.033664000000,
    0.770000000000, 0.027964500000, 0.033893500000, 0.780000000000, 0.028042000000, 0.034126000000, 0.790000000000,
    0.028120500000, 0.034361500000, 0.800000000000, 0.028200000000, 0.034600000000, 0.810000000000, 0.028280500000,
    0.034841500000, 0.820000000000, 0.028362000000, 0.035086000000, 0.830000000000, 0.028444500000, 0.035333500000,
    0.840000000000, 0.028528000000, 0.035584000000, 0.850000000000, 0.028612500000, 0.035837500000, 0.860000000000,
    0.028698000000, 0.036094000000, 0.870000000000, 0.028784500000, 0.036353500000, 0.880000000000, 0.028872000000,
    0.036616000000, 0.890000000000, 0.028960500000, 0.036881500000, 0.900000000000, 0.029050000000, 0.037150000000,
    0.910000000000, 0.029140500000, 0.037421500000, 0.920000000000, 0.029232000000, 0.037696000000, 0.930000000000,
    0.029324500000, 0.037973500000, 0.940000000000, 0.029418000000, 0.038254000000, 0.950000000000, 0.029512500000,
    0.038537500000, 0.960000000000, 0.029608000000, 0.038824000000, 0.970000000000, 0.029704500000, 0.039113500000,
    0.980000000000, 0.029802000000, 0.039406000000, 0.990000000000, 0.029900500000, 0.039701500000, 1.000000000000,
    0.030000000000, 0.040000000000, 1.010000000000, 0.030100441089, 0.040288650000, 1.020000000000, 0.030201729412,
    0.040574600000, 1.030000000000, 0.030303814078, 0.040857850000, 1.040000000000, 0.030406646154, 0.041138400000,
    1.050000000000, 0.030510178571, 0.041416250000, 1.060000000000, 0.030614366038, 0.041691400000, 1.070000000000,
    0.030719164953, 0.041963850000, 1.080000000000, 0.030824533333, 0.042233600000, 1.090000000000, 0.030930430734,
    0.042500650000, 1.100000000000, 0.031036818182, 0.042765000000, 1.110000000000, 0.031143658108, 0.043026650000,
    1.120000000000, 0.031250914286, 0.043285600000, 1.130000000000, 0.031358551770, 0.043541850000, 1.140000000000,
    0.031466536842, 0.043795400000, 1.150000000000, 0.031574836957, 0.044046250000, 1.160000000000, 0.031683420690,
    0.044294400000, 1.170000000000, 0.031792257692, 0.044539850000, 1.180000000000, 0.031901318644, 0.044782600000,
    1.190000000000, 0.032010575210, 0.045022650000, 1.200000000000, 0.032120000000, 0.045260000000, 1.210000000000,
    0.032229566529, 0.045494650000, 1.220000000000, 0.032339249180, 0.045726600000, 1.230000000000, 0.032449023171,
    0.045955850000, 1.240000000000, 0.032558864516, 0.046182400000, 1.250000000000, 0.032668750000, 0.046406250000,
    1.260000000000, 0.032778657143, 0.046627400000, 1.270000000000, 0.032888564173, 0.046845850000, 1.280000000000,
    0.032998450000, 0.047061600000, 1.290000000000, 0.033108294186, 0.047274650000, 1.300000000000, 0.033218076923,
    0.047485000000, 1.310000000000, 0.033327779008, 0.047692650000, 1.320000000000, 0.033437381818, 0.047897600000,
    1.330000000000, 0.033546867293, 0.048099850000, 1.340000000000, 0.033656217910, 0.048299400000, 1.350000000000,
    0.033765416667, 0.048496250000, 1.360000000000, 0.033874447059, 0.048690400000, 1.370000000000, 0.033983293066,
    0.048881850000, 1.380000000000, 0.034091939130, 0.049070600000, 1.390000000000, 0.034200370144, 0.049256650000,
    1.400000000000, 0.034308571429, 0.049440000000, 1.410000000000, 0.034416528723, 0.049620650000, 1.420000000000,
    0.034524228169, 0.049798600000, 1.430000000000, 0.034631656294, 0.049973850000, 1.440000000000, 0.034738800000,
    0.050146400000, 1.450000000000, 0.034845646552, 0.050316250000, 1.460000000000, 0.034952183562, 0.050483400000,
    1.470000000000, 0.035058398980, 0.050647850000, 1.480000000000, 0.035164281081, 0.050809600000, 1.490000000000,
    0.035269818456, 0.050968650000, 1.500000000000, 0.035375000000, 0.051125000000, 1.510000000000, 0.035479814901,
    0.051278650000, 1.520000000000, 0.035584252632, 0.051429600000, 1.530000000000, 0.035688302941, 0.051577850000,
    1.540000000000, 0.035791955844, 0.051723400000, 1.550000000000, 0.035895201613, 0.051866250000, 1.560000000000,
    0.035998030769, 0.052006400000, 1.570000000000, 0.036100434076, 0.052143850000, 1.580000000000, 0.036202402532,
    0.052278600000, 1.590000000000, 0.036303927358, 0.052410650000, 1.600000000000, 0.036405000000, 0.052540000000,
    1.610000000000, 0.036505612112, 0.052666650000, 1.620000000000, 0.036605755556, 0.052790600000, 1.630000000000,
    0.036705422393, 0.052911850000, 1.640000000000, 0.036804604878, 0.053030400000, 1.650000000000, 0.036903295455,
    0.053146250000, 1.660000000000, 0.037001486747, 0.053259400000, 1.670000000000, 0.037099171557, 0.053369850000,
    1.680000000000, 0.037196342857, 0.053477600000, 1.690000000000, 0.037292993787, 0.053582650000, 1.700000000000,
    0.037389117647, 0.053685000000, 1.710000000000, 0.037484707895, 0.053784650000, 1.720000000000, 0.037579758140,
    0.053881600000, 1.730000000000, 0.037674262139, 0.053975850000, 1.740000000000, 0.037768213793, 0.054067400000,
    1.750000000000, 0.037861607143, 0.054156250000, 1.760000000000, 0.037954436364, 0.054242400000, 1.770000000000,
    0.038046695763, 0.054325850000, 1.780000000000, 0.038138379775, 0.054406600000, 1.790000000000, 0.038229482961,
    0.054484650000, 1.800000000000, 0.038320000000, 0.054560000000, 1.810000000000, 0.038409925691, 0.054632650000,
    1.820000000000, 0.038499254945, 0.054702600000, 1.830000000000, 0.038587982787, 0.054769850000, 1.840000000000,
    0.038676104348, 0.054834400000, 1.850000000000, 0.038763614865, 0.054896250000, 1.860000000000, 0.038850509677,
    0.054955400000, 1.870000000000, 0.038936784225, 0.055011850000, 1.880000000000, 0.039022434043, 0.055065600000,
    1.890000000000, 0.039107454762, 0.055116650000, 1.900000000000, 0.039191842105, 0.055165000000, 1.910000000000,
    0.039275591885, 0.055210650000, 1.920000000000, 0.039358700000, 0.055253600000, 1.930000000000, 0.039441162435,
    0.055293850000, 1.940000000000, 0.039522975258, 0.055331400000, 1.950000000000, 0.039604134615, 0.055366250000,
    1.960000000000, 0.039684636735, 0.055398400000, 1.970000000000, 0.039764477919, 0.055427850000, 1.980000000000,
    0.039843654545, 0.055454600000, 1.990000000000, 0.039922163065, 0.055478650000, 2.000000000000, 0.040000000000,
    0.055500000000, 2.010000000000, 0.040077114428, 0.055500000000, 2.020000000000, 0.040153465347, 0.055500000000,
    2.030000000000, 0.040229064039, 0.055500000000, 2.040000000000, 0.040303921569, 0.055500000000, 2.050000000000,
    0.040378048780, 0.055500000000, 2.060000000000, 0.040451456311, 0.055500000000, 2.070000000000, 0.040524154589,
    0.055500000000, 2.080000000000, 0.040596153846, 0.055500000000, 2.090000000000, 0.040667464115, 0.055500000000,
    2.100000000000, 0.040738095238, 0.055500000000, 2.110000000000, 0.040808056872, 0.055500000000, 2.120000000000,
    0.040877358491, 0.055500000000, 2.130000000000, 0.040946009390, 0.055500000000, 2.140000000000, 0.041014018692,
    0.055500000000, 2.150000000000, 0.041081395349, 0.055500000000, 2.160000000000, 0.041148148148, 0.055500000000,
    2.170000000000, 0.041214285714, 0.055500000000, 2.180000000000, 0.041279816514, 0.055500000000, 2.190000000000,
    0.041344748858, 0.055500000000, 2.200000000000, 0.041409090909, 0.055500000000, 2.210000000000, 0.041472850679,
    0.055500000000, 2.220000000000, 0.041536036036, 0.055500000000, 2.230000000000, 0.041598654709, 0.055500000000,
    2.240000000000, 0.041660714286, 0.055500000000, 2.250000000000, 0.041722222222, 0.055500000000, 2.260000000000,
    0.041783185841, 0.055500000000, 2.270000000000, 0.041843612335, 0.055500000000, 2.280000000000, 0.041903508772,
    0.055500000000, 2.290000000000, 0.041962882096, 0.055500000000, 2.300000000000, 0.042021739130, 0.055500000000,
    2.310000000000, 0.042080086580, 0.055500000000, 2.320000000000, 0.042137931034, 0.055500000000, 2.330000000000,
    0.042195279293, 0.055500650184, 2.340000000000, 0.042252155556, 0.055509800000, 2.350000000000, 0.042308609686,
    0.055529753214, 2.360000000000, 0.042364692829, 0.055560509826, 2.370000000000, 0.042420455262, 0.055602069835,
    2.380000000000, 0.042475946422, 0.055654433242, 2.390000000000, 0.042531214915, 0.055717600046, 2.400000000000,
    0.042586308540, 0.055791570248, 2.410000000000, 0.042641274300, 0.055876343848, 2.420000000000, 0.042696158423,
    0.055971920845, 2.430000000000, 0.042751006375, 0.056078301240, 2.440000000000, 0.042805862877, 0.056195485032,
    2.450000000000, 0.042860771920, 0.056323472222, 2.460000000000, 0.042915776779, 0.056462262810, 2.470000000000,
    0.042970920029, 0.056611856795, 2.480000000000, 0.043026243558, 0.056772254178, 2.490000000000, 0.043081788581,
    0.056943454959, 2.500000000000, 0.043137595654, 0.057125459137, 2.510000000000, 0.043193704685, 0.057318266713,
    2.520000000000, 0.043250154952, 0.057521877686, 2.530000000000, 0.043306985109, 0.057736292057, 2.540000000000,
    0.043364233201, 0.057961509826, 2.550000000000, 0.043421936680, 0.058197530992, 2.560000000000, 0.043480132407,
    0.058444355556, 2.570000000000, 0.043538856675, 0.058701983517, 2.580000000000, 0.043598145211, 0.058970414876,
    2.590000000000, 0.043658033191, 0.059249649633, 2.600000000000, 0.043718555249, 0.059539687787, 2.610000000000,
    0.043779745491, 0.059840529339, 2.620000000000, 0.043841637499, 0.060152174288, 2.630000000000, 0.043904264345,
    0.060474622635, 2.640000000000, 0.043967658603, 0.060807874380, 2.650000000000, 0.044031852350, 0.061151929522,
    2.660000000000, 0.044096877185, 0.061506788062, 2.670000000000, 0.044162764232, 0.061872450000, 2.680000000000,
    0.044229544151, 0.062248915335, 2.690000000000, 0.044297247144, 0.062636184068, 2.700000000000, 0.044365902969,
    0.063034256198, 2.710000000000, 0.044435540942, 0.063443131726, 2.720000000000, 0.044506189949, 0.063862810652,
    2.730000000000, 0.044577878453, 0.064293292975, 2.740000000000, 0.044650634501, 0.064734578696, 2.750000000000,
    0.044724485732, 0.065186667815, 2.760000000000, 0.044799459384, 0.065649560331, 2.770000000000, 0.044875582304,
    0.066123256244, 2.780000000000, 0.044952880949, 0.066607755556, 2.790000000000, 0.045031381399, 0.067103058264,
    2.800000000000, 0.045111109362, 0.067609164371, 2.810000000000, 0.045192090179, 0.068126073875, 2.820000000000,
    0.045274348831, 0.068653786777, 2.830000000000, 0.045357909947, 0.069192303076, 2.840000000000, 0.045442797808,
    0.069741622773, 2.850000000000, 0.045529036356, 0.070301745868, 2.860000000000, 0.045616649197, 0.070872672360,
    2.870000000000, 0.045705659608, 0.071454402250, 2.880000000000, 0.045796090542, 0.072046935537, 2.890000000000,
    0.045887964635, 0.072650272222, 2.900000000000, 0.045981304212, 0.073264412305, 2.910000000000, 0.046076131290,
    0.073889355785, 2.920000000000, 0.046172467583, 0.074525102663, 2.930000000000, 0.046270334511, 0.075171652938,
    2.940000000000, 0.046369753202, 0.075829006612, 2.950000000000, 0.046470744496, 0.076497163682, 2.960000000000,
    0.046573328952, 0.077176124151, 2.970000000000, 0.046677526854, 0.077865888017, 2.980000000000, 0.046783358211,
    0.078566455280, 2.990000000000, 0.046890842767, 0.079277825941, 3.000000000000, 0.047000000000, 0.080000000000,
    3.010000000000, 0.047111499698, 0.081119098751, 3.020000000000, 0.047225929666, 0.082215882183, 3.030000000000,
    0.047343187244, 0.083290350296, 3.040000000000, 0.047463171122, 0.084342503090, 3.050000000000, 0.047585781320,
    0.085372340565, 3.060000000000, 0.047710919163, 0.086379862722, 3.070000000000, 0.047838487263, 0.087365069559,
    3.080000000000, 0.047968389495, 0.088327961078, 3.090000000000, 0.048100530982, 0.089268537278, 3.100000000000,
    0.048234818067, 0.090186798159, 3.110000000000, 0.048371158300, 0.091082743721, 3.120000000000, 0.048509460416,
    0.091956373964, 3.130000000000, 0.048649634315, 0.092807688889, 3.140000000000, 0.048791591048, 0.093636688494,
    3.150000000000, 0.048935242791, 0.094443372781, 3.160000000000, 0.049080502836, 0.095227741749, 3.170000000000,
    0.049227285566, 0.095989795398, 3.180000000000, 0.049375506442, 0.096729533728, 3.190000000000, 0.049525081985,
    0.097446956739, 3.200000000000, 0.049675929761, 0.098142064431, 3.210000000000, 0.049827968361, 0.098814856805,
    3.220000000000, 0.049981117387, 0.099465333859, 3.230000000000, 0.050135297440, 0.100093495595, 3.240000000000,
    0.050290430097, 0.100699342012, 3.250000000000, 0.050446437904, 0.101282873110, 3.260000000000, 0.050603244354,
    0.101844088889, 3.270000000000, 0.050760773878, 0.102382989349, 3.280000000000, 0.050918951828, 0.102899574490,
    3.290000000000, 0.051077704464, 0.103393844313, 3.300000000000, 0.051236958938, 0.103865798817, 3.310000000000,
    0.051396643286, 0.104315438001, 3.320000000000, 0.051556686407, 0.104742761867, 3.330000000000, 0.051717018057,
    0.105147770414, 3.340000000000, 0.051877568831, 0.105530463642, 3.350000000000, 0.052038270155, 0.105890841552,
    3.360000000000, 0.052199054269, 0.106228904142, 3.370000000000, 0.052359854219, 0.106544651414, 3.380000000000,
    0.052520603842, 0.106838083366, 3.390000000000, 0.052681237758, 0.107109200000, 3.400000000000, 0.052841691354,
    0.107358001315, 3.410000000000, 0.053001900775, 0.107584487311, 3.420000000000, 0.053161802914, 0.107788657988,
    3.430000000000, 0.053321335399, 0.107970513346, 3.440000000000, 0.053480436583, 0.108130053386, 3.450000000000,
    0.053639045536, 0.108267278107, 3.460000000000, 0.053797102030, 0.108382187508, 3.470000000000, 0.053954546532,
    0.108474781591, 3.480000000000, 0.054111320193, 0.108545060355, 3.490000000000, 0.054267364839, 0.108593023800,
    3.500000000000, 0.054422622961, 0.108618671926, 3.510000000000, 0.054577037672, 0.108621906337, 3.520000000000,
    0.054730550995, 0.108601562327, 3.530000000000, 0.054883103047, 0.108557302458, 3.540000000000, 0.055034634416,
    0.108489126731, 3.550000000000, 0.055185086358, 0.108397035145, 3.560000000000, 0.055334400791, 0.108281027701,
    3.570000000000, 0.055482520281, 0.108141104398, 3.580000000000, 0.055629388040, 0.107977265235, 3.590000000000,
    0.055774947908, 0.107789510215, 3.600000000000, 0.055919144352, 0.107577839335, 3.610000000000, 0.056061922453,
    0.107342252597, 3.620000000000, 0.056203227901, 0.107082750000, 3.630000000000, 0.056343006980, 0.106799331544,
    3.640000000000, 0.056481206569, 0.106491997230, 3.650000000000, 0.056617774127, 0.106160747057, 3.660000000000,
    0.056752657687, 0.105805581025, 3.670000000000, 0.056885805848, 0.105426499134, 3.680000000000, 0.057017167771,
    0.105023501385, 3.690000000000, 0.057146693163, 0.104596587777, 3.700000000000, 0.057274332279, 0.104145758310,
    3.710000000000, 0.057400035908, 0.103671012985, 3.720000000000, 0.057523755369, 0.103172351801, 3.730000000000,
    0.057645442503, 0.102649774758, 3.740000000000, 0.057765049665, 0.102103281856, 3.750000000000, 0.057882529721,
    0.101532873096, 3.760000000000, 0.057997836035, 0.100938548476, 3.770000000000, 0.058110922468, 0.100320307999,
    3.780000000000, 0.058221743368, 0.099678151662, 3.790000000000, 0.058330253566, 0.099012079467, 3.800000000000,
    0.058436408369, 0.098322091413, 3.810000000000, 0.058540163550, 0.097608187500, 3.820000000000, 0.058641475348,
    0.096870367729, 3.830000000000, 0.058740300460, 0.096108632098, 3.840000000000, 0.058836596030, 0.095322980609,
    3.850000000000, 0.058930319650, 0.094513413262, 3.860000000000, 0.059021429352, 0.093679930055, 3.870000000000,
    0.059109883601, 0.092822530990, 3.880000000000, 0.059195641289, 0.091941216066, 3.890000000000, 0.059278661732,
    0.091035985284, 3.900000000000, 0.059358904663, 0.090106838643, 3.910000000000, 0.059436330227, 0.089153776143,
    3.920000000000, 0.059510898977, 0.088176797784, 3.930000000000, 0.059582571864, 0.087175903566, 3.940000000000,
    0.059651310239, 0.086151093490, 3.950000000000, 0.059717075842, 0.085102367555, 3.960000000000, 0.059779830801,
    0.084029725762, 3.970000000000, 0.059839537625, 0.082933168109, 3.980000000000, 0.059896159200, 0.081812694598,
    3.990000000000, 0.059949658783, 0.080668305229, 4.000000000000, 0.060000000000, 0.079500000000, 4.010000000000,
    0.060047901434, 0.078917925000, 4.020000000000, 0.060094123881, 0.078341700000, 4.030000000000, 0.060138694355,
    0.077771325000, 4.040000000000, 0.060181639604, 0.077206800000, 4.050000000000, 0.060222986111, 0.076648125000,
    4.060000000000, 0.060262760099, 0.076095300000, 4.070000000000, 0.060300987531, 0.075548325000, 4.080000000000,
    0.060337694118, 0.075007200000, 4.090000000000, 0.060372905318, 0.074471925000, 4.100000000000, 0.060406646341,
    0.073942500000, 4.110000000000, 0.060438942153, 0.073418925000, 4.120000000000, 0.060469817476, 0.072901200000,
    4.130000000000, 0.060499296792, 0.072389325000, 4.140000000000, 0.060527404348, 0.071883300000, 4.150000000000,
    0.060554164157, 0.071383125000, 4.160000000000, 0.060579600000, 0.070888800000, 4.170000000000, 0.060603735432,
    0.070400325000, 4.180000000000, 0.060626593780, 0.069917700000, 4.190000000000, 0.060648198150, 0.069440925000,
    4.200000000000, 0.060668571429, 0.068970000000, 4.210000000000, 0.060687736283, 0.068504925000, 4.220000000000,
    0.060705715166, 0.068045700000, 4.230000000000, 0.060722530319, 0.067592325000, 4.240000000000, 0.060738203774,
    0.067144800000, 4.250000000000, 0.060752757353, 0.066703125000, 4.260000000000, 0.060766212676, 0.066267300000,
    4.270000000000, 0.060778591159, 0.065837325000, 4.280000000000, 0.060789914019, 0.065413200000, 4.290000000000,
    0.060800202273, 0.064994925000, 4.300000000000, 0.060809476744, 0.064582500000, 4.310000000000, 0.060817758063,
    0.064175925000, 4.320000000000, 0.060825066667, 0.063775200000, 4.330000000000, 0.060831422806, 0.063380325000,
    4.340000000000, 0.060836846544, 0.062991300000, 4.350000000000, 0.060841357759, 0.062608125000, 4.360000000000,
    0.060844976147, 0.062230800000, 4.370000000000, 0.060847721224, 0.061859325000, 4.380000000000, 0.060849612329,
    0.061493700000, 4.390000000000, 0.060850668622, 0.061133925000, 4.400000000000, 0.060850909091, 0.060780000000,
    4.410000000000, 0.060850352551, 0.060431925000, 4.420000000000, 0.060849017647, 0.060089700000, 4.430000000000,
    0.060846922856, 0.059753325000, 4.440000000000, 0.060844086486, 0.059422800000, 4.450000000000, 0.060840526685,
    0.059098125000, 4.460000000000, 0.060836261435, 0.058779300000, 4.470000000000, 0.060831308557, 0.058466325000,
    4.480000000000, 0.060825685714, 0.058159200000, 4.490000000000, 0.060819410412, 0.057857925000, 4.500000000000,
    0.060812500000, 0.057562500000, 4.510000000000, 0.060804971674, 0.057272925000, 4.520000000000, 0.060796842478,
    0.056989200000, 4.530000000000, 0.060788129305, 0.056711325000, 4.540000000000, 0.060778848899, 0.056439300000,
    4.550000000000, 0.060769017857, 0.056173125000, 4.560000000000, 0.060758652632, 0.055912800000, 4.570000000000,
    0.060747769530, 0.055658325000, 4.580000000000, 0.060736384716, 0.055409700000, 4.590000000000, 0.060724514216,
    0.055166925000, 4.600000000000, 0.060712173913, 0.054930000000, 4.610000000000, 0.060699379555, 0.054698925000,
    4.620000000000, 0.060686146753, 0.054473700000, 4.630000000000, 0.060672490983, 0.054254325000, 4.640000000000,
    0.060658427586, 0.054040800000, 4.650000000000, 0.060643971774, 0.053833125000, 4.660000000000, 0.060629138627,
    0.053631300000, 4.670000000000, 0.060613943094, 0.053435325000, 4.680000000000, 0.060598400000, 0.053245200000,
    4.690000000000, 0.060582524041, 0.053060925000, 4.700000000000, 0.060566329787, 0.052882500000, 4.710000000000,
    0.060549831688, 0.052709925000, 4.720000000000, 0.060533044068, 0.052543200000, 4.730000000000, 0.060515981131,
    0.052382325000, 4.740000000000, 0.060498656962, 0.052227300000, 4.750000000000, 0.060481085526, 0.052078125000,
    4.760000000000, 0.060463280672, 0.051934800000, 4.770000000000, 0.060445256132, 0.051797325000, 4.780000000000,
    0.060427025523, 0.051665700000, 4.790000000000, 0.060408602349, 0.051539925000, 4.800000000000, 0.060390000000,
    0.051420000000, 4.810000000000, 0.060371231757, 0.051305925000, 4.820000000000, 0.060352310788, 0.051197700000,
    4.830000000000, 0.060333250155, 0.051095325000, 4.840000000000, 0.060314062810, 0.050998800000, 4.850000000000,
    0.060294761598, 0.050908125000, 4.860000000000, 0.060275359259, 0.050823300000, 4.870000000000, 0.060255868429,
    0.050744325000, 4.880000000000, 0.060236301639, 0.050671200000, 4.890000000000, 0.060216671319, 0.050603925000,
    4.900000000000, 0.060196989796, 0.050542500000, 4.910000000000, 0.060177269297, 0.050486925000, 4.920000000000,
    0.060157521951, 0.050437200000, 4.930000000000, 0.060137759787, 0.050393325000, 4.940000000000, 0.060117994737,
    0.050355300000, 4.950000000000, 0.060098238636, 0.050323125000, 4.960000000000, 0.060078503226, 0.050296800000,
    4.970000000000, 0.060058800151, 0.050276325000, 4.980000000000, 0.060039140964, 0.050261700000, 4.990000000000,
    0.060019537124, 0.050252925000, 5.000000000000, 0.060000000000, 0.050250000000, 5.010000000000, 0.059980538922,
    0.050250000000, 5.020000000000, 0.059961155378, 0.050250000000, 5.030000000000, 0.059941848907, 0.050250000000,
    5.040000000000, 0.059922619048, 0.050250000000, 5.050000000000, 0.059903465347, 0.050250000000, 5.060000000000,
    0.059884387352, 0.050250000000, 5.070000000000, 0.059865384615, 0.050250000000, 5.080000000000, 0.059846456693,
    0.050250000000, 5.090000000000, 0.059827603143, 0.050250000000, 5.100000000000, 0.059808823529, 0.050250000000,
    5.110000000000, 0.059790117417, 0.050250000000, 5.120000000000, 0.059771484375, 0.050250000000, 5.130000000000,
    0.059752923977, 0.050250000000, 5.140000000000, 0.059734435798, 0.050250000000, 5.150000000000, 0.059716019417,
    0.050250000000, 5.160000000000, 0.059697674419, 0.050250000000, 5.170000000000, 0.059679400387, 0.050250000000,
    5.180000000000, 0.059661196911, 0.050250000000, 5.190000000000, 0.059643063584, 0.050250000000, 5.200000000000,
    0.059625000000, 0.050250000000};

int test_monotone_convex_iterpolator()
{
	double data_x[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
	double data_y[6] = {0.0, 0.03, 0.04, 0.047, 0.06, 0.06};
	FixedRegionAllocator *alloc = get_threadspecific_allocators()->tempspace_allocator;
	int failure_count = 0;
	const int size = 6;
	MonotoneConvexInterpolator<redukti_adouble_t> interp(data_x, data_y, size, false, false,
							     InterpolatorType::MONOTONE_CONVEX, 1, true);
	auto p = std::begin(testdata_monotone_convex);
	auto endp = std::end(testdata_monotone_convex);
	for (double *ptr = p; ptr < endp; ptr += 3) {
		FixedRegionAllocatorGuard g(alloc);
		auto rate = interp.interpolate_with_sensitivities(ptr[0], alloc);
		if (unequal(redukti_adouble_get_value(rate.get()), ptr[1])) {
			failure_count++;
		}
		auto forward = interp.forward_with_sensitivities(ptr[0], alloc);
		if (unequal(redukti_adouble_get_value(forward.get()), ptr[2])) {
			failure_count++;
		}
	}
	{
		FixedRegionAllocatorGuard g(alloc);
		auto rptr = interp.interpolate_with_sensitivities(2.21456, alloc);
		if (unequal(redukti_adouble_get_value(rptr.get()), 0.041501733978)) {
			failure_count++;
		}
	}
	// fprintf(stderr, "Status\tTerm\tPillar\tNumeric\tAutodiff\n");
	int d_failure_count = 0;
#if 0
	int d_success_count = 0;
	for (size_t i = 0; i < size; i++) {
		for (int t = 10; t > 0; t--) {
			for (size_t pillar = 1; pillar < size; pillar++) {
				double term;
				if (i == 0) {
					term = data_x[1] / (double(t) / 10.0);
				} else {
					term = data_x[i - 1] +
					       (data_x[i] - data_x[i - 1]) *
						   (double(t) / 10.0);
				}
				if (term < 1.0 || term > 5.0)
					continue;
				FixedRegionAllocatorGuard g(alloc);
				double expected =
				    interp.sensitivity_to(term, pillar);
				auto autodiff =
				    interp.interpolate_with_sensitivities(
					term, alloc);
				double result = redukti_adouble_get_derivative1(
				    autodiff.get(), pillar);
				double err = ::fabs(expected - result);
				if (err > 1e-5) {
					fprintf(
					    stderr,
					    "Break\t%.12f\t%d\t%.12f\t%.12f\n",
					    term, (int)pillar, expected,
					    result);
					d_failure_count++;
				} else {
					d_success_count++;
				}
			}
		}
	}
#endif
	// fprintf(stderr,
	//        "Monotone Convex: comparison of autodiff versus numeric %d
	//        matches "
	//        "%d breaks in total\n",
	//        d_success_count, d_failure_count);

	failure_count += d_failure_count;

	return failure_count;
}

int test_interpolators()
{
	int rc = 0;
	rc += test_simple_cubicspline();
	rc += testSplineErrorOnGaussianValues();
	rc += testSplineOnRPN15AValues();
	rc += testSplineOnGenericValues();
	rc += testSimmetricEndConditions();
	rc += testDerivativeEndConditions();
	rc += test_linear();
	rc += test_linear_sensitivities();
	rc += test_cspline();
	rc += test_cspline2();
	rc += test_cspline3();
	rc += test_monotone_convex_iterpolator();
	if (rc == 0)
		printf("Interpolator Tests OK\n");
	else
		printf("Interpolator Tests FAILED\n");
	return rc;
}

} // namespace redukti
