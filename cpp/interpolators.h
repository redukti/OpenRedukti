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
/**
 * Portions derived from Quantlib.
 * License: http://quantlib.org/license.shtml
 */

#ifndef _REDUKTI_INTERPOLATORS_H
#define _REDUKTI_INTERPOLATORS_H

#include <allocators.h>
#include <autodiff.h>
#include <enums.pb.h>

#include <memory>
#include <vector>

namespace redukti
{

struct InterpolationOptions;
typedef std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>> SensitivitiesPointerType; // To help Cython

class Interpolator
{
	public:
	virtual ~Interpolator() {}

	// Interpolate at x
	virtual double interpolate(double x) = 0;
	double operator()(double x) { return interpolate(x); }

	// Interpolate at x
	// And also compute sensitivities of value at x
	// to the various terms in the data set.
	// Both first order and second order sensitivities
	// can be computed depending upon how the
	// the interpolator was created.
	// Uses automatic differentiation
	virtual SensitivitiesPointerType interpolate_with_sensitivities(double x, FixedRegionAllocator *A) = 0;

	// Interpolate at x
	// And also compute sensitivities of value at x
	// to the various terms in the data set.
	// Both first order and second order sensitivities
	// can be computed depending upon how the
	// the interpolator was created.
	// Uses numeric differentiation
	virtual SensitivitiesPointerType interpolate_with_numeric_sensitivities(double x, FixedRegionAllocator *A) = 0;

	// If underlying values have changed, this
	// method can be called to reinitialise the
	// interpolator.
	virtual void update() = 0;

	// Only available on Monotone Convex interpolator as it is an
	// interest rate aware interpolator - for everything else
	// an exception will be thrown.
	virtual double forward(double x)
	{
		assert(false);
		return 0.0;
	}
	virtual bool coefficients(std::vector<double> &cvec) { return false; }

	// Following are for testing only - not guaranteed
	// to give meaningful results for all interpolators
	virtual double primitive(double) { return 0.0; }
	virtual double derivative(double) { return 0.0; }
	virtual double derivative2(double) { return 0.0; }

	// Return the interpolator type
	virtual InterpolatorType type() const = 0;

	// Returns 0 if derivatives are not enabled
	// Returns 1 if first order derivatives are enabled
	// Returns 2 if both first and second order derivatives are enabled
	virtual int order() const = 0;

	// Returns the options that are enabled

	virtual void get_options(InterpolationOptions &optons) const = 0;
};

struct InterpolationOptions {
	bool monotoneconvex_inputs_are_forwards;
	double cubic_left_condition_value;
	double cubic_right_condition_value;
	bool extrapolate;
	int differentiation_order;

	InterpolationOptions()
	    : monotoneconvex_inputs_are_forwards(false), cubic_left_condition_value(0.0),
	      cubic_right_condition_value(0.0), extrapolate(false), differentiation_order(0)
	{
	}
};

typedef std::unique_ptr<Interpolator, Deleter<Interpolator>> InterpolatorPointerType; // To help Cython
// Return an Interpolator of the desired type.
// The x and y arrays will be referenced by the Interpolator,
// and therefore the caller must carefully manage
// changes.
extern InterpolatorPointerType make_interpolator(InterpolatorType type, double *x, double *y, unsigned int size,
						 Allocator *A = get_default_allocator(),
						 const InterpolationOptions &options = InterpolationOptions());

extern int test_interpolators();

} // namespace redukti

#endif
