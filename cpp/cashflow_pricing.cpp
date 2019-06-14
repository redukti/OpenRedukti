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

#include <cashflow_pricing.h>

#include <autodiff.h>
#include <cashflow.h>
#include <index.h>
#include <schedule.h>

#include <cashflow.pb.h>
#include <internal/autodiff_helpers.h>
#include <internal/cashflow_internal.h>
#include <internal/cashflow_pricing_internal.h>

#include <map>
#include <sstream>
#include <vector>

#include <assert.h>
#include <cctype>
#include <cmath>
#include <cstdio>

namespace redukti
{

// Helper to map from variable type to
// aggregation function
template <typename V> struct ValuationContextHelper;

// Specialization for double that is a no-op.
template <> struct ValuationContextHelper<double> {
	// No need to aggregate as sensitivities are not calculated
	void aggregateSensitivity(FixedRegionAllocator *A, const VariableResolver *discount_resolver,
				  const CurveReference *discount_curve, const VariableResolver *forward1_resolver,
				  const CurveReference *forward_curve1, const VariableResolver *forward2_resolver,
				  const CurveReference *forward_curve2, Sensitivities *container, double data,
				  StatusCode &status)
	{
	}
};

// Specialization for autodiff
template <> struct ValuationContextHelper<redukti_adouble_t> {
	void aggregateSensitivity(FixedRegionAllocator *A, const VariableResolver *discount_resolver,
				  const CurveReference *discount_curve, const VariableResolver *forward1_resolver,
				  const CurveReference *forward_curve1, const VariableResolver *forward2_resolver,
				  const CurveReference *forward_curve2, Sensitivities *container,
				  redukti_adouble_t &data, StatusCode &status)
	{
		aggregate_sensitivity(A, discount_resolver, discount_curve, forward1_resolver, forward_curve1,
				      forward2_resolver, forward_curve2, container, data, status);
	}
};

// Utility for setting manually computed derivatives for the sake
// of efficiency
static void adouble_exp_init(redukti_adouble_t *v, size_t num_variables, int order, int index, double v1,
			     double c1) noexcept
{
	assert(index >= 0);
	assert(index < (int)num_variables);
	assert(order >= 0);
	assert(order <= 2);
	double scalar = std::exp(v1 * c1);
	redukti_adouble_init(v, num_variables, order, index, scalar);
	if (order && index >= 0 && num_variables > 0) {
		double d1 = c1 * scalar;
		redukti_adouble_set_derivative1(v, index, d1);
		if (order == 2)
			redukti_adouble_set_derivative2(v, index, index, c1 * d1);
	}
}

// Utility for setting manually computed derivatives for the sake
// of efficiency
static void adouble_set_derivatives(redukti_adouble_t *v, size_t num_variables, int order, double value, int index1,
				    double d1, int index2, double d2, double d1d1, double d1d2, double d2d1,
				    double d2d2)
{
	assert(index1 >= 0);
	assert(index1 != index2);
	assert(index1 < (int)num_variables);
	assert(index2 < (int)num_variables);
	assert(order >= 0);
	assert(order <= 2);
	redukti_adouble_init(v, num_variables, order, -1, value);
	if (order && index1 >= 0 && num_variables > 0)
		redukti_adouble_set_derivative1(v, index1, d1);
	if (order && index2 >= 0 && num_variables > 0)
		redukti_adouble_set_derivative1(v, index2, d2);
	if (order == 2 && index1 >= 0 && index2 >= 0 && num_variables > 0) {
		redukti_adouble_set_derivative2(v, index1, index1, d1d1);
		redukti_adouble_set_derivative2(v, index1, index2, d1d2);
		redukti_adouble_set_derivative2(v, index2, index1, d2d1);
		redukti_adouble_set_derivative2(v, index2, index2, d2d2);
	}
}

// discount function using adouble
struct discount_adtype {
	static void discount(FixedRegionAllocator *A, redukti_adouble_t &v, size_t num_vars, int order, int var,
			     const CurveReference *curve, Date date, StatusCode &status)
	{
		auto t = curve->get()->time_from_reference(date);
		double r = curve->get()->zero_rate(t);
		// return adouble2(num_vars, order, var, std::exp(-r*t), -t);
		adouble_exp_init(&v, num_vars, order, var, r, -t);
	}
};

// discount function using double
struct discount_double {
	static void discount(FixedRegionAllocator *A, double &v, size_t num_vars, int order, int var,
			     const CurveReference *curve, Date date, StatusCode &status)
	{
		v = curve->get()->discount(date);
	}
};

// forward function using adouble
struct forward_adtype {
	static void forward(FixedRegionAllocator *A, redukti_adouble_t &v, size_t num_vars, int order,
			    const CurveReference *curve, Date valueDt, int varValueDt, Date maturityDt,
			    int varMaturityDt, double tau, IndexId fixing_key, Date fixingDt,
			    const ValuationContext &ctx, StatusCode &status)
	{
		double fixing = 0.0;
		if (ctx.get_fixing(fixing_key, fixingDt, fixing, status)) {
			redukti_adouble_init(&v, num_vars, order, -1, fixing);
			return;
		}
		if (status != StatusCode::kOk) {
			redukti_adouble_init(&v, num_vars, order, -1, 0.0);
			return;
		}

		// here we implement an optimized version of
		// (df(value date) / df(maturity date) - 1.0) / tau

		auto t1 = curve->get()->time_from_reference(valueDt);
		auto t2 = curve->get()->time_from_reference(maturityDt);
		auto r1 = curve->get()->zero_rate(t1);
		auto r2 = curve->get()->zero_rate(t2);

		t1 = -t1;

		auto x = r2 * t2 + r1 * t1;
		auto z = std::exp(x);
		auto value = (z - 1.0) / tau;

		auto d1 = (t1 * z) / tau;
		auto d2 = (t2 * z) / tau;

		auto d1d1 = d1 * t1;
		auto d1d2 = d1 * t2;
		auto d2d1 = d1d2;
		auto d2d2 = d2 * t2;

		adouble_set_derivatives(&v, num_vars, order, value, varValueDt, d1, varMaturityDt, d2, d1d1, d1d2, d2d1,
					d2d2);
	}
};

// forward function using double
struct forward_double {
	static void forward(FixedRegionAllocator *A, double &v, size_t num_vars, int order, const CurveReference *curve,
			    Date valueDt, int varValueDt, Date maturityDt, int varMaturityDt, double tau,
			    IndexId fixing_key, Date fixingDt, const ValuationContext &ctx, StatusCode &status)
	{
		double fixing = 0.0;
		if (ctx.get_fixing(fixing_key, fixingDt, fixing, status)) {
			v = fixing;
			return;
		}
		if (status != StatusCode::kOk) {
			v = 0.0;
			return;
		}

		auto t1 = curve->get()->time_from_reference(valueDt);
		auto t2 = curve->get()->time_from_reference(maturityDt);
		auto r2 = curve->get()->zero_rate(t2);
		auto r1 = curve->get()->zero_rate(t1);
		v = (std::exp(r2 * t2 - r1 * t1) - 1.0) / tau;
	}
};

// telescopic forward function using double
struct telescope_double {
	static void telescopic_forward(FixedRegionAllocator *A, double &v, size_t num_vars, int order,
				       const CurveReference *curve, Date valueDt, int varValueDt, Date maturityDt,
				       int varMaturityDt, StatusCode &status)
	{
		// this is doing discount(valueDt) / disount(maturityDt)
		auto t1 = curve->get()->time_from_reference(valueDt);
		auto t2 = curve->get()->time_from_reference(maturityDt);
		v = std::exp(curve->get()->zero_rate(t2) * t2 - curve->get()->zero_rate(t1) * t1);
	}
};

struct telescope_adtype {
	static void telescopic_forward(FixedRegionAllocator *A, redukti_adouble_t &v, size_t num_vars, int order,
				       const CurveReference *curve, Date valueDt, int varValueDt, Date maturityDt,
				       int varMaturityDt, StatusCode &status)
	{
		// this is doing discount(valueDt) / disount(maturityDt)
		auto t1 = curve->get()->time_from_reference(valueDt);
		double r1 = curve->get()->zero_rate(t1);
		auto t2 = curve->get()->time_from_reference(maturityDt);
		double r2 = curve->get()->zero_rate(t2);
		// We need some temp memory
		FixedRegionAllocatorGuard g(A);
		auto required_size = redukti_adouble_alloc_size(num_vars, order);
		redukti_adouble_t *temp = (redukti_adouble_t *)A->allocate(required_size);
		redukti_adouble_init(temp, num_vars, order, varValueDt, r1);
		redukti_adouble_init(&v, num_vars, order, varMaturityDt, r2);
		redukti_adouble_scalar_multiply(&v, t2);
		redukti_adouble_add(&v, temp, -t1);
		redukti_adouble_exp(&v, temp);
	}
};

// General case not implemented
template <typename ad_type> struct adtype_traits;

template <> struct adtype_traits<double> {
	typedef discount_double discount_type;
	typedef forward_double forward_type;
	typedef telescope_double telescope_type;
};

template <> struct adtype_traits<redukti_adouble_t> {
	typedef discount_adtype discount_type;
	typedef forward_adtype forward_type;
	typedef telescope_adtype telescope_type;
};

// Compute PV of a simple cashflow
// Only discount sensitivity to payment date
template <typename V>
void pv(FixedRegionAllocator *A, V &accrued_amount, int vars, int order, const NSimpleCashflow *cf,
	const CurveReference *discount_curve, const CurveReference *, const CurveReference *, int,
	const ValuationContext &ctx, StatusCode &status)
{
	typedef typename adtype_traits<V>::discount_type discount_type;
	typedef AutoDiffTraits<V> AD;

	discount_type::discount(A, accrued_amount, vars, order, cf->var, discount_curve, cf->paydate, status);
	AD::scalar_multiply(accrued_amount, cf->amount);
}

// OIS Cashflow pricing
template <typename V>
void pv(FixedRegionAllocator *A, V &rate, int vars, int order, const NOISCashflow *cf,
	const CurveReference *discount_curve, const CurveReference *forward_curve, const CurveReference *, int,
	const ValuationContext &ctx, StatusCode &status)
{
	typedef typename adtype_traits<V>::telescope_type telescope_type;
	typedef typename adtype_traits<V>::discount_type discount_type;
	typedef AutoDiffTraits<V> AD;
	typedef typename AD::ptr_type ptr_type;

	FixedRegionAllocatorGuard G(A);
	StackRegionWithFallbackAllocator<64> buf(A);
	// compound factor is a constant
	if (cf->value_date_var != -1) {
		telescope_type::telescopic_forward(A, rate, vars, order, forward_curve, cf->value_date,
						   cf->value_date_var, cf->maturity_date, cf->maturity_date_var,
						   status);
		AD::scalar_multiply(rate, cf->compound_factor);
	} else {
		AD::init(rate, vars, order, -1, cf->compound_factor);
	}
	// V rate = (compoundFactor - 1.0) / cf->period_tau;
	// TODO round OIS rate
	// rate *= cf->notional;
	// rate *= cf->period_tau;
	AD::scalar_add(rate, -1.0);
	AD::scalar_multiply(rate, cf->notional);

	ptr_type discount_factor = AD::allocate(&buf, vars, order, -1, 0.0);
	ptr_type tempbuf1 = AD::allocate(&buf, vars, order, -1, 0.0); // For multiply
	discount_type::discount(A, AD::reference(discount_factor), vars, order, cf->var, discount_curve, cf->paydate,
				status);
	// rate *= discount(vars, ctx.derivative_order, cf->var, discount_curve,
	//  cf->paydate);
	// return rate;
	AD::multiply(rate, AD::reference(discount_factor), AD::reference(tempbuf1));
}

template <typename V>
void pv(FixedRegionAllocator *A, V &accrued_amount, int vars, int order, const NFloatingCashflow *cf,
	const CurveReference *discount_curve, const CurveReference *forward_curve, const CurveReference *forward_curve2,
	int, const ValuationContext &ctx, StatusCode &status)
{
	typedef typename adtype_traits<V>::discount_type discount_type;
	typedef typename adtype_traits<V>::forward_type forward_type;
	typedef AutoDiffTraits<V> AD;
	typedef typename AD::ptr_type ptr_type;

	FixedRegionAllocatorGuard G(A);
	StackRegionWithFallbackAllocator<64> buf(A);
	ptr_type rate = AD::allocate(&buf, vars, order, -1, 0.0);
	ptr_type rate2 = AD::allocate(&buf, vars, order, -1, 0.0);
	ptr_type accrual_rate = AD::allocate(&buf, vars, order, -1, 0.0);
	ptr_type tempbuf1 = AD::allocate(&buf, vars, order, -1, 0.0); // for multiply
	AD::init(accrued_amount, vars, order, -1, 0.0);
	bool flag = true;
	for (FloatingPeriod *fp = cf->head; fp != nullptr; fp = fp->next) {
		// Calculate first rate
		forward_type::forward(A, AD::reference(rate), vars, order, forward_curve, fp->float1.value_date,
				      fp->float1.value_date_var, fp->float1.maturity_date, fp->float1.maturity_date_var,
				      fp->float1.index_tau, fp->float1.index->id(), fp->float1.fixing_date, ctx,
				      status);
		if (status != StatusCode::kOk)
			return;

		if (fp->float2.index != nullptr) {
			assert(forward_curve2 != nullptr);
			forward_type::forward(
			    A, AD::reference(rate2), vars, order, forward_curve2, fp->float2.value_date,
			    fp->float2.value_date_var, fp->float2.maturity_date, fp->float2.maturity_date_var,
			    fp->float2.index_tau, fp->float2.index->id(), fp->float2.fixing_date, ctx, status);
			if (status != StatusCode::kOk)
				return;

			// TODO must round
			// interpolated_rate = rate1 + (rate2 - rate1) * weight;
			AD::add(AD::reference(rate2), AD::reference(rate), -1.0);
			AD::scalar_multiply(AD::reference(rate2), fp->weight);
			AD::add(AD::reference(rate), AD::reference(rate2), 1.0);
		}

		if (!flag) {
			if (cf->compounding_method != CompoundingMethod::COMPOUNDING_METHOD_NONE) {
				if (cf->compounding_method == CompoundingMethod::FLAT) { // Flat
					// accrual_rate = rate + spread;
					// accrued_amount += fp->notional *
					// fp->period_tau * rate +
					//    accrued_amount * fp->period_tau *
					//    accrual_rate;
					// V accrual_rate(rate);
					AD::assign(AD::reference(accrual_rate), AD::reference(rate));
					// accrual_rate *= fp->period_tau;
					// accrual_rate *= accrued_amount;
					AD::scalar_multiply(AD::reference(accrual_rate), fp->period_tau);
					AD::multiply(AD::reference(accrual_rate), accrued_amount,
						     AD::reference(tempbuf1));

					// rate += fp->spread;
					AD::scalar_add(AD::reference(rate), fp->spread);
					// rate *= (fp->notional *
					// fp->period_tau);
					AD::scalar_multiply(AD::reference(rate), (fp->notional * fp->period_tau));
					// accrued_amount += rate;
					AD::add(accrued_amount, AD::reference(rate), 1.0);
					// accrued_amount += accrual_rate;
					AD::add(accrued_amount, AD::reference(accrual_rate), 1.0);
				} else if (cf->compounding_method == CompoundingMethod::STRAIGHT) { // Straight
					// accrued_amount += fp->notional *
					// fp->period_tau * rate +
					//    accrued_amount * fp->period_tau *
					//    rate;
					// rate += fp->spread;
					AD::scalar_add(AD::reference(rate), fp->spread);
					// accrued_amount *= rate;
					AD::assign(AD::reference(accrual_rate), accrued_amount);
					AD::multiply(AD::reference(accrual_rate), AD::reference(rate),
						     AD::reference(tempbuf1));
					// accrued_amount *= fp->period_tau;
					AD::scalar_multiply(AD::reference(accrual_rate), fp->period_tau);
					// rate *= (fp->notional *
					// fp->period_tau);
					AD::scalar_multiply(AD::reference(rate), (fp->notional * fp->period_tau));
					// accrued_amount += rate;
					AD::add(accrued_amount, AD::reference(rate), 1.0);
					AD::add(accrued_amount, AD::reference(accrual_rate), 1.0);
				}
			} else {
				// accrued_amount += fp->notional *
				// fp->period_tau * rate rate += fp->spread;
				AD::scalar_add(AD::reference(rate), fp->spread);
				// rate *= (fp->notional * fp->period_tau);
				AD::scalar_multiply(AD::reference(rate), (fp->notional * fp->period_tau));
				// accrued_amount += rate;
				AD::add(accrued_amount, AD::reference(rate), 1.0);
			}
		} else {
			flag = false;
			// rate += fp->spread;
			AD::scalar_add(AD::reference(rate), fp->spread);
			// rate *= fp->notional;
			// rate *= fp->period_tau;
			AD::scalar_multiply(AD::reference(rate), (fp->notional * fp->period_tau));
			// accrued_amount = rate;
			AD::assign(accrued_amount, AD::reference(rate));
		}
	}
	discount_type::discount(A, AD::reference(rate2), vars, order, cf->var, discount_curve, cf->paydate, status);
	AD::multiply(accrued_amount, AD::reference(rate2), AD::reference(tempbuf1));
}

template <typename V>
void pv(FixedRegionAllocator *A, V &accrued_amount, int vars, int order, const NFraCashflow *cf,
	const CurveReference *discount_curve, const CurveReference *forward_curve, const CurveReference *forward_curve2,
	int, const ValuationContext &ctx, StatusCode &status)
{
	typedef typename adtype_traits<V>::discount_type discount_type;
	typedef typename adtype_traits<V>::forward_type forward_type;
	typedef AutoDiffTraits<V> AD;
	typedef typename AD::ptr_type ptr_type;

	FixedRegionAllocatorGuard G(A);
	StackRegionWithFallbackAllocator<64> buf(A);
	ptr_type rate = AD::allocate(&buf, vars, order, -1, 0.0);

	FloatingPeriod *fp = cf->floating_period;
	forward_type::forward(A, AD::reference(rate), vars, order, forward_curve, fp->float1.value_date,
			      fp->float1.value_date_var, fp->float1.maturity_date, fp->float1.maturity_date_var,
			      fp->float1.index_tau, fp->float1.index->id(), fp->float1.fixing_date, ctx, status);
	if (status != StatusCode::kOk)
		return;
	if (fp->float2.index != nullptr) {
		assert(forward_curve2 != nullptr);
		forward_type::forward(A, accrued_amount, vars, order, forward_curve2, fp->float2.value_date,
				      fp->float2.value_date_var, fp->float2.maturity_date, fp->float2.maturity_date_var,
				      fp->float2.index_tau, fp->float2.index->id(), fp->float2.fixing_date, ctx,
				      status);
		if (status != StatusCode::kOk)
			return;

		// TODO must round
		// interpolated_rate = rate1 + (rate2 - rate1) * weight;
		AD::add(accrued_amount, AD::reference(rate), -1.0);
		AD::scalar_multiply(accrued_amount, fp->weight);
		AD::add(AD::reference(rate), accrued_amount, 1.0);
	}

	discount_type::discount(A, accrued_amount, vars, order, cf->var, discount_curve, cf->paydate, status);

	ptr_type temp = AD::allocate(&buf, vars, order, -1, cf->fixed_rate);
	ptr_type tempbuf1 = AD::allocate(&buf, vars, order, -1, 0.0); // For multiply
	ptr_type tempbuf2 = AD::allocate(&buf, vars, order, -1, 0.0); // For divide

	// fixedRate - rate
	AD::add(AD::reference(temp), AD::reference(rate), -1.0);
	// (fixedRate - rate) * period * doubleDiscountingFactor * notional
	AD::scalar_multiply(AD::reference(temp), fp->period_tau * cf->double_discounting_factor * fp->notional);
	// DF * (fixedRate - rate) * period * doubleDiscountingFactor * notional
	AD::multiply(accrued_amount, AD::reference(temp), AD::reference(tempbuf1));
	// rate * period
	AD::scalar_multiply(AD::reference(rate), fp->period_tau);
	// 1 + rate * period
	AD::scalar_add(AD::reference(rate), 1.0);
	// npv = doubleDiscountingFactor * DF * notional * (fixedRate - rate) *
	// period / (1 + rate * period);
	AD::divide(accrued_amount, AD::reference(rate), AD::reference(tempbuf1), AD::reference(tempbuf2));
}

template <typename V, typename CF>
double pv_specific(FixedRegionAllocator *A, const CF *cf, VariableResolvers *resolvers,
		   const CurveReference *discount_curve, const CurveReference *forward_curve,
		   const CurveReference *forward_curve2, double sign, const ValuationContext &ctx,
		   Sensitivities &container, StatusCode &status, const VariableResolver *discount_resolver = nullptr,
		   const VariableResolver *forward1_resolver = nullptr,
		   const VariableResolver *forward2_resolver = nullptr)
{
	typedef AutoDiffTraits<V> AD;
	typedef typename AD::ptr_type ptr_type;
	if (cf->paydate < ctx.payment_cutoff_date())
		return 0.0;
	int vars = resolvers ? resolvers->num_vars() : 0;
	int order = ctx.derivative_order();
	FixedRegionAllocatorGuard G(A);
	StackRegionWithFallbackAllocator<64> buf(A);
	ptr_type pv_amount = AD::allocate(&buf, vars, order, -1, 0.0);
	pv<V>(A, AD::reference(pv_amount), vars, order, cf, discount_curve, forward_curve, forward_curve2, 0, ctx,
	      status);
	if (status != StatusCode::kOk)
		return 0.0;
	AD::scalar_multiply(AD::reference(pv_amount), sign);
	double value = AD::value(pv_amount);
	ValuationContextHelper<V> helper;
	helper.aggregateSensitivity(A, discount_resolver, discount_curve, forward1_resolver, forward_curve,
				    forward2_resolver, forward_curve2, &container, AD::reference(pv_amount), status);
	return value;
}

template <typename V, typename CF>
double pv_specific(FixedRegionAllocator *A, const CF *cf, VariableResolvers *resolvers,
		   const CurveProvider *mapping_provider, double sign, const ValuationContext &ctx,
		   Sensitivities &container, StatusCode &status)
{
	assert(resolvers);
	if (resolvers->num_resolvers() > 3) {
		status = StatusCode::kVAL_DetectedTooManyRiskFactorsInCashflow;
		return 0.0;
	}
	PricingCurve discountcurveid = resolvers->get_curve_of_type(PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT, 0);
	PricingCurve forwardcurveid = resolvers->get_curve_of_type(PricingCurveType::PRICING_CURVE_TYPE_FORWARD, 0);
	PricingCurve forwardcurveid2 = resolvers->get_curve_of_type(PricingCurveType::PRICING_CURVE_TYPE_FORWARD, 1);
	if (forwardcurveid.id() == 0)
		forwardcurveid = discountcurveid;
	if (forwardcurveid2.id() == 0)
		forwardcurveid2 = forwardcurveid;
	const CurveReference *discount_curve = mapping_provider->get_curve(discountcurveid);
	const CurveReference *forward_curve = mapping_provider->get_curve(forwardcurveid);
	const CurveReference *forward_curve2 = mapping_provider->get_curve(forwardcurveid2);
	if (discount_curve == nullptr || forward_curve == nullptr || forward_curve2 == nullptr) {
		status = StatusCode::kVAL_MappedCurveNotFound;
		return 0.0;
	}
	auto resolverdiscount = resolvers->get_resolver(discountcurveid);
	auto resolverforward1 = resolvers->get_resolver(forwardcurveid);
	auto resolverforward2 = resolvers->get_resolver(forwardcurveid2);
	//  if (redukti_debug) {
	//    fprintf(stdout, "INFO: Using %s for discounting, %s and %s for
	//    forward rates\n",
	//      curve_id_to_string(discount_curve->get()->id()).c_str(),
	//      curve_id_to_string(forward_curve->get()->id()).c_str(),
	//      curve_id_to_string(forward_curve2->get()->id()).c_str());
	//  }
	return pv_specific<V>(A, cf, resolvers, discount_curve, forward_curve, forward_curve2, sign, ctx, container,
			      status, &resolverdiscount, &resolverforward1, &resolverforward2);
}

template <typename V>
double pv(FixedRegionAllocator *A, const Cashflow *flow, const CurveReference *discount_curve,
	  const CurveReference *forward_curve, const CurveReference *forward_curve2, double sign,
	  const ValuationContext &ctx, Sensitivities &container, StatusCode &status)
{
	double value = 0.0;
	switch (flow->type) {
	case CashflowType::CF_SIMPLE: {
		auto cf = flow->simple;
		value = pv_specific<V>(A, cf, flow->resolvers, discount_curve, forward_curve, forward_curve2, sign, ctx,
				       container, status);
		break;
	}
	case CashflowType::CF_FLOATING: {
		auto cf = flow->floating;
		value = pv_specific<V>(A, cf, flow->resolvers, discount_curve, forward_curve, forward_curve2, sign, ctx,
				       container, status);
		break;
	}
	case CashflowType::CF_OIS: {
		auto cf = flow->ois;
		value = pv_specific<V>(A, cf, flow->resolvers, discount_curve, forward_curve, forward_curve2, sign, ctx,
				       container, status);
		break;
	}
	case CashflowType::CF_FRA: {
		auto cf = flow->fra;
		value = pv_specific<V>(A, cf, flow->resolvers, discount_curve, forward_curve, forward_curve2, sign, ctx,
				       container, status);
		break;
	}
	case CashflowType::CF_NONE:
		break;
	}
	return value;
}

template <typename V>
double pv(FixedRegionAllocator *A, const Cashflow *flow, const CurveProvider *mapping_provider, double sign,
	  const ValuationContext &ctx, Sensitivities &container, StatusCode &status)
{
	double value = 0.0;
	switch (flow->type) {
	case CashflowType::CF_SIMPLE: {
		auto cf = flow->simple;
		value = pv_specific<V>(A, cf, flow->resolvers, mapping_provider, sign, ctx, container, status);
		break;
	}
	case CashflowType::CF_FLOATING: {
		auto cf = flow->floating;
		value = pv_specific<V>(A, cf, flow->resolvers, mapping_provider, sign, ctx, container, status);
		break;
	}
	case CashflowType::CF_OIS: {
		auto cf = flow->ois;
		value = pv_specific<V>(A, cf, flow->resolvers, mapping_provider, sign, ctx, container, status);
		break;
	}
	case CashflowType::CF_FRA: {
		auto cf = flow->fra;
		value = pv_specific<V>(A, cf, flow->resolvers, mapping_provider, sign, ctx, container, status);
		break;
	}
	case CashflowType::CF_NONE:
		break;
	}
	return value;
}

template <typename V>
double pv(FixedRegionAllocator *A, const CashflowStream *leg, const CurveReference *discount_curve,
	  const CurveReference *forward_curve, const CurveReference *forward_curve2, double sign,
	  const ValuationContext &ctx, Sensitivities &container, StatusCode &status)
{
	double value = 0.0;
	for (Cashflow *flow = leg->first; flow != nullptr; flow = flow->next) {
		value += pv<V>(A, flow, discount_curve, forward_curve, forward_curve2, sign, ctx, container, status);
		if (status != StatusCode::kOk)
			break;
	}
	return value;
}

template <typename V>
double pv(FixedRegionAllocator *A, const CashflowStream *leg, const CurveProvider *mapping_provider, double sign,
	  const ValuationContext &ctx, Sensitivities &container, StatusCode &status)
{
	double value = 0.0;
	for (Cashflow *flow = leg->first; flow != nullptr; flow = flow->next) {
		value += pv<V>(A, flow, mapping_provider, sign, ctx, container, status);
		if (status != StatusCode::kOk)
			break;
	}
	return value;
}

// FIXME below is problematic as it cannot handle more than 2 curve references
template <typename V>
double pv(FixedRegionAllocator *A, const Cashflows *cf, const CurveReference *discount_curve,
	  const CurveReference *forward_curve1, const CurveReference *forward_curve2, const ValuationContext &ctx,
	  Sensitivities &container, StatusCode &status)
{
	status = StatusCode::kOk;
	if (discount_curve == nullptr || forward_curve1 == nullptr || forward_curve2 == nullptr) {
		status = StatusCode::kVAL_MappedCurveNotFound;
		return 0.0;
	}
	if (cf->n_streams() == 2) {
		// FIXME there is an assumption below that first leg is based on
		// first forward curve and second leg on second forward curve
		double leg1 = pv<V>(A, cf->first_, discount_curve, forward_curve1, forward_curve1, cf->first_->factor,
				    ctx, container, status);
		if (status != StatusCode::kOk)
			return 0.0;
		double leg2 = pv<V>(A, cf->last_, discount_curve, forward_curve2, forward_curve2, cf->last_->factor,
				    ctx, container, status);
		if (status != StatusCode::kOk)
			return 0.0;
		leg1 += leg2;
		return leg1;
	} else if (cf->n_streams() == 1) {
		return pv<V>(A, cf->first_, discount_curve, forward_curve1, forward_curve1, cf->first_->factor, ctx,
			     container, status);
	}
	status = StatusCode::kVAL_UnsupportedCashflowType;
	return 0;
}

// FIXME below is problematic as it cannot handle more than 2 curve references
template <typename V>
double pv(FixedRegionAllocator *A, const Cashflows *cf, const CurveProvider *mapping_provider,
	  const ValuationContext &ctx, Sensitivities &container, StatusCode &status)
{
	status = StatusCode::kOk;
	double v = 0.0;
	for (CashflowStream *stream = cf->first_; stream != nullptr; stream = stream->next) {
		v += pv<V>(A, stream, mapping_provider, stream->factor, ctx, container, status);
		if (status != StatusCode::kOk)
			break;
	}
	return v;
}

template <typename V> struct PresentValueFunction {
	double present_value(FixedRegionAllocator *A, const Cashflows *cf, const CurveReference *discount_curve,
			     const CurveReference *forward_curve1, const CurveReference *forward_curve2,
			     ValuationContext &ctx, Sensitivities &sens, StatusCode &status)
	{
		FixedRegionAllocatorGuard guard(A);
		return pv<V>(A, cf, discount_curve, forward_curve1, forward_curve2, ctx, sens, status);
	}
};

double compute_present_value(FixedRegionAllocator *A, const Cashflows *flows, const CurveReference *discount_curve,
			     const CurveReference *forward_curve1, const CurveReference *forward_curve2,
			     const ValuationContext &ctx, Sensitivities &sensitivities, StatusCode &status)
{
	FixedRegionAllocatorGuard guard(A);
	if (ctx.derivative_order() == 0) {
		return pv<double>(A, flows, discount_curve, forward_curve1, forward_curve2, ctx, sensitivities, status);
	} else {
		return pv<redukti_adouble_t>(A, flows, discount_curve, forward_curve1, forward_curve2, ctx,
					     sensitivities, status);
	}
}

double compute_present_value(FixedRegionAllocator *A, const ValuationContext &ctx, const Cashflows *flows,
			     const CurveProvider *mapping_provider, Sensitivities &sensitivities, StatusCode &status)
{
	FixedRegionAllocatorGuard guard(A);
	if (ctx.derivative_order() == 0) {
		double v = pv<double>(A, flows, mapping_provider, ctx, sensitivities, status);
		return v;
	} else {
		return pv<redukti_adouble_t>(A, flows, mapping_provider, ctx, sensitivities, status);
	}
}

inline bool is_zero(double v)
{
	// return ::fabs(v) < 1e-13;
	return v == 0.0;
}

class SimpleCurveReference : public CurveReference
{
	private:
	std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve_;

	public:
	SimpleCurveReference(std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve) noexcept : curve_(std::move(curve))
	{
	}
	YieldCurve *get() const noexcept override final { return curve_.get(); }
};

void compute_sensitivity_numerically(FixedRegionAllocator *A, const Cashflows *flows,
				     const CurveReference *discount_curve, const CurveReference *forward_curve1,
				     const CurveReference *forward_curve2, Sensitivities *container, StatusCode &status,
				     double h)
{
	if (discount_curve == nullptr || forward_curve1 == nullptr || forward_curve2 == nullptr) {
		status = StatusCode::kVAL_MappedCurveNotFound;
		return;
	}

	ValuationContextImpl ctx(0, 0, nullptr);
	PresentValueFunction<double> pv_function;

	const double up = h;
	const double down = -h;

	// There are three possibilities
	// 1) Only one curve is involved - i.e. discount_curve == forward_curve1
	// == forward_curve2 2) discount_curve != forward_curve1 but
	// forward_curve1 == forward_curve2 (vanilla swap) 3) discount_curve !=
	// forward_curve1 and forward_curve1 != forward_curve2 (basis swap)

	// In the first case, we need to compute gamma for DD
	// In the second case, we need to compute gamma for F1F1, F1D, DD
	//   F1F1 F1D
	//   -DF1 DD
	// In the third case, we need to compute gamma for
	// F1F1, F1F2, F1D, F2F2, F2D, DD
	//   F1F1  F1F2  F1D
	//   -F2F1 F2F2  F2D
	//   -DF1  -DF2  DD

	auto discount_delta = container->first_order_sensitivities(discount_curve->get());
	auto discount_discount_gamma =
	    container->second_order_sensitivities(discount_curve->get(), discount_curve->get());

	DynamicRegionAllocator A2(32 * 1024);

	auto up_discount_curves = discount_curve->get()->get_bumped_curves(&A2, up);
	auto down_discount_curves = discount_curve->get()->get_bumped_curves(&A2, down);

	// discount delta
	for (int i = 1; i <= discount_curve->get()->last_pillar(); i++) {
		CurveWrapper up_curve(up_discount_curves[i - 1].get());
		CurveWrapper down_curve(down_discount_curves[i - 1].get());

		double npv_up = pv_function.present_value(A, flows, &up_curve, forward_curve1, forward_curve2, ctx,
							  *container, status);
		if (status != StatusCode::kOk)
			return;
		double npv_down = pv_function.present_value(A, flows, &down_curve, forward_curve1, forward_curve2, ctx,
							    *container, status);
		if (status != StatusCode::kOk)
			return;
		double tmp = npv_up - npv_down;
		if (!is_zero(tmp)) {
			double delta = tmp / (2 * h);
			discount_delta->at(i) += delta;
		}
	}
	// discount discount gamma
	for (int i = 1; i <= discount_curve->get()->last_pillar(); i++) {
		for (int j = 1; j <= discount_curve->get()->last_pillar(); j++) {
			DynamicRegionAllocator Atemp(32 * 1024);
			// in = i negative, ip = i positive
			// jn = j negative, jp = j positive
			SimpleCurveReference in_jn_curve(
			    down_discount_curves[i - 1]->get_bumped_curve(&Atemp, j, down));
			SimpleCurveReference ip_jn_curve(up_discount_curves[i - 1]->get_bumped_curve(&Atemp, j, down));
			SimpleCurveReference in_jp_curve(down_discount_curves[i - 1]->get_bumped_curve(&Atemp, j, up));
			SimpleCurveReference ip_jp_curve(up_discount_curves[i - 1]->get_bumped_curve(&Atemp, j, up));

			double npv_in_jn = pv_function.present_value(A, flows, &in_jn_curve, forward_curve1,
								     forward_curve2, ctx, *container, status);
			if (status != StatusCode::kOk)
				return;
			double npv_ip_jn = pv_function.present_value(A, flows, &ip_jn_curve, forward_curve1,
								     forward_curve2, ctx, *container, status);
			if (status != StatusCode::kOk)
				return;
			double npv_in_jp = pv_function.present_value(A, flows, &in_jp_curve, forward_curve1,
								     forward_curve2, ctx, *container, status);
			if (status != StatusCode::kOk)
				return;
			double npv_ip_jp = pv_function.present_value(A, flows, &ip_jp_curve, forward_curve1,
								     forward_curve2, ctx, *container, status);
			if (status != StatusCode::kOk)
				return;

			double tmp = npv_in_jn - npv_ip_jn - npv_in_jp + npv_ip_jp;
			if (!is_zero(tmp)) {
				double gamma = tmp / (4.0 * h * h);
				discount_discount_gamma->at(i, j) += gamma;
			}
		}
	}

	if (discount_curve->get()->id() != forward_curve1->get()->id()) {
		// second case
		auto forward1_delta = container->first_order_sensitivities(forward_curve1->get());
		auto forward1_forward1_gamma =
		    container->second_order_sensitivities(forward_curve1->get(), forward_curve1->get());
		auto forward1_discount_gamma =
		    container->second_order_sensitivities(forward_curve1->get(), discount_curve->get());

		DynamicRegionAllocator A2(32 * 1024);
		auto up_forward_curves = forward_curve1->get()->get_bumped_curves(&A2, up);
		auto down_forward_curves = forward_curve1->get()->get_bumped_curves(&A2, down);

		// forward1 delta
		for (int i = 1; i <= forward_curve1->get()->last_pillar(); i++) {
			CurveWrapper up_curve(up_forward_curves[i - 1].get());
			CurveWrapper down_curve(down_forward_curves[i - 1].get());

			double npv_up = pv_function.present_value(A, flows, discount_curve, &up_curve, forward_curve2,
								  ctx, *container, status);
			if (status != StatusCode::kOk)
				return;
			double npv_down = pv_function.present_value(A, flows, discount_curve, &down_curve,
								    forward_curve2, ctx, *container, status);
			if (status != StatusCode::kOk)
				return;
			double tmp = npv_up - npv_down;
			if (!is_zero(tmp)) {
				double delta = tmp / (2 * h);
				forward1_delta->at(i) += delta;
			}
		}

		// forward1 forward1 gamma
		for (int i = 1; i <= forward_curve1->get()->last_pillar(); i++) {
			for (int j = 1; j <= forward_curve1->get()->last_pillar(); j++) {
				DynamicRegionAllocator Atemp;
				SimpleCurveReference in_jn_curve(
				    down_forward_curves[i - 1]->get_bumped_curve(&Atemp, j, down));
				SimpleCurveReference ip_jn_curve(
				    up_forward_curves[i - 1]->get_bumped_curve(&Atemp, j, down));
				SimpleCurveReference in_jp_curve(
				    down_forward_curves[i - 1]->get_bumped_curve(&Atemp, j, up));
				SimpleCurveReference ip_jp_curve(
				    up_forward_curves[i - 1]->get_bumped_curve(&Atemp, j, up));
				double npv_in_jn = pv_function.present_value(A, flows, discount_curve, &in_jn_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double npv_ip_jn = pv_function.present_value(A, flows, discount_curve, &ip_jn_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double npv_in_jp = pv_function.present_value(A, flows, discount_curve, &in_jp_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double npv_ip_jp = pv_function.present_value(A, flows, discount_curve, &ip_jp_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double tmp = npv_in_jn - npv_ip_jn - npv_in_jp + npv_ip_jp;
				if (!is_zero(tmp)) {
					double gamma = tmp / (4.0 * h * h);
					forward1_forward1_gamma->at(i, j) += gamma;
				}
			}
		}

		// forward1 discount delta
		for (int i = 1; i <= forward_curve1->get()->last_pillar(); i++) {
			for (int j = 1; j <= discount_curve->get()->last_pillar(); j++) {
				CurveWrapper f_up_curve(up_forward_curves[i - 1].get());
				CurveWrapper f_dn_curve(down_forward_curves[i - 1].get());
				CurveWrapper d_up_curve(up_discount_curves[j - 1].get());
				CurveWrapper d_dn_curve(down_discount_curves[j - 1].get());

				double npv_in_jn = pv_function.present_value(A, flows, &d_dn_curve, &f_dn_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double npv_ip_jn = pv_function.present_value(A, flows, &d_up_curve, &f_dn_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double npv_in_jp = pv_function.present_value(A, flows, &d_dn_curve, &f_up_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double npv_ip_jp = pv_function.present_value(A, flows, &d_up_curve, &f_up_curve,
									     forward_curve2, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double tmp = npv_in_jn - npv_ip_jn - npv_in_jp + npv_ip_jp;
				if (!is_zero(tmp)) {
					double gamma = tmp / (4.0 * h * h);
					forward1_discount_gamma->at(i, j) += gamma;
				}
			}
		}

		if (discount_curve->get()->id() != forward_curve2->get()->id() &&
		    forward_curve1->get()->id() != forward_curve2->get()->id()) {
			// third case

			auto forward2_delta = container->first_order_sensitivities(forward_curve2->get());
			auto forward1_forward2_gamma =
			    container->second_order_sensitivities(forward_curve1->get(), forward_curve2->get());
			auto forward2_forward2_gamma =
			    container->second_order_sensitivities(forward_curve2->get(), forward_curve2->get());
			auto forward2_discount_gamma =
			    container->second_order_sensitivities(forward_curve2->get(), discount_curve->get());

			DynamicRegionAllocator A3;
			auto up_forward2_curves = forward_curve2->get()->get_bumped_curves(&A3, up);
			auto down_forward2_curves = forward_curve2->get()->get_bumped_curves(&A3, down);

			// forward2 delta
			for (int i = 1; i <= forward_curve2->get()->last_pillar(); i++) {
				CurveWrapper up_curve(up_forward2_curves[i - 1].get());
				CurveWrapper down_curve(down_forward2_curves[i - 1].get());

				double npv_up = pv_function.present_value(A, flows, discount_curve, forward_curve1,
									  &up_curve, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double npv_down = pv_function.present_value(A, flows, discount_curve, forward_curve1,
									    &down_curve, ctx, *container, status);
				if (status != StatusCode::kOk)
					return;
				double tmp = npv_up - npv_down;
				if (!is_zero(tmp)) {
					double delta = tmp / (2 * h);
					forward2_delta->at(i) += delta;
				}
			}

			// forward2 forward2 gamma
			for (int i = 1; i <= forward_curve2->get()->last_pillar(); i++) {
				for (int j = 1; j <= forward_curve2->get()->last_pillar(); j++) {
					DynamicRegionAllocator Atemp(32 * 1024);
					SimpleCurveReference in_jn_curve(
					    down_forward2_curves[i - 1]->get_bumped_curve(&Atemp, j, down));
					SimpleCurveReference ip_jn_curve(
					    up_forward2_curves[i - 1]->get_bumped_curve(&Atemp, j, down));
					SimpleCurveReference in_jp_curve(
					    down_forward2_curves[i - 1]->get_bumped_curve(&Atemp, j, up));
					SimpleCurveReference ip_jp_curve(
					    up_forward2_curves[i - 1]->get_bumped_curve(&Atemp, j, up));

					double npv_in_jn =
					    pv_function.present_value(A, flows, discount_curve, forward_curve1,
								      &in_jn_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_ip_jn =
					    pv_function.present_value(A, flows, discount_curve, forward_curve1,
								      &ip_jn_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_in_jp =
					    pv_function.present_value(A, flows, discount_curve, forward_curve1,
								      &in_jp_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_ip_jp =
					    pv_function.present_value(A, flows, discount_curve, forward_curve1,
								      &ip_jp_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double tmp = npv_in_jn - npv_ip_jn - npv_in_jp + npv_ip_jp;
					if (!is_zero(tmp)) {
						double gamma = tmp / (4.0 * h * h);
						forward2_forward2_gamma->at(i, j) += gamma;
					}
				}
			}

			// forward2 discount gamma
			for (int i = 1; i <= forward_curve2->get()->last_pillar(); i++) {
				for (int j = 1; j <= discount_curve->get()->last_pillar(); j++) {
					CurveWrapper f_up_curve(up_forward2_curves[i - 1].get());
					CurveWrapper f_dn_curve(down_forward2_curves[i - 1].get());
					CurveWrapper d_up_curve(up_discount_curves[j - 1].get());
					CurveWrapper d_dn_curve(down_discount_curves[j - 1].get());

					double npv_in_jn =
					    pv_function.present_value(A, flows, &d_dn_curve, forward_curve1,
								      &f_dn_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_ip_jn =
					    pv_function.present_value(A, flows, &d_up_curve, forward_curve1,
								      &f_dn_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_in_jp =
					    pv_function.present_value(A, flows, &d_dn_curve, forward_curve1,
								      &f_up_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_ip_jp =
					    pv_function.present_value(A, flows, &d_up_curve, forward_curve1,
								      &f_up_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double tmp = npv_in_jn - npv_ip_jn - npv_in_jp + npv_ip_jp;
					if (!is_zero(tmp)) {
						double gamma = tmp / (4.0 * h * h);
						forward2_discount_gamma->at(i, j) += gamma;
					}
				}
			}

			// forward1 forward2 gamma
			for (int i = 1; i <= forward_curve1->get()->last_pillar(); i++) {
				for (int j = 1; j <= forward_curve2->get()->last_pillar(); j++) {
					CurveWrapper f1_up_curve(up_forward_curves[i - 1].get());
					CurveWrapper f1_dn_curve(down_forward_curves[i - 1].get());
					CurveWrapper f2_up_curve(up_forward2_curves[i - 1].get());
					CurveWrapper f2_dn_curve(down_forward2_curves[i - 1].get());

					double npv_in_jn =
					    pv_function.present_value(A, flows, discount_curve, &f1_dn_curve,
								      &f2_dn_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_ip_jn =
					    pv_function.present_value(A, flows, discount_curve, &f1_up_curve,
								      &f2_dn_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_in_jp =
					    pv_function.present_value(A, flows, discount_curve, &f1_dn_curve,
								      &f2_up_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double npv_ip_jp =
					    pv_function.present_value(A, flows, discount_curve, &f1_up_curve,
								      &f2_up_curve, ctx, *container, status);
					if (status != StatusCode::kOk)
						return;
					double tmp = npv_in_jn - npv_ip_jn - npv_in_jp + npv_ip_jp;
					if (!is_zero(tmp)) {
						double gamma = tmp / (4.0 * h * h);
						forward1_forward2_gamma->at(i, j) += gamma;
					}
				}
			}
		}
	}
}

int test_pricing()
{
	// CODE not released
	printf("CashflowPricing tests skipped\n");
	return 0;
}

} // namespace redukti
