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
#include <converters.h>
#include <curve.h>
#include <curve.pb.h>

#include <assert.h>
#include <cmath>
#include <cstring>

#include <algorithm>

namespace redukti
{

// ccy - 5 bits (max value 31)
// index family - 6 bits (max value 63)
// tenor - 7 bits (max value 127)
// as of date - 24 bits
// cycle - 7 bits (max value 127)
// qual - 1 bit
// scenario - 12 bits (max value 4095)
// curvetype - 2 bits

CurveId make_curve_id(PricingCurveType type, Currency ccy, IndexFamily index_family, Tenor tenor, Date as_of_date,
		      short int cycle, MarketDataQualifier qual, short int scenario)
{
	assert(cycle >= 0 && cycle <= 0x7F);
	assert(scenario >= 0 && scenario <= 0xFFF);
	assert(qual >= 0 && qual <= 1);
	assert(tenor >= 0 && tenor <= 0x7F);
	assert(ccy >= 0 && ccy <= Currency::PLN && ccy <= 0x1F);
	assert(index_family >= 0 && index_family <= IndexFamily::TIIE && index_family <= 0x3F);
	assert(type >= 0 && type <= 3);
	CurveId id = ((uint64_t)ccy & 0x1F) | (((uint64_t)index_family & 0x3F) << 5) |
		     (((uint64_t)tenor & 0x7F) << 11) | (((uint64_t)as_of_date & 0xFFFFFF) << 18) |
		     (((uint64_t)cycle & 0x7F) << 42) | (((uint64_t)qual & 0x1) << 49) |
		     (((uint64_t)scenario & 0xFFF) << 50) | (((uint64_t)type & 0x3) << 62);
	return id;
}
bool curve_id_components(CurveId id, PricingCurveType &type, Currency &ccy, IndexFamily &index_family, Tenor &tenor,
			 Date &as_of_date, short int &cycle, MarketDataQualifier &qual, short int &scenario)
{
	ccy = (Currency)(id & 0x1F);
	index_family = (IndexFamily)((id >> 5) & 0x3F);
	tenor = (Tenor)((id >> 11) & 0x7F);
	as_of_date = (Date)((id >> 18) & 0xFFFFFF);
	cycle = (short int)((id >> 42) & 0x7F);
	qual = (MarketDataQualifier)((id >> 49) & 0x1);
	scenario = (short int)((id >> 50) & 0xFFF);
	type = (PricingCurveType)((id >> 62) & 0x3);
	if (ccy < 0 || ccy > Currency::PLN || index_family < 0 || index_family > IndexFamily::TIIE || tenor < 0 ||
	    tenor > TENOR_1T || !is_valid_date(as_of_date) || type < 0 ||
	    type > PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT)
		return false;
	return true;
}
std::string curve_id_to_string(CurveId id)
{
	Currency ccy;
	IndexFamily index_family;
	Tenor tenor;
	Date as_of_date;
	short int cycle;
	MarketDataQualifier qual;
	short int scenario;
	PricingCurveType type;
	if (!curve_id_components(id, type, ccy, index_family, tenor, as_of_date, cycle, qual, scenario))
		return std::string("Bad curveid");
	char buf[100];
	auto ymd = date_components(as_of_date);
	int d = ymd.d, m = ymd.m, y = ymd.y;
	auto converter = get_default_converter();
	snprintf(buf, sizeof buf, "Curve[%s-%s-%s-%s-%d%02d%02d-%s-%d-%d]",
		 type == PRICING_CURVE_TYPE_DISCOUNT ? "Discount"
						     : (type == PRICING_CURVE_TYPE_FORWARD ? "Forward" : ""),
		 converter->currency_to_string(ccy), converter->index_family_to_string(index_family),
		 converter->tenor_to_string(tenor).c_str(), y, m, d,
		 qual == MarketDataQualifier::MDQ_CLOSING ? "Closing" : "Normal", cycle, scenario);
	return std::string(buf);
}

std::unique_ptr<YieldCurve, Deleter<YieldCurve>> make_curve(Date as_of_date, const IRCurveDefinition *defn,
							    const ZeroCurve &curve, int deriv_order,
							    PricingCurveType type, MarketDataQualifier mdq,
							    short int cycle, short int scenario)
{
	auto n = curve.maturities_size();
	int reqsize = n * sizeof(double) * 2;
	FixedRegionAllocator *tempalloc = get_threadspecific_allocators()->tempspace_allocator;
	FixedRegionAllocatorGuard guard(tempalloc);
	StackRegionWithFallbackAllocator<1024> buf(tempalloc);
	Date *maturities = (Date *)buf.safe_allocate(sizeof(Date) * n);
	double *values = (double *)buf.safe_allocate(sizeof(double) * n);
	for (int i = 0; i < n; i++) {
		maturities[i] = curve.maturities(i);
		values[i] = curve.values(i);
	}
	if (defn->interpolated_on() == IRRateType::DISCOUNT_FACTOR) {
		// Convert the values to discount factors as we need to
		// interpolate on discount factors
		// FIXME the day count fraction ought to be a parameter in curve
		// definition
		auto fraction = get_day_fraction(DayCountFraction::ACT_365_FIXED);
		for (int i = 0; i < n; i++) {
			double t = fraction->year_fraction(as_of_date, maturities[i]);
			values[i] = std::exp(-values[i] * t);
		}
	}
	CurveId curveId = make_curve_id(type, defn->currency(), defn->index_family(), defn->tenor(), as_of_date, cycle,
					mdq, scenario);
	return make_curve(&GlobalAllocator, curveId, as_of_date, maturities, values, n, defn->interpolator_type(),
			  defn->interpolated_on(), deriv_order);
}
/**
 * Curve that takes zero rates as input
 *
 */
template <typename Derived> class InterestRateCurveImpl : public YieldCurve
{
	public:
	InterestRateCurveImpl(CurveId id, Date as_of_date, Date maturities[], double values[], size_t n,
			      Allocator *alloc, InterpolatorType interpolator,
			      DayCountFraction fraction = DayCountFraction::ACT_365_FIXED,
			      IRRateType rateType = IRRateType::ZERO_RATE, int deriv_order = 0) noexcept;
	InterestRateCurveImpl(CurveId id, Date as_of_date, Date maturities[], double values[], size_t n,
			      Allocator *alloc, InterpolatorType interpolator, bool special,
			      DayCountFraction fraction = DayCountFraction::ACT_365_FIXED,
			      IRRateType rateType = IRRateType::ZERO_RATE, int deriv_order = 0) noexcept;
	virtual ~InterestRateCurveImpl() noexcept;
	virtual double discount(double time) const noexcept final;
	virtual double forward(double time) const noexcept final;
	virtual double discount(Date d) const noexcept final;
	virtual double zero_rate(Date d) const noexcept final;
	virtual double zero_rate(double t) const noexcept final;
	virtual double forward_rate(Date d1, Date d2) const noexcept final;
	virtual double forward_rate(double t1, double t2) const noexcept final;
	virtual std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
	get_sensitivities(double x, FixedRegionAllocator *A) const noexcept final;
	virtual const DayFraction &day_fraction() const noexcept final { return *day_fraction_; }
	Date as_of_date() const noexcept override final { return as_of_date_; }
	int last_pillar() const noexcept override final { return n_ - 1; }
	Date maturity_date(int pillar) const override final
	{
		if (pillar < 0 || pillar > last_pillar())
			return 0;
		return maturities_[pillar];
	}
	double value(int pillar) const noexcept override final
	{
		if (pillar >= 1 && pillar < (int)n_)
			return values_[pillar];
		assert(pillar >= 1 && pillar < (int)n_);
		return values_[1];
	}
	double maturity_time(int pillar) const noexcept override final
	{
		if (pillar >= 1 && pillar < (int)n_)
			return times_[pillar];
		assert(pillar >= 1 && pillar < (int)n_);
		return times_[1];
	}
	virtual Date last_maturity() const noexcept final;
	virtual void update_rates(const double *rates, size_t n) noexcept;
	virtual std::vector<std::unique_ptr<YieldCurve, Deleter<YieldCurve>>>
	get_bumped_curves(Allocator *A, double h = 0.00001) const noexcept;
	virtual std::unique_ptr<YieldCurve, Deleter<YieldCurve>> get_bumped_curve(Allocator *A, int pillar,
										  double h = 0.00001) const noexcept;
	virtual void dump(FILE *fp = stderr) const noexcept;
	virtual bool is_valid() const noexcept final { return valid_; }
	virtual InterpolatorType interpolator_type() const noexcept final { return interpolator_->type(); }

	protected:
	size_t n_;
	Date *maturities_;
	double *values_;
	Date as_of_date_;
	double *times_;
	const DayFraction *day_fraction_;
	std::unique_ptr<Interpolator, Deleter<Interpolator>> interpolator_;
	Allocator *alloc_;
	IRRateType valueType_;
	int deriv_order_;
	bool valid_;

	private:
	InterestRateCurveImpl(const InterestRateCurveImpl &) = delete;
	InterestRateCurveImpl &operator=(const InterestRateCurveImpl &) = delete;

	private:
	IRRateType apply_values(InterpolatorType itype, const double *rates, size_t n) noexcept;
};

Curve::Curve(CurveId id) noexcept : curve_id_(id) {}

template <typename Derived>
InterestRateCurveImpl<Derived>::InterestRateCurveImpl(CurveId id, Date as_of_date, Date maturities[], double values[],
						      size_t n, Allocator *alloc, InterpolatorType itype,
						      DayCountFraction fraction, IRRateType rateType,
						      int deriv_order) noexcept
    : YieldCurve(id), as_of_date_(as_of_date), alloc_(alloc), valueType_(rateType), deriv_order_(deriv_order)
{
	// FIXME will not work for bus_252
	day_fraction_ = get_day_fraction(fraction);
	assert(day_fraction_ != nullptr);
	// We need an extra point for the zeroeth position
	maturities_ = new (*alloc) Date[n + 1];
	maturities_[0] = as_of_date;
	std::copy(maturities, maturities + n, maturities_ + 1);
	times_ = new (*alloc) double[n + 1];
	times_[0] = 0;
	for (unsigned int i = 1; i < n + 1; i++) {
		times_[i] = day_fraction_->year_fraction(maturities_[0], maturities_[i]);
	}
	n_ = n + 1;
	values_ = new (*alloc) double[n + 1];
	valueType_ = rateType = apply_values(itype, values, n);
	InterpolationOptions data;
	data.differentiation_order = deriv_order;
	if (itype == InterpolatorType::MONOTONE_CONVEX) {
		data.monotoneconvex_inputs_are_forwards = (rateType == IRRateType::FORWARD_RATE);
	}
	interpolator_ = make_interpolator(itype, times_, values_, n_, alloc, data);
	valid_ = true;
}

// Special constructor used when bumping curves only
template <typename Derived>
InterestRateCurveImpl<Derived>::InterestRateCurveImpl(CurveId id, Date as_of_date, Date maturities[], double values[],
						      size_t n, Allocator *alloc, InterpolatorType itype, bool special,
						      DayCountFraction fraction, IRRateType rateType,
						      int deriv_order) noexcept
    : YieldCurve(id), as_of_date_(as_of_date), alloc_(alloc), valueType_(rateType), deriv_order_(deriv_order)
{
	// FIXME will not work for bus_252
	day_fraction_ = get_day_fraction(fraction);
	assert(day_fraction_ != nullptr);
	// This is a private constructor called with all points
	// So no need for special 0 value
	maturities_ = new (*alloc) Date[n];
	std::copy(maturities, maturities + n, maturities_);
	times_ = new (*alloc) double[n];
	times_[0] = 0;
	for (unsigned int i = 1; i < n; i++) {
		times_[i] = day_fraction_->year_fraction(maturities_[0], maturities_[i]);
	}
	n_ = n;
	values_ = new (*alloc) double[n];
	std::copy(values, &values[n], values_);
	valueType_ = rateType;
	InterpolationOptions data;
	data.differentiation_order = deriv_order;
	if (itype == InterpolatorType::MONOTONE_CONVEX) {
		data.monotoneconvex_inputs_are_forwards = (rateType == IRRateType::FORWARD_RATE);
	}
	interpolator_ = make_interpolator(itype, times_, values_, n_, alloc, data);
	valid_ = true;
}

template <typename Derived> void InterestRateCurveImpl<Derived>::dump(FILE *fp) const noexcept
{
	for (size_t i = 1; i < n_; i++) {
		fprintf(fp, "%.15f\t%.15f\n", times_[i], values_[i]);
	}
}

template <typename Derived>
std::vector<std::unique_ptr<YieldCurve, Deleter<YieldCurve>>>
InterestRateCurveImpl<Derived>::get_bumped_curves(Allocator *A, double bp) const noexcept
{
	std::vector<std::unique_ptr<YieldCurve, Deleter<YieldCurve>>> bumped_curves;
	const double h = bp; // 10000.0;
	for (size_t i = 1; i < n_; i++) {
		double *tmp_y = reinterpret_cast<double *>(alloca(sizeof(double) * n_));
		std::copy(values_, values_ + n_, tmp_y);
		tmp_y[i] = tmp_y[i] + h;
		// The constructor will copy the tmp_y values so we are okay
		bumped_curves.push_back(std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
		    new (*A) InterestRateCurveImpl(id(), as_of_date_, maturities_, tmp_y, n_, A, interpolator_->type(),
						   true, day_fraction_->id(), valueType_, 0),
		    Deleter<YieldCurve>(A)));
	}
	return bumped_curves;
}

template <typename Derived>
std::unique_ptr<YieldCurve, Deleter<YieldCurve>> InterestRateCurveImpl<Derived>::get_bumped_curve(Allocator *A, int i,
												  double bp) const
    noexcept
{
	const double h = bp; // 10000.0;
	assert(i >= 1 && i < n_);
	double *tmp_y = reinterpret_cast<double *>(alloca(sizeof(double) * n_));
	std::copy(values_, values_ + n_, tmp_y);
	tmp_y[i] = tmp_y[i] + h;
	return std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
	    new (*A) InterestRateCurveImpl(id(), as_of_date_, maturities_, tmp_y, n_, A, interpolator_->type(), true,
					   day_fraction_->id(), valueType_, 0),
	    Deleter<YieldCurve>(A));
}

template <typename Derived>
IRRateType InterestRateCurveImpl<Derived>::apply_values(InterpolatorType itype, const double *values, size_t n) noexcept
{
	assert(n == n_ - 1);
	if (this->valueType_ == IRRateType::DISCOUNT_FACTOR) {
		if (itype == InterpolatorType::MONOTONE_CONVEX) {
			// Monotone Convex does not support interpolation
			// over discount factors so we convert to ZeroRate
			for (int i = 1; i < n_; i++) {
				double t = times_[i];
				values_[i] = -std::log(values[i - 1]) / t;
			}
			values_[0] = values_[1];
			return IRRateType::ZERO_RATE;
		} else {
			// discount factors - first one must be 1.0
			values_[0] = 1.0;
			std::copy(values, &values[n], values_ + 1);
		}
	} else {
		// The zeroeth value is the same as the first
		// pillar
		values_[0] = values[0];
		std::copy(values, &values[n], values_ + 1);
	}
	return valueType_;
}

template <typename Derived> void InterestRateCurveImpl<Derived>::update_rates(const double *values, size_t n) noexcept
{
	apply_values(interpolator_->type(), values, n);
	interpolator_->update();
}

template <typename Derived> InterestRateCurveImpl<Derived>::~InterestRateCurveImpl() noexcept
{
	alloc_->deallocate(times_);
	alloc_->deallocate(values_);
	alloc_->deallocate(maturities_);
	// if (redukti_debug)
	//  std::cerr << "curve " << name() << " destroyed\n";
}

// Compute the discount function
template <typename Derived> double InterestRateCurveImpl<Derived>::discount(double t) const noexcept
{
	const Derived *implementation = static_cast<const Derived *>(this);
	return implementation->discount_impl(t);
}

template <typename Derived> double InterestRateCurveImpl<Derived>::discount(Date d) const noexcept
{
	double t = time_from_reference(d);
	return discount(t);
}

// Compute the zero rate
// Only continuous compounding supported right now
template <typename Derived> double InterestRateCurveImpl<Derived>::zero_rate(double t) const noexcept
{
	const Derived *implementation = static_cast<const Derived *>(this);
	return implementation->zero_rate_impl(t);
}

template <typename Derived> double InterestRateCurveImpl<Derived>::zero_rate(Date d) const noexcept
{
	double t = time_from_reference(d);
	return zero_rate(t);
}

// Compute instantaneous forward
template <typename Derived> double InterestRateCurveImpl<Derived>::forward(double t) const noexcept
{
	const Derived *implementation = static_cast<const Derived *>(this);
	return implementation->forward_impl(t);
}

template <typename Derived> double InterestRateCurveImpl<Derived>::forward_rate(double t1, double t2) const noexcept
{
	const Derived *implementation = static_cast<const Derived *>(this);
	return implementation->forward_rate_impl(t1, t2);
}

template <typename Derived> double InterestRateCurveImpl<Derived>::forward_rate(Date d1, Date d2) const noexcept
{
	double t = time_from_reference(d1);
	if (d2 == d1) {
		return forward_rate(t, t);
	} else {
		return forward_rate(t, time_from_reference(d2));
	}
}

template <typename Derived>
std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>>
InterestRateCurveImpl<Derived>::get_sensitivities(double x, FixedRegionAllocator *A) const noexcept
{
	return interpolator_->interpolate_with_sensitivities(x, A);
}

template <typename Derived> Date InterestRateCurveImpl<Derived>::last_maturity() const noexcept
{
	return maturities_[n_ - 1];
}

// Implementation where interpolating on zero rates
// for interpolators other than monotone convex
// Continuous compounded rates assumed
class ZeroRateImpl : public InterestRateCurveImpl<ZeroRateImpl>
{
	public:
	double discount_impl(double t) const noexcept { return std::exp(-zero_rate_impl(t) * t); }

	double zero_rate_impl(double t) const noexcept { return interpolator_->interpolate(t); }

	double forward_impl(double t) const noexcept { return forward_rate_impl(t, t); }

	double forward_rate_impl(double t1, double t2) const noexcept
	{
		double delta = 0.0001;
		if (t1 == t2) {
			// instantaneous forward rate is
			// computed based upon a short interval of time
			t1 = std::max(t1 - delta / 2.0, 0.0);
			t2 = t1 + delta;
		}
		double r1 = zero_rate_impl(t1);
		double r2 = zero_rate_impl(t2);
		return (r2 * t2 - r1 * t1) / (t2 - t1);
	}
};

// Monotone convex on zero rates
// Special handling of forwards as the interpolator
// natively supports these
// Continuous compounded rates assumed
class MonotoneConvexImpl : public InterestRateCurveImpl<MonotoneConvexImpl>
{
	public:
	double discount_impl(double t) const noexcept { return std::exp(-zero_rate_impl(t) * t); }

	double zero_rate_impl(double t) const noexcept { return interpolator_->interpolate(t); }

	double forward_impl(double t) const noexcept { return interpolator_->forward(t); }

	double forward_rate_impl(double t1, double t2) const noexcept
	{
		if (t1 == t2) {
			return interpolator_->forward(t1);
		}
		double r1 = zero_rate_impl(t1);
		double r2 = zero_rate_impl(t2);
		return (r2 * t2 - r1 * t1) / (t2 - t1);
	}
};

// Implementation where interpolating on forward rates
// And interpolator other than monotone convex
// Continuous compounded rates assumed
class ForwardRateImpl : public InterestRateCurveImpl<ForwardRateImpl>
{
	public:
	double discount_impl(double t) const noexcept { return std::exp(-zero_rate_impl(t) * t); }

	double zero_rate_impl(double t) const noexcept
	{
		if (t == 0.0)
			return forward_impl(0.0);
		// we need the integral here - underlying
		double integral = interpolator_->primitive(t);
		return integral / t;
	}

	double forward_impl(double t) const noexcept { return interpolator_->interpolate(t); }

	double forward_rate_impl(double t1, double t2) const noexcept
	{
		double delta = 0.0001;
		if (t1 == t2) {
			// instantaneous forward rate is
			// computed based upon a short interval of time
			t1 = std::max(t1 - delta / 2.0, 0.0);
			t2 = t1 + delta;
		}
		double r1 = zero_rate_impl(t1);
		double r2 = zero_rate_impl(t2);
		return (r2 * t2 - r1 * t1) / (t2 - t1);
	}
};

// Implementation where interpolating on discount factors
// And interpolator other than monotone convex
// Continuous compounded rates assumed
class DiscountFactorImpl : public InterestRateCurveImpl<DiscountFactorImpl>
{
	public:
	double discount_impl(double t) const noexcept { return interpolator_->interpolate(t); }

	double zero_rate_impl(double t) const noexcept { return -std::log(discount_impl(t)) / t; }

	double forward_impl(double t) const noexcept { return forward_rate_impl(t, t); }

	double forward_rate_impl(double t1, double t2) const noexcept
	{
		double delta = 0.0001;
		if (t1 == t2) {
			// instantaneous forward rate is
			// computed based upon a short interval of time
			t1 = std::max(t1 - delta / 2.0, 0.0);
			t2 = t1 + delta;
		}
		double r1 = zero_rate_impl(t1);
		double r2 = zero_rate_impl(t2);
		return (r2 * t2 - r1 * t1) / (t2 - t1);
	}
};

/**
 * Construct a curve
 */
std::unique_ptr<YieldCurve, Deleter<YieldCurve>> make_curve(Allocator *A, CurveId id, Date as_of_date,
							    Date maturities[], double values[], size_t n,
							    InterpolatorType interpolator, IRRateType type,
							    int deriv_order, DayCountFraction fraction) noexcept
{
	switch (interpolator) {
	case InterpolatorType::MONOTONE_CONVEX:
		return std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
		    new (*A) InterestRateCurveImpl<MonotoneConvexImpl>(id, as_of_date, maturities, values, n, A,
								       interpolator, fraction, type, deriv_order),
		    Deleter<YieldCurve>(A));
	default:
		switch (type) {
		case IRRateType::DISCOUNT_FACTOR:
			return std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
			    new (*A) InterestRateCurveImpl<DiscountFactorImpl>(
				id, as_of_date, maturities, values, n, A, interpolator, fraction, type, deriv_order),
			    Deleter<YieldCurve>(A));
		case IRRateType::FORWARD_RATE:
			return std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
			    new (*A) InterestRateCurveImpl<ForwardRateImpl>(id, as_of_date, maturities, values, n, A,
									    interpolator, fraction, type, deriv_order),
			    Deleter<YieldCurve>(A));
		case IRRateType::ZERO_RATE:
		default:
			return std::unique_ptr<YieldCurve, Deleter<YieldCurve>>(
			    new (*A) InterestRateCurveImpl<ZeroRateImpl>(id, as_of_date, maturities, values, n, A,
									 interpolator, fraction, type, deriv_order),
			    Deleter<YieldCurve>(A));
		}
	}
}

} // namespace redukti