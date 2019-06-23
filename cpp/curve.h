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
#ifndef _REDUKTI_CURVE_H_
#define _REDUKTI_CURVE_H_

#include <enums.pb.h>

#include <allocators.h>
#include <autodiff.h>
#include <date.h>
#include <dayfractions.h>
#include <interpolators.h>

#include <array>
#include <inttypes.h>

namespace redukti
{
class IRCurveDefinition;
class ZeroCurve;

// Curve identifier
typedef uint64_t CurveId;

// Constructs a curve id by combining the constituents
extern CurveId make_curve_id(PricingCurveType type, Currency ccy, IndexFamily index_family, Tenor tenor,
			     Date as_of_date, short int cycle = 0,
			     MarketDataQualifier qual = MarketDataQualifier::MDQ_NORMAL, short int scenario = 0);
// Extracts the constituents from a curve id
extern bool curve_id_components(CurveId id, PricingCurveType &type, Currency &ccy, IndexFamily &index_family,
				Tenor &tenor, Date &as_of_date, short int &cycle, MarketDataQualifier &qual,
				short int &scenario);

// Gets a string representation of the curve Id,
// note that this is an expensive operation so use only for
// debugging
extern std::string curve_id_to_string(CurveId id);

class Curve
{
	public:
	virtual ~Curve() noexcept {}
	double time_from_reference(Date d) const noexcept { return day_fraction().year_fraction(as_of_date(), d); }
	virtual const DayFraction &day_fraction() const noexcept = 0;
	virtual Date as_of_date() const noexcept = 0;
	CurveId id() const noexcept { return curve_id_; }
	std::string name() const noexcept { return curve_id_to_string(curve_id_); }
	virtual bool is_valid() const noexcept { return false; }
	virtual bool curve_sensitivities_supported() const noexcept
	{
		return curve_type() == CurveType::CURVE_TYPE_INTERPOLATED;
	}
	virtual CurveType curve_type() const noexcept = 0;

	private:
	Curve(const Curve &) = delete;
	Curve &operator=(const Curve &) = delete;

	protected:
	Curve(CurveId id) noexcept;
	CurveId curve_id_;
};

typedef std::unique_ptr<redukti_adouble_t, Deleter<redukti_adouble_t>> CurveSensitivitiesPointerType;

class YieldCurve : public Curve
{
	public:
	virtual ~YieldCurve() noexcept {}

	virtual double discount(double time) const noexcept = 0;

	// Discount factors
	// These methods return the discount factor from a given date or time
	// to the reference date.  In the latter case, the time is calculated
	// as a fraction of year from the reference date.
	virtual double discount(Date d) const noexcept = 0;

	// Zero-yield rates
	// These methods return the implied zero-yield rate for a
	// given date or time.  In the former case, the time is
	// calculated as a fraction of year from the reference date.
	virtual double zero_rate(Date d) const noexcept = 0;

	// The resulting interest rate has the same day-counting rule
	// used by the term structure. The same rule should be used
	// for calculating the passed time t.
	virtual double zero_rate(double t) const noexcept = 0;

	// Forward rates
	// These methods returns the forward interest rate between two dates
	// or times.  In the former case, times are calculated as fractions
	// of year from the reference date.
	// If both dates (times) are equal the instantaneous forward rate is
	// returned.
	virtual double forward_rate(Date d1, Date d2) const noexcept = 0;

	// The resulting interest rate has the same day-counting rule
	// used by the term structure. The same rule should be used
	// for calculating the passed times t1 and t2.
	virtual double forward_rate(double t1, double t2) const noexcept = 0;

	// Instantaneous forward rate
	virtual double forward(double t) const noexcept = 0;

	// The offset of the last pillar for interpolated curves
	// The first pillar is numbered 1.
	// For parametric curves offset of last parameter
	virtual int last_pillar() const noexcept = 0;

	// Update the rates or parameters
	virtual void update_rates(const double *rates, size_t n) noexcept = 0;

	// Value at pillar point for interpolated curves
	// Parameter value for parametric curves
	virtual double value(int pillar) const noexcept = 0;

	// maturity time from ref date
	// Only available on interpolated curves
	virtual double maturity_time(int pillar) const = 0;

	// maturity date for a pillar
	// Only available on interpolated curves
	virtual Date maturity_date(int pillar) const = 0;

	// Only available on interpolated curves
	double last_maturity_time() const { return maturity_time(last_pillar()); }

	// Gets the sensitivities to pillars using the underlying
	// interpolator.
	// Note: only available on interolated curves
	virtual CurveSensitivitiesPointerType get_sensitivities(double x, FixedRegionAllocator *A) const noexcept = 0;

	// Only available on interpolated curves
	virtual std::vector<std::unique_ptr<YieldCurve, Deleter<YieldCurve>>>
	get_bumped_curves(Allocator *A, double h = 0.00001) const noexcept = 0;
	virtual std::unique_ptr<YieldCurve, Deleter<YieldCurve>> get_bumped_curve(Allocator *A, int pillar,
										  double h = 0.00001) const
	    noexcept = 0;

	virtual void dump(FILE *fp = stderr) const noexcept = 0;

	virtual InterpolatorType interpolator_type() const noexcept = 0;

	// Only available on interpolated curves
	virtual Date last_maturity() const noexcept = 0;

	protected:
	YieldCurve(CurveId id) noexcept : Curve(id) {}

	private:
	YieldCurve(const YieldCurve &) = delete;
	YieldCurve &operator=(const YieldCurve &) = delete;
};

// When referencing a curve it is useful to have some
// indirection as this allows the curve to be modified without
// affecting the client code. This is particularly needed when
// bootstrapping curves. The CurveReference interface provides this
// indirection.
class CurveReference
{
	public:
	virtual ~CurveReference() noexcept {}
	virtual YieldCurve *get() const noexcept = 0;
};

// Wraps a curve pointer
class CurveWrapper : public CurveReference
{
	private:
	YieldCurve *curve_;

	public:
	CurveWrapper(YieldCurve *curve = nullptr) noexcept : curve_(curve) {}
	virtual YieldCurve *get() const noexcept { return curve_; }
	void set(YieldCurve *c) { curve_ = c; }
};

// Construct a curve
// @param A - Memory allocator
// @param id - ID of the curve
// @param as_of_date - As of date
// @param maturities - Curve pillar points
// @param values - interpretation depends upon type below
// @param n - Size of the arrays above
// @param interpolator - Type of interpolator to be used
// @param rateType - ZeroRate, DiscountFactor or FowardRate
// @param derive_order - the order to which node sensitivities are to be
// computed
// @fraction - day count fraction
//
// Note that the curve object will copy the maturities and values arrays
// so caller need not retain these arrays. Since the arrays are copied
// changes to original values do not impact the curve. You can invoke
// the method update_rates() to update the values after the curve is
// created.

typedef std::unique_ptr<YieldCurve, Deleter<YieldCurve>> YieldCurvePointerType;
extern YieldCurvePointerType make_curve(Allocator *A, CurveId id, Date as_of_date, Date maturities[], double values[],
					size_t n, InterpolatorType interpolator,
					IRRateType type = IRRateType::ZERO_RATE, int deriv_order = 0,
					DayCountFraction fraction = DayCountFraction::ACT_365_FIXED) noexcept;

extern YieldCurvePointerType make_svensson_curve(Allocator *A, CurveId id, Date as_of_date, double parameters[], size_t n,
						 DayCountFraction fraction = DayCountFraction::ACT_365_FIXED) noexcept;

extern YieldCurvePointerType make_curve(Date as_of_date, const IRCurveDefinition *defn, const ZeroCurve &curve,
					int deriv_order, PricingCurveType type = PRICING_CURVE_TYPE_UNSPECIFIED,
					MarketDataQualifier mdq = MDQ_NORMAL, short int cycle = 0,
					short int scenario = 0);

// extern int test_curves();

} // namespace redukti

#endif
