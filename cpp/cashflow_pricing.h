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
#ifndef _REDUKTI_CASHFLOW_PRICING_H
#define _REDUKTI_CASHFLOW_PRICING_H

#include <cashflow.h>
#include <enums.pb.h>

namespace redukti
{

class Sensitivities;

// First order sensitivities (i.e. delta)
struct Sensitivities1D {
	friend class Sensitivities;

      private:
	Allocator *A_;
	const CurveReference *curve_;
	Sensitivities1D *next_;
	double *delta_;

      private:
	Sensitivities1D(const Sensitivities1D &) = delete;
	Sensitivities1D &operator=(const Sensitivities1D &) = delete;

      public:
	Sensitivities1D(const CurveReference *curve, Allocator *A);
	~Sensitivities1D();
	YieldCurve *curve() const { return curve_->get(); }
	double at(size_t i) const
	{
		assert(i > 0);
		return delta_[i];
	}
	double &at(size_t i)
	{
		assert(i > 0);
		return delta_[i];
	}
	int count() const { return curve_->get()->last_pillar() + 1; }
	int count_differences(Sensitivities1D *other, FILE *f = stdout, const char *desc = "") const;
	Sensitivities1D *next() const { return next_; }
};

// Second order sensitivities (i.e. gamma)
struct Sensitivities2D {
	friend class Sensitivities;

      private:
	Allocator *A_;
	const CurveReference *curve1_;
	const CurveReference *curve2_;
	int m_, n_;
	Sensitivities2D *next_;
	double *gamma_;

      private:
	Sensitivities2D(const Sensitivities2D &) = delete;
	Sensitivities2D &operator=(const Sensitivities2D &) = delete;

      public:
	Sensitivities2D(const CurveReference *curve1, const CurveReference *curve2, Allocator *A);
	~Sensitivities2D();
	YieldCurve *curve1() const { return curve1_->get(); }
	YieldCurve *curve2() const { return curve2_->get(); }
	double at(size_t i, size_t j) const
	{
		assert(i > 0);
		assert(j > 0);
		return gamma_[i * n_ + j];
	}
	double &at(size_t i, size_t j)
	{
		assert(i > 0);
		assert(j > 0);
		return gamma_[i * n_ + j];
	}
	int count1() const { return m_; }
	int count2() const { return n_; }
	int count_differences(Sensitivities2D *other, FILE *f = stdout, const char *desc = "") const;
	Sensitivities2D *next() const { return next_; }
};

class Sensitivities
{
      private:
	Allocator *A_;
	Sensitivities1D *first_1d_;
	Sensitivities1D *last_1d_;
	Sensitivities2D *first_2d_;
	Sensitivities2D *last_2d_;
	std::array<CurveWrapper, 20> curves_;

      private:
	void add_first_order_sensitivities(Sensitivities1D *d)
	{
		if (!first_1d_)
			first_1d_ = d;
		if (last_1d_)
			last_1d_->next_ = d;
		last_1d_ = d;
	}
	void add_second_order_sensitivities(Sensitivities2D *g)
	{
		if (!first_2d_)
			first_2d_ = g;
		if (last_2d_)
			last_2d_->next_ = g;
		last_2d_ = g;
	}

      public:
	Sensitivities(Allocator *A);
	~Sensitivities();
	// Find or add
	Sensitivities1D *first_order_sensitivities(YieldCurve *curve);
	// Find or add
	Sensitivities2D *second_order_sensitivities(YieldCurve *curve1, YieldCurve *curve2);
	Sensitivities1D *find_first_order_sensitivities(YieldCurve *curve) const
	{
		return find_first_order_sensitivities(curve->id());
	}
	Sensitivities1D *find_first_order_sensitivities(CurveId id) const;
	Sensitivities2D *find_second_order_sensitivities(YieldCurve *curve1, YieldCurve *curve2) const
	{
		return find_second_order_sensitivities(curve1->id(), curve2->id());
	}
	Sensitivities2D *find_second_order_sensitivities(CurveId id1, CurveId id2) const;

	void dump(FILE *fp) const;
	int count_differences(Sensitivities *other, FILE *f = stdout, const char *prefix = "") const;
	void reset();
	// Find or add
	const CurveReference *get(YieldCurve *curve);
	const CurveReference *find(YieldCurve *curve) const;
	Sensitivities1D *first_1d_sensitivities() const { return first_1d_; }
	Sensitivities2D *first_2d_sensitivities() const { return first_2d_; }
	void get_curve_ids(std::vector<CurveId> &ids) const;
};

class Cashflows;

// Calculate sensitivities (delta and gamma) numerically
// and store in supplied container
extern void compute_sensitivity_numerically(FixedRegionAllocator *allocator, const Cashflows *flows,
					    const CurveReference *discount_curve, const CurveReference *forward_curve1,
					    const CurveReference *forward_curve2, Sensitivities *sensitivities,
					    StatusCode &status, double h = 0.00001);

// When cashflows are defined, they reference logical curves via
// PricingCurve identifiers. At the time of valuation these logical curves
// must be mapped to physical instances of curves - the CurveProvider
// interfaces= defines such a component.
class CurveProvider
{
      public:
	virtual ~CurveProvider() {}
	virtual const CurveReference *get_curve(PricingCurve curve) const = 0;
};

// Calculate PV and if ValuationContext.derivative_order > 0
// then also delta and gamma
extern double compute_present_value(FixedRegionAllocator *A, const ValuationContext &ctx, const Cashflows *flows,
				    const CurveProvider *mapping_provider, Sensitivities &sensitivities,
				    StatusCode &status);

extern double compute_present_value(FixedRegionAllocator *A, const Cashflows *flows,
				    const CurveReference *discount_curve, const CurveReference *forward_curve1,
				    const CurveReference *forward_curve2, const ValuationContext &ctx,
				    Sensitivities &sensitivities, StatusCode &status);

int test_pricing();

} // namespace redukti

#endif