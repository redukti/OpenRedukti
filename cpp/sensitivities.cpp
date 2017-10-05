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

#include <allocators.h>
#include <cashflow.h>
#include <index.h>
#include <logger.h>
#include <schedule.h>

#include <cashflow_pricing.h>
#include <internal/cashflow_internal.h>
#include <internal/cashflow_pricing_internal.h>

#include <map>
#include <sstream>
#include <vector>

#include <cctype>
#include <cmath>
#include <cstdio>

namespace redukti
{

inline bool is_zero(double v) { return v == 0.0; }

Sensitivities1D::Sensitivities1D(const CurveReference *curve, Allocator *A) : A_(A), curve_(curve), next_(nullptr)
{
	auto size = (curve->get()->last_pillar() + 1) * sizeof(double);
	delta_ = (double *)A->allocate(size);
	memset(delta_, 0, size);
}

Sensitivities1D::~Sensitivities1D() { A_->deallocate(delta_); }

int Sensitivities1D::count_differences(Sensitivities1D *o, FILE *f, const char *desc) const
{
	int errors = 0;
	if (o->curve_ == curve_ || o->curve_->get()->id() == curve_->get()->id()) {
		const CurveReference *yc = curve_;
		for (int i = 1; i <= yc->get()->last_pillar(); i++) {
			double v1 = at(i);
			double v2 = o->at(i);
			double diff = v1 - v2;
			if (std::fabs(v1) > 10.0 || std::fabs(v2) > 10.0) {
				if (std::fabs(diff) > 1.0 && std::fabs(diff / v1) > 0.01) {
					// fprintf(stderr, "delta mismatch
					// %s\t[%d]\t%f\t%f\n",
					// yc->get()->name(), i, at(i),
					// o->at(i));
					errors++;
				}
			}
			// fprintf(f,
			// "%s\tdelta\t%s\t\t%d\t\t%.12f\t%.12f\t%.12f\n", desc,
			//        yc->get()->name().c_str(), i, v1, v2, diff);
		}
	} else {
		errors++;
	}
	return errors;
}

Sensitivities2D::Sensitivities2D(const CurveReference *curve1, const CurveReference *curve2, Allocator *A)
    : A_(A), curve1_(curve1), curve2_(curve2), m_(curve1->get()->last_pillar() + 1),
      n_(curve2->get()->last_pillar() + 1), next_(nullptr)
{
	auto size = m_ * n_ * sizeof(double);
	gamma_ = (double *)A->allocate(size);
	memset(gamma_, 0, size);
}

Sensitivities2D::~Sensitivities2D() { A_->deallocate(gamma_); }

int Sensitivities2D::count_differences(Sensitivities2D *o, FILE *f, const char *desc) const
{
	int errors = 0;
	if (o->curve1_->get()->id() == curve1_->get()->id() && o->curve2_->get()->id() == curve2_->get()->id()) {
		YieldCurve *yc1 = curve1_->get();
		YieldCurve *yc2 = curve2_->get();
		for (int i = 1; i <= yc1->last_pillar(); i++) {
			for (int j = 1; j <= yc2->last_pillar(); j++) {
				double v1 = at(i, j);
				double v2 = o->at(i, j);
				double diff = v1 - v2;
				if (std::fabs(v1) > 50.0 || std::fabs(v2) > 50.0) {
					if (std::fabs(diff) > 500.0 && std::fabs(diff / v1) > 0.1) {
						// fprintf(stderr, "gamma
						// mismatch\t[%d][%d]\t%f\t%f\n",
						// i, j, at(i, j), o->at(i, j));
						errors++;
					}
				}
				// fprintf(f,
				// "%s\tgamma\t%s\t%s\t%d\t%d\t%.12f\t%.12f\t%.12f\n",
				//        desc, yc1->name().c_str(),
				//        yc2->name().c_str(), i, j, v1, v2,
				//        diff);
			}
		}
	} else if (o->curve1_->get()->id() == curve2_->get()->id() && o->curve2_->get()->id() == curve1_->get()->id()) {
		YieldCurve *yc1 = curve1_->get();
		YieldCurve *yc2 = curve2_->get();
		for (int i = 1; i <= yc1->last_pillar(); i++) {
			for (int j = 1; j <= yc2->last_pillar(); j++) {
				double v1 = at(i, j);
				double v2 = o->at(j, i);
				double diff = v1 - v2;
				if (std::fabs(v1) > 50.0 || std::fabs(v2) > 50.0) {
					if (std::fabs(diff) > 500.0 && std::fabs(diff / v1) > 0.1) {
						// fprintf(stderr, "gamma
						// mismatch\t[%d][%d]\t%f\t%f\n",
						// i, j, at(i, j), o->at(j, i));
						errors++;
					}
				}
				// fprintf(f,
				// "%s\tgamma\t%s\t%s\t%d\t%d\t%.12f\t%.12f\t%.12f\n",
				//        desc, yc1->name().c_str(),
				//        yc2->name().c_str(), i, j, v1, v2,
				//        diff);
			}
		}
	} else {
		errors++;
	}
	return errors;
}

Sensitivities::Sensitivities(Allocator *A)
    : A_(A), first_1d_(nullptr), last_1d_(nullptr), first_2d_(nullptr), last_2d_(nullptr)
{
}

Sensitivities::~Sensitivities() { reset(); }

void Sensitivities::reset()
{
	while (first_1d_) {
		auto next = first_1d_->next();
		first_1d_->~Sensitivities1D();
		A_->deallocate(reinterpret_cast<void *>(first_1d_));
		first_1d_ = next;
	}
	while (first_2d_) {
		auto next = first_2d_->next();
		first_2d_->~Sensitivities2D();
		A_->deallocate(reinterpret_cast<void *>(first_2d_));
		first_2d_ = next;
	}
	first_1d_ = last_1d_ = nullptr;
	first_2d_ = last_2d_ = nullptr;
	for (int i = 0; i < curves_.size(); i++)
		curves_[i].set(nullptr);
}

const CurveReference *Sensitivities::find(YieldCurve *curve) const
{
	for (size_t i = 0; i < curves_.size(); i++) {
		auto cr = curves_[i].get();
		if (!cr)
			break;
		if (curve->id() == cr->id()) {
			return &curves_[i];
		}
	}
	return nullptr;
}

const CurveReference *Sensitivities::get(YieldCurve *curve)
{
	for (size_t i = 0; i < curves_.size(); i++) {
		if (curves_[i].get() == nullptr) {
			curves_[i].set(curve);
			return &curves_[i];
		}
		if (curve->id() == curves_[i].get()->id()) {
			return &curves_[i];
		}
	}
	die("Attempt to add more curves than currently allowed limit of %d\n", (int)curves_.size());
	return nullptr;
}

Sensitivities1D *Sensitivities::first_order_sensitivities(YieldCurve *curve)
{
	auto cr = get(curve);
	for (Sensitivities1D *d = first_1d_sensitivities(); d != nullptr; d = d->next()) {
		if (d->curve()->id() == cr->get()->id()) {
			return d;
		}
	}
	Sensitivities1D *d = new (*A_) Sensitivities1D(cr, A_);
	add_first_order_sensitivities(d);
	return d;
}

Sensitivities2D *Sensitivities::second_order_sensitivities(YieldCurve *c1, YieldCurve *c2)
{
	auto curve1 = get(c1);
	auto curve2 = get(c2);
	for (Sensitivities2D *g = first_2d_sensitivities(); g != nullptr; g = g->next()) {
		if (g->curve1()->id() == curve1->get()->id() && g->curve2()->id() == curve2->get()->id()) {
			return g;
		}
	}
	Sensitivities2D *g = new (*A_) Sensitivities2D(curve1, curve2, A_);
	add_second_order_sensitivities(g);
	return g;
}

Sensitivities1D *Sensitivities::find_first_order_sensitivities(CurveId id) const
{
	for (Sensitivities1D *d = first_1d_sensitivities(); d != nullptr; d = d->next()) {
		if (d->curve()->id() == id) {
			return d;
		}
	}
	return nullptr;
}

Sensitivities2D *Sensitivities::find_second_order_sensitivities(CurveId id1, CurveId id2) const
{
	for (Sensitivities2D *g = first_2d_sensitivities(); g != nullptr; g = g->next()) {
		if ((g->curve1()->id() == id1 && g->curve2()->id() == id2) ||
		    (g->curve1()->id() == id2 && g->curve2()->id() == id1)) {
			return g;
		}
	}
	return nullptr;
}

void Sensitivities::get_curve_ids(std::vector<CurveId> &ids) const
{
	for (int i = 0; i < curves_.size(); i++) {
		auto curve = curves_[i].get();
		if (curve)
			ids.push_back(curve->id());
		else
			break;
	}
}

void Sensitivities::dump(FILE *fp) const
{
	for (Sensitivities1D *ds = first_1d_sensitivities(); ds != nullptr; ds = ds->next()) {
		YieldCurve *yc = ds->curve();
		fprintf(fp, "Delta Sensitivies: %s\n", yc->name().c_str());
		for (int i = 1; i <= yc->last_pillar(); i++) {
			if (!is_zero(ds->at(i)))
				fprintf(fp, "delta\t[%d]\t%f\n", i, ds->at(i));
		}
	}
	for (Sensitivities2D *gs = first_2d_sensitivities(); gs != nullptr; gs = gs->next()) {
		YieldCurve *yc1 = gs->curve1();
		YieldCurve *yc2 = gs->curve2();
		fprintf(fp, "Gamma Sensitivies: %s %s\n", yc1->name().c_str(), yc2->name().c_str());
		for (int i = 1; i <= yc1->last_pillar(); i++) {
			for (int j = 1; j <= yc2->last_pillar(); j++) {
				if (!is_zero(gs->at(i, j)))
					fprintf(fp, "gamma\t[%d]\t[%d]\t%f\n", i, j, gs->at(i, j));
			}
		}
	}
}

int Sensitivities::count_differences(Sensitivities *other, FILE *f, const char *prefix) const
{
	int errors = 0;
	for (Sensitivities1D *ds = first_1d_sensitivities(); ds != nullptr; ds = ds->next()) {
		auto ds1 = other->first_order_sensitivities(ds->curve());
		errors += ds->count_differences(ds1, f, prefix);
	}
	for (Sensitivities2D *gs = first_2d_sensitivities(); gs != nullptr; gs = gs->next()) {
		auto gs1 = other->find_second_order_sensitivities(gs->curve1()->id(), gs->curve2()->id());
		if (gs1 == nullptr)
			errors++;
		else {
			errors += gs->count_differences(gs1, f, prefix);
		}
	}
	return errors;
}

// aggregate delta sensitivities to curve pillar points
// input is sensitivities against cashflow dates (data)
// resolver provides the variables that are specific to the curve
static void aggregate_delta(FixedRegionAllocator *A, const VariableResolver *resolver, const CurveReference *curve,
			    Sensitivities1D *delta, redukti_adouble_t &data)
{
	// Loop through the variables that are associated with curve
	for (auto ptr = resolver->begin(); ptr != resolver->end(); ++ptr) {
		// cashflow date
		auto date = ptr.date();
		// variable number
		auto var = ptr.var();
		// get the first order derivative of PV to cashflow date on the
		// curve
		double d = redukti_adouble_get_derivative1(&data, var);
		if (!is_zero(d)) {
			FixedRegionAllocatorGuard g(A);
			// project to curve pillars
			// first get sensitivities of cashflow date to pillars
			auto node_sensitivities =
			    curve->get()->get_sensitivities(curve->get()->time_from_reference(date), A);
			// loop through the curve variables i.e. pillars
			for (int x = 1; x < node_sensitivities->vars_; x++) {
				double node_sens = redukti_adouble_get_derivative1(node_sensitivities.get(), x);
				if (!is_zero(node_sens))
					delta->at(x) += node_sens * d;
			}
		}
	}
}

// Vector (first order) cashflow sens * Second order curve node sens
static void aggregate_gamma2(FixedRegionAllocator *A, const VariableResolver *resolver, const CurveReference *curve,
			     Sensitivities2D *gamma, redukti_adouble_t &data)
{
	if (curve->get()->interpolator_type() == InterpolatorType::LINEAR)
		return;
	// loop through cashflow dates
	for (auto ptr = resolver->begin(); ptr != resolver->end(); ++ptr) {
		// cashflow date
		auto date = ptr.date();
		// variable number
		auto var = ptr.var();
		// Get first derivative with respect to cashflow date (from
		// vector)
		double first_order_sens = redukti_adouble_get_derivative1(&data, var);
		if (!is_zero(first_order_sens)) {
			FixedRegionAllocatorGuard g(A);
			// Get node sensitivities
			auto node_sens = curve->get()->get_sensitivities(curve->get()->time_from_reference(date), A);
			// for each i, j pillars
			for (int i = 1; i < node_sens->vars_; i++) {
				for (int j = 1; j < node_sens->vars_; j++) {
					// get second order sens
					double second_order_node_sens =
					    redukti_adouble_get_derivative2(node_sens.get(), i, j);
					if (!is_zero(second_order_node_sens)) {
						// aggregate
						gamma->at(i, j) += second_order_node_sens * first_order_sens;
					}
				}
			}
		}
	}
}

// second order cashflow sens * first order curve sens
static void aggregate_gamma(FixedRegionAllocator *A, const VariableResolver *resolver1,
			    const VariableResolver *resolver2, const CurveReference *curve1,
			    const CurveReference *curve2, Sensitivities2D *gamma, redukti_adouble_t &data)
{
	// Loop through the variables / cashflow dates that are associated with
	// curve1
	for (auto ptr1 = resolver1->begin(); ptr1 != resolver1->end(); ++ptr1) {
		// cashflow date
		auto date1 = ptr1.date();
		// variable number
		auto var1 = ptr1.var();
		// Loop through the variables / dates on the second curve.
		for (auto ptr2 = resolver2->begin(); ptr2 != resolver2->end(); ++ptr2) {
			// second cashflow date
			auto date2 = ptr2.date();
			// variable number
			auto var2 = ptr2.var();
			// get second order cashflow sens
			double second_order_cashflow_sens = redukti_adouble_get_derivative2(&data, var1, var2);
			if (!is_zero(second_order_cashflow_sens)) {
				FixedRegionAllocatorGuard g(A);
				auto node_sens1 =
				    curve1->get()->get_sensitivities(curve1->get()->time_from_reference(date1), A);
				auto node_sens2 =
				    curve2->get()->get_sensitivities(curve2->get()->time_from_reference(date2), A);
				// loop through each of the curve pillars
				for (int n1 = 1; n1 < node_sens1->vars_; n1++) {
					for (int n2 = 1; n2 < node_sens2->vars_; n2++) {
						double s1 = redukti_adouble_get_derivative1(node_sens1.get(), n1);
						double s2 = redukti_adouble_get_derivative1(node_sens2.get(), n2);
						if (!is_zero(s1) && !is_zero(s2)) {
							gamma->at(n1, n2) += second_order_cashflow_sens * s1 * s2;
						}
					}
				}
			}
		}
	}

	// If the curves are the same then we need to aggregate
	// first order cashflow sens * second order curve node sens
	if (curve1->get()->id() == curve2->get()->id()) {
		aggregate_gamma2(A, resolver1, curve1, gamma, data);
		aggregate_gamma2(A, resolver2, curve2, gamma, data);
	}
}

PricingCurve make_pricing_curve(PricingCurveType type, CurveId id)
{
	Currency ccy;
	IndexFamily index_family;
	Tenor tenor;
	Date as_of_date;
	MarketDataQualifier qual;
	short int cycle;
	short int scenario;
	PricingCurveType curve_type;
	bool status = curve_id_components(id, curve_type, ccy, index_family, tenor, as_of_date, cycle, qual, scenario);
	assert(status);
	return PricingCurve(type, ccy, index_family, tenor);
}

// aggregate and project sensitivities to curve pillars
void aggregate_sensitivity(FixedRegionAllocator *A, const VariableResolver *discount_resolver, const CurveReference *dc,
			   const VariableResolver *forward1_resolver, const CurveReference *f1,
			   const VariableResolver *forward2_resolver, const CurveReference *f2,
			   Sensitivities *container, redukti_adouble_t &data, StatusCode &status)
{
	assert(discount_resolver);
	assert(forward1_resolver);
	assert(forward2_resolver);

	const CurveReference *discount_curve = container->get(dc->get());
	const CurveReference *forward_curve1 = container->get(f1->get());
	const CurveReference *forward_curve2 = container->get(f2->get());

	if (discount_curve == nullptr || forward_curve1 == nullptr || forward_curve2 == nullptr) {
		fprintf(stderr, "aggregate_sensitivity: Error: one or more input "
				"curves are null\n");
		status = StatusCode::kInternalError;
		return;
	}

	auto discount_delta = container->first_order_sensitivities(discount_curve->get());
	auto discount_discount_gamma =
	    container->second_order_sensitivities(discount_curve->get(), discount_curve->get());
	auto forward1_delta = container->first_order_sensitivities(forward_curve1->get());
	auto forward1_forward1_gamma =
	    container->second_order_sensitivities(forward_curve1->get(), forward_curve1->get());
	auto forward1_discount_gamma =
	    container->second_order_sensitivities(forward_curve1->get(), discount_curve->get());

	aggregate_delta(A, discount_resolver, discount_curve, discount_delta, data);
	if (discount_curve->get()->id() != forward_curve1->get()->id()) {
		aggregate_delta(A, forward1_resolver, forward_curve1, forward1_delta, data);
	}
	aggregate_gamma(A, discount_resolver, discount_resolver, discount_curve, discount_curve,
			discount_discount_gamma, data);
	if (discount_curve->get()->id() != forward_curve1->get()->id()) {
		aggregate_gamma(A, forward1_resolver, forward1_resolver, forward_curve1, forward_curve1,
				forward1_forward1_gamma, data);
		aggregate_gamma(A, forward1_resolver, discount_resolver, forward_curve1, discount_curve,
				forward1_discount_gamma, data);
	}
	if (forward_curve1->get()->id() != forward_curve2->get()->id()) {
		// second forward curve is different from first forward curve
		auto forward2_delta = container->first_order_sensitivities(forward_curve2->get());
		auto forward2_forward2_gamma =
		    container->second_order_sensitivities(forward_curve2->get(), forward_curve2->get());
		auto forward2_discount_gamma =
		    container->second_order_sensitivities(forward_curve2->get(), discount_curve->get());
		auto forward1_forward2_gamma =
		    container->second_order_sensitivities(forward_curve1->get(), forward_curve2->get());

		aggregate_delta(A, forward2_resolver, forward_curve2, forward2_delta, data);
		aggregate_gamma(A, forward2_resolver, forward2_resolver, forward_curve2, forward_curve2,
				forward2_forward2_gamma, data);
		aggregate_gamma(A, forward2_resolver, discount_resolver, forward_curve2, discount_curve,
				forward2_discount_gamma, data);
		aggregate_gamma(A, forward1_resolver, forward2_resolver, forward_curve1, forward_curve2,
				forward1_forward2_gamma, data);
	}
}

void dump_sensitivity(redukti_adouble_t &v, VariableResolvers &vars_)
{
	for (size_t i = 0; i < vars_.num_resolvers(); i++) {
		auto resolver = vars_.get_resolver_at(i);
		// loop through cashflow dates
		for (auto ptr = resolver.begin(); ptr != resolver.end(); ++ptr) {
			// cashflow date
			auto date = ptr.date();
			// variable number
			auto var = ptr.var();
			// Get first derivative with respect to cashflow date
			// (from vector)
			double first_order_sens = redukti_adouble_get_derivative1(&v, var);
			if (!is_zero(first_order_sens)) {
				auto ymd = date_components(date);
				fprintf(stderr, "%s\tdelta[%02d-%02d-%d]\t%f\n", resolver.curve().name().c_str(),
					(int)ymd.d, (int)ymd.m, ymd.y, first_order_sens);
			}
		}
	}
}

} // namespace redukti
