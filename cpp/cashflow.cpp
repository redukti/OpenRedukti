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

#include <cashflow.h>

#include <converters.h>
#include <index.h>
#include <logger.h>
#include <schedule.h>
#include <status.h>

#include <cashflow.pb.h>
#include <internal/cashflow_internal.h>

#include <map>
#include <sstream>
#include <vector>

#include <assert.h>
#include <cctype>
#include <cstdio>

using namespace google::protobuf;

namespace redukti
{

// The conversion of a PricingCurve to string
// is expensive so should only be used for
// debugging
std::string PricingCurve::name() const
{
	auto converter = get_default_converter();
	PricingCurveType type = curve_type();
	Currency ccy = currency();
	IndexFamily family = index_family();
	Tenor t = tenor();
	char buf[100];
	snprintf(buf, sizeof buf, "PricingCurve[%s-%s-%s-%s]",
		 type == PRICING_CURVE_TYPE_DISCOUNT ? "Discount"
						     : (type == PRICING_CURVE_TYPE_FORWARD ? "Forward" : ""),
		 converter->currency_to_string(ccy), converter->index_family_to_string(family),
		 converter->tenor_to_string(t).c_str());
	return std::string(buf);
}

void VariableResolvers::reset()
{
	for (int i = 0; i < MAX_CURVES_PER_RESOLVER; i++)
		curves_[i] = PricingCurve();
	memset(&data_[0], 0, sizeof(data_[0]) * datasize_);
}

// how many variables so we have?
int VariableResolvers::num_vars() const
{
	int i = 0;
	for (int j = 0; j < datasize_; j++)
		if (!data_[j].is_null())
			i++;
		else
			break;
	return i;
}

// Search and add a curve to the resolvers set
// If successful then the returned resolver will
// have give curve's location
// A value of -1 indicates failure
VariableResolver VariableResolvers::add_resolver(PricingCurve curve)
{
	for (int i = 0; i < MAX_CURVES_PER_RESOLVER; i++) {
		if (curves_[i] == curve) {
			return VariableResolver{this, i};
		} else if (curves_[i].id() == 0) {
			curves_[i] = curve;
			return VariableResolver{this, i};
		}
	}
	return VariableResolver{this, -1};
}

// search for a var resolver by curve
// If not found returns resolver with curve = -1
VariableResolver VariableResolvers::get_resolver(PricingCurve curve)
{
	for (int i = 0; i < MAX_CURVES_PER_RESOLVER; i++) {
		if (curves_[i] == curve) {
			return VariableResolver{this, i};
		}
	}
	return VariableResolver{this, -1};
}

PricingCurve VariableResolvers::get_curve(int32_t id) const
{
	if (id < 0 || id >= MAX_CURVES_PER_RESOLVER)
		return PricingCurve();
	return curves_[id];
}

PricingCurve VariableResolver::curve() const { return resolvers_->get_curve(curve_); }

// Obtain a registered curve of specified type
// Optional index specifies which relative curve to get
// E.g. whether to get the first (index = 0) forward curve
// or second (index = 1) - note that first or second do not
// have any meaning other than that may return different
// curves
PricingCurve VariableResolvers::get_curve_of_type(PricingCurveType type, int index) const
{
	int n = -1;
	for (int i = 0; i < MAX_CURVES_PER_RESOLVER; i++) {
		if (curves_[i].curve_type() == type) {
			n++;
			if (n == index)
				return curves_[i];
		}
	}
	return PricingCurve();
}

// How many resolvers, i.e. curves
size_t VariableResolvers::num_resolvers() const
{
	size_t n = 0;
	for (int i = 0; i < MAX_CURVES_PER_RESOLVER; i++)
		if (curves_[i].id() != 0)
			n++;
		else
			break;
	return n;
}

VariableResolver VariableResolvers::get_resolver_at(int32_t curve_id)
{
	if (curve_id < 0 || curve_id >= MAX_CURVES_PER_RESOLVER)
		return VariableResolver(this, -1);
	PricingCurve curve = curves_[curve_id];
	if (curve.id() == 0)
		return VariableResolver(this, -1);
	return VariableResolver(this, curve_id);
}

// Searches for a var by curve id and date
// if found returns 1, and *varp and *curvep are set
// If error (no more space or invalid curve_id) return -1, and *varp/*curvep are
// not set
// If not found returns 0 and *var is set to the next available slot
// This interface  allows this to be reused in add_var()
int VariableResolvers::search_var(int curve_id, Date date, PricingCurve *curvep, int *varp) const
{
	if (curve_id < 0 || curve_id >= MAX_CURVES_PER_RESOLVER)
		return -1;
	PricingCurve curve = curves_[curve_id];
	if (curve.id() == 0)
		return -1;
	// printf("Searching for curveid %d date %d\n", curve_id, date);
	for (int i = 0; i < datasize_; i++) {
		Date thisdate = data_[i].date;
		if (thisdate == 0) {
			*varp = i;
			*curvep = curve;
			return 0;
		} else {
			// printf("[%d] date = %d\n", i, thisdate);
			if (thisdate == date) {
				uint32_t cid = data_[i].curveid;
				// printf("[%d] curveid = %d\n", i, cid);
				if (cid == curve_id) {
					*varp = i;
					*curvep = curve;
					return 1;
				}
			}
		}
	}
	return -1;
}

// Searches for next date/var for given curve
// search starts as pos
// returns new pos which can be used to iterate
int VariableResolvers::next_var(int curve_id, Date *date, int var) const
{
	if (curve_id < 0 || curve_id >= MAX_CURVES_PER_RESOLVER)
		return datasize_;
	PricingCurve curve = curves_[curve_id];
	if (curve.id() == 0)
		return datasize_;
	for (int i = var; i < datasize_; i++) {
		if (!data_[i].is_null() && data_[i].curveid == curve_id) {
			*date = data_[i].date;
			return i;
		}
	}
	return datasize_;
}

VariableIterator VariableResolver::begin() const
{
	Date d = 0;
	int var = resolvers_->next_var(curve_, &d, 0);
	return VariableIterator(resolvers_, curve_, var, d);
}

VariableIterator VariableResolver::end() const
{
	return VariableIterator(resolvers_, curve_, resolvers_->capacity(), 0);
}

VariableIterator &VariableIterator::operator++()
{
	if (pos_ < resolvers_->capacity()) {
		pos_ = resolvers_->next_var(curve_, &date_, pos_ + 1);
	}
	return *this;
}

// Adds a var and returns the var id
// if var id is -1 then error occurred - there was
// no more space available
int VariableResolvers::add_var(int curve_id, Date date)
{
	int var = 0;
	PricingCurve curve;
	int result = search_var(curve_id, date, &curve, &var);
	if (result == 1)
		return var;
	if (result < 0)
		return -1;
	data_[var].date = date;
	data_[var].curveid = curve_id;
	return var;
}

// Add or get the variable assignment for date
int VariableResolver::add_var(Date date) { return resolvers_->add_var(curve_, date); }

// Retrieve the variable assignment for a date, -1 if no assignment found
int VariableResolver::get_var(Date date) const
{
	int var;
	PricingCurve curve;
	if (resolvers_->search_var(curve_, date, &curve, &var) == 1)
		return var;
	return -1;
}

FloatingCalculation::FloatingCalculation() : index(nullptr), index_tau(0.0), value_date_var(-1), maturity_date_var(-1)
{
}

FloatingPeriod::FloatingPeriod()
    : next(nullptr), notional(0.0), period_tau(0.0), spread(0.0), weight(0.0), distance(0.0)
{
}

NFloatingCashflow::NFloatingCashflow()
    : head(nullptr), tail(nullptr), var(-1), compounding_method(CompoundingMethod::CM_NONE)
{
}
void NFloatingCashflow::push_back(FloatingPeriod *fp)
{
	assert(fp != nullptr);
	if (!head)
		head = fp;
	if (tail)
		tail->next = fp;
	tail = fp;
}
void NFloatingCashflow::set_paydate(Date dt, int variable)
{
	assert(is_valid_date(dt));
	paydate = dt;
	var = variable;
}

// Retrieve a fixing.
// If the fixing date is < evaluation date then the absence of a fixing
// will be an error reported via status. If the fixing date is == evaluation
// date
// then a missing fixing is not treated as error - instead the method will
// return false; in all cases a true return value indicates that
// the fixing was found and is set
bool ValuationContextImpl::get_fixing(IndexId fixing_key, Date fixingDt, double &fixing, StatusCode &status) const
{
	fixing = 0.0;
	bool lookup_fixing = (fixingDt < evaluation_date_ || (fixingDt == evaluation_date_ && include_todays_fixing_));
	if (lookup_fixing && fixing_data_service_ == nullptr) {
		status = StatusCode::kVAL_FixingDataUnavailable;
		return false;
	}

	if (lookup_fixing) {
		auto ts = fixing_data_service_->get_fixings(fixing_key);
		if (ts == nullptr) {
			status = StatusCode::kVAL_FixingDataUnavailable;
			return false;
		}
		if (ts->find(fixingDt, fixing)) {
			//      fprintf(stderr, "fixing for %s %d-%d-%d =
			//      %.6f\n", fixing_key,
			//              fixingDt.year(), fixingDt.month(),
			//              fixingDt.day_of_month(),
			//              fixing_value->value);
			return true;
		}
		// fixing not found - if the fixing date is
		// in the past we are in trouble
		if (fixingDt < evaluation_date_) {
			status = StatusCode::kVAL_FixingNotFound;
			auto ymd = date_components(fixingDt);
			error("Fixing for %d-%d-%d was not found\n", ymd.y, (int)ymd.m, (int)ymd.d);
			return false;
		}
	}
	return false;
}

PricingCurve SimpleCurveMapper::map_index_tenor(PricingCurveType curve_type, Currency currency, IndexFamily family,
						Tenor tenor) const
{
	PricingCurve x(curve_type, currency, family, tenor);
	auto iter = mappings_.find(x);
	if (iter != mappings_.end())
		return iter->second;
	if (family != INDEX_FAMILY_UNSPECIFIED && tenor != TENOR_UNSPECIFIED) {
		PricingCurve x2(curve_type, currency, family);
		iter = mappings_.find(x2);
		if (iter != mappings_.end())
			return iter->second;
	}
	if (family != INDEX_FAMILY_UNSPECIFIED) {
		PricingCurve x3(curve_type, currency);
		iter = mappings_.find(x3);
		if (iter != mappings_.end())
			return iter->second;
	}
	return x;
}

void SimpleCurveMapper::add_mapping(PricingCurveType curve_type, Currency currency, IndexFamily family, Tenor tenor,
				    PricingCurve curve)
{
	PricingCurve x(curve_type, currency, family, tenor);
	add_mapping(x, curve);
}

void SimpleCurveMapper::add_mapping(PricingCurve from_curve, PricingCurve to_curve)
{
	// fprintf(stdout, "Added curve mapping from %s to %s\n",
	//        from_curve.name().c_str(), to_curve.name().c_str());
	mappings_.insert({from_curve, to_curve});
}

Cashflow *create_simple_cashflow(RegionAllocator *A, const ValuationContext &ctx, VariableResolvers *resolvers,
				 Currency currency, double amount, Date payment_date, IsdaIndex trade_index,
				 const CurveMapper *curve_mapper)
{
	IndexFamily index_family = INDEX_FAMILY_UNSPECIFIED;
	if (trade_index != ISDA_INDEX_UNSPECIFIED) {
		auto index = get_default_index_service()->get_index(trade_index, TENOR_UNSPECIFIED);
		if (!index) {
			error("Index definition not found for %s\n",
			      get_default_converter()->isda_index_to_string(trade_index));
			return nullptr;
		}
		index_family = index->family();
	}
	Cashflow *cashflow = new (*A) Cashflow(resolvers);
	NSimpleCashflow *simple_cashflow = new (*A) NSimpleCashflow();
	cashflow->set(simple_cashflow);
	simple_cashflow->amount = amount;
	simple_cashflow->paydate = payment_date;
	PricingCurve discount =
	    curve_mapper->map_index_tenor(PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT, currency, index_family);
	simple_cashflow->var = -1;
	if (resolvers) {
		auto resolver = resolvers->add_resolver(discount);
		if (resolvers->capacity())
			simple_cashflow->var = resolver.add_var(simple_cashflow->paydate);
	}
	return cashflow;
}

Cashflow *create_simple_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFSimple *cf,
				 VariableResolvers *resolvers, const CurveMapper *curve_mapper)
{
	return create_simple_cashflow(A, ctx, resolvers, cf->currency(), cf->amount(), cf->payment_date(),
				      cf->trade_index(), curve_mapper);
}

static NOISCashflow *create_ois_cashflow(RegionAllocator *A, const ValuationContext &ctx, VariableResolvers *resolvers,
					 IsdaIndex isda_index, Date accrual_start, Date accrual_end, Date paymentDate,
					 DayCountFraction day_count_fraction, double notional,
					 const CurveMapper *curve_mapper)
{
	const DayFraction *fraction = get_day_fraction(day_count_fraction);
	if (!fraction) {
		error("unexpected error: day fraction is null\n");
		return nullptr;
	}
	auto index = get_default_index_service()->get_index(isda_index, TENOR_1D);
	if (!index) {
		error("Index definition not found for %s\n", get_default_converter()->isda_index_to_string(isda_index));
		return nullptr;
	}

	NOISCashflow *flow = new (*A) NOISCashflow();

	// first fixing date
	auto fixing_date = index->fixing_date(accrual_start);
	// and corresponding value date
	auto value_date = accrual_start;
	auto fixing_key = index->id();
	auto fixing_calendar = index->fixing_calendar();

	double compoundFactor = 1.0;
	Period oneday(1, PeriodUnit::DAYS);

	// we need to accrue interest for past fixings
	StatusCode status = StatusCode::kOk;
	while (fixing_date <= ctx.evaluation_date()) {
		double fixing = 0.0;
		if (!ctx.get_fixing(fixing_key, fixing_date, fixing, status))
			break;
		auto next_value_date = fixing_calendar->advance(value_date, oneday);
		auto accrual_period = fraction->year_fraction(value_date, next_value_date);
		compoundFactor *= (1 + fixing * accrual_period);
		value_date = next_value_date;
		fixing_date = fixing_calendar->advance(fixing_date, oneday);
	}
	if (status != StatusCode::kOk) {
		error("unexpected error: fixings not found\n");
		return nullptr;
	}

	flow->compound_factor = compoundFactor;
	flow->fixing_date = fixing_date;
	flow->value_date = value_date;
	flow->maturity_date = accrual_end;
	flow->index = index;
	flow->period_tau = fraction->year_fraction(accrual_start, accrual_end);
	flow->tele_tau = fraction->year_fraction(value_date, accrual_end);
	flow->notional = notional;
	// each OIS cash flow contributes one variable on
	// discount curve
	flow->paydate = paymentDate;
	flow->value_date_var = -1;
	flow->maturity_date_var = -1;
	flow->var = -1;

	if (resolvers) {
		// Need a logical name for the forward curve and
		// discount curve as the curves are not known at this stage
		PricingCurve forward_curve_name = curve_mapper->map_index_tenor(
		    PricingCurveType::PRICING_CURVE_TYPE_FORWARD, index->currency(), index->family(), index->tenor());
		PricingCurve discount_curve_name = curve_mapper->map_index_tenor(
		    PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT, index->currency(), index->family());

		auto resolver_discount = resolvers->add_resolver(discount_curve_name);
		auto resolver_forward = resolvers->add_resolver(forward_curve_name);
		if (resolvers->capacity() && value_date < accrual_end) {
			flow->value_date_var = resolver_forward.add_var(value_date);
			flow->maturity_date_var = resolver_forward.add_var(accrual_end);
		}
		if (resolvers->capacity())
			flow->var = resolver_discount.add_var(paymentDate);
	}
	return flow;
}

Cashflow *create_ois_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFOis *ois,
			      VariableResolvers *resolvers, const CurveMapper *curve_mapper)
{
	Cashflow *cashflow = new (*A) Cashflow(resolvers);
	NOISCashflow *oisflow =
	    create_ois_cashflow(A, ctx, resolvers, ois->index(), ois->accrual_start_date(), ois->accrual_end_date(),
				ois->payment_date(), ois->day_count_fraction(), ois->notional(), curve_mapper);
	if (!oisflow)
		return nullptr;
	cashflow->set(oisflow);
	return cashflow;
}

static FloatingPeriod *create_floating_period(RegionAllocator *A, const ValuationContext &ctx,
					      VariableResolvers *resolvers, Date accrual_start, Date accrual_end,
					      IsdaIndex isda_index, Tenor tenor, IsdaIndex isda_index2, Tenor tenor2,
					      DayCountFraction day_count_fraction, double notional, double spread,
					      const CurveMapper *curve_mapper)
{
	const DayFraction *fraction = get_day_fraction(day_count_fraction);
	if (fraction == nullptr) {
		error("unexpected error: day fraction is null\n");
		return nullptr;
	}
	double year_fraction = fraction->year_fraction(accrual_start, accrual_end);

	auto index = get_default_index_service()->get_index(isda_index, tenor);
	auto fixing_dt = index->fixing_date(accrual_start);
	auto value_dt = index->value_date(fixing_dt);
	auto maturity_dt = index->maturity_date(value_dt);
	double index_tau = fraction->year_fraction(value_dt, maturity_dt);
	// Need a logical name for the forward curve and
	// discount curve as the curves are not known at this stage
	PricingCurve forward1 = curve_mapper->map_index_tenor(PricingCurveType::PRICING_CURVE_TYPE_FORWARD,
							      index->currency(), index->family(), index->tenor());

	FloatingPeriod *fp = new (*A) FloatingPeriod();
	fp->float1.index = index;
	fp->float1.index_tau = index_tau;
	fp->float1.fixing_date = fixing_dt;
	fp->float1.value_date = value_dt;
	fp->float1.maturity_date = maturity_dt;
	fp->float1.value_date_var = -1;
	fp->float1.maturity_date_var = -1;
	if (resolvers) {
		auto resolver_forward1 = resolvers->add_resolver(forward1);
		if (resolvers->capacity()) {
			fp->float1.value_date_var = resolver_forward1.add_var(value_dt);
			fp->float1.maturity_date_var = resolver_forward1.add_var(maturity_dt);
		}
	}
	if (isda_index2 != IsdaIndex::ISDA_INDEX_UNSPECIFIED && tenor2 != TENOR_UNSPECIFIED) {
		// second floating rate may exist in stubs
		if (isda_index != isda_index2) {
			error("Both tenors must have the same index\n");
			return nullptr;
		}
		auto index2 = get_default_index_service()->get_index(isda_index2, tenor2);

		PricingCurve forward2 =
		    curve_mapper->map_index_tenor(PricingCurveType::PRICING_CURVE_TYPE_FORWARD, index2->currency(),
						  index2->family(), index2->tenor());

		auto fixing_dt2 = index2->fixing_date(accrual_start);
		auto value_dt2 = index2->value_date(fixing_dt2);
		auto maturity_dt2 = index2->maturity_date(value_dt2);
		double index_tau2 = fraction->year_fraction(value_dt2, maturity_dt2);

		auto end = accrual_end;
		auto first_maturity = maturity_dt;
		auto second_maturity = maturity_dt2;
		auto distance_first_second = second_maturity - first_maturity;
		// when two floating rates are used we need to
		// assign linear weight to each rate
		double weight1 = (double(end) - double(first_maturity)) / double(distance_first_second);

		fp->weight = weight1;
		fp->distance = distance_first_second;

		fp->float2.index = index2;
		fp->float2.index_tau = index_tau2;
		fp->float2.fixing_date = fixing_dt2;
		fp->float2.value_date = value_dt2;
		fp->float2.maturity_date = maturity_dt2;

		fp->float2.value_date_var = -1;
		fp->float2.maturity_date_var = -1;
		if (resolvers) {
			auto resolver_forward2 = resolvers->add_resolver(forward2);
			if (resolvers->capacity()) {
				fp->float2.value_date_var = resolver_forward2.add_var(value_dt2);
				fp->float2.maturity_date_var = resolver_forward2.add_var(maturity_dt2);
			}
		}
	}
	fp->notional = notional;
	fp->period_tau = year_fraction;
	fp->spread = spread;
	return fp;
}

/**
 * Converts a CFFloating to internal format suitable for valuation
 */
Cashflow *create_floating_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFFloating *in,
				   VariableResolvers *resolvers, const CurveMapper *curve_mapper)
{
	assert(in->floating_periods_size() > 0);
	auto index =
	    get_default_index_service()->get_index(in->floating_periods(0).index(), in->floating_periods(0).tenor());
	if (!index) {
		error("Index definition not found for %s\n",
		      get_default_converter()->isda_index_to_string(in->floating_periods(0).index()));
		return nullptr;
	}
	IndexFamily index_family = index->family();
	PricingCurve discount_curve =
	    curve_mapper->map_index_tenor(PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT, in->currency(), index_family);
	Cashflow *cashflow = new (*A) Cashflow(resolvers);
	NFloatingCashflow *flow = new (*A) NFloatingCashflow();
	cashflow->set(flow);
	for (int i = 0; i < in->floating_periods_size(); i++) {
		const CFFloatingPeriod &p = in->floating_periods(i);
		auto fp = create_floating_period(A, ctx, resolvers, p.accrual_start_date(), p.accrual_end_date(),
						 p.index(), p.tenor(), p.index2(), p.tenor2(), in->day_count_fraction(),
						 p.notional(), p.spread(), curve_mapper);
		if (!fp)
			return nullptr;
		flow->push_back(fp);
	}
	flow->paydate = in->payment_date();
	flow->var = -1;
	if (resolvers) {
		auto resolver = resolvers->add_resolver(discount_curve);
		if (resolvers->capacity())
			flow->var = resolver.add_var(flow->paydate);
	}
	flow->compounding_method = in->compounding_method();
	return cashflow;
}

// Count variables
static int count_vars_in_floating_cashflows(const CFFloating *in)
{
	int n = 0;
	for (int i = 0; i < in->floating_periods_size(); i++) {
		const CFFloatingPeriod &p = in->floating_periods(i);
		// If a floating period has two different tenors then
		// we need to allow for 2 extra vars
		if (p.index2() != IsdaIndex::ISDA_INDEX_UNSPECIFIED && p.tenor2() != Tenor::TENOR_UNSPECIFIED) {
			// value and maturity date for second tenor
			n += 2;
		}
		// value and maturity date for first tenor
		n += 2;
	}
	// payment date
	n += 1;
	return n;
}

Cashflow *create_fra_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFFra *in,
			      VariableResolvers *resolvers, const CurveMapper *curve_mapper)
{
	auto index =
	    get_default_index_service()->get_index(in->floating_period().index(), in->floating_period().tenor());
	if (!index) {
		error("Index definition not found for %s\n",
		      get_default_converter()->isda_index_to_string(in->floating_period().index()));
		return nullptr;
	}
	IndexFamily index_family = index->family();
	Cashflow *cashflow = new (*A) Cashflow(resolvers);
	NFraCashflow *flow = new (*A) NFraCashflow();
	cashflow->set(flow);
	flow->floating_period = create_floating_period(
	    A, ctx, resolvers, in->floating_period().accrual_start_date(), in->floating_period().accrual_end_date(),
	    in->floating_period().index(), in->floating_period().tenor(), in->floating_period().index2(),
	    in->floating_period().tenor2(), in->day_count_fraction(), in->floating_period().notional(), 0.0,
	    curve_mapper);
	if (!flow->floating_period)
		return nullptr;
	flow->fixed_rate = in->fixed_rate();
	if (in->currency() == Currency::AUD || in->currency() == Currency::NZD) {
		flow->double_discounting_factor = 1.0 / (1.0 + flow->fixed_rate * flow->floating_period->period_tau);
	} else {
		flow->double_discounting_factor = 1.0;
	}
	flow->paydate = in->payment_date();
	flow->var = -1;
	if (resolvers) {
		PricingCurve discount_curve = curve_mapper->map_index_tenor(
		    PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT, in->currency(), index_family);
		auto resolver_discount = resolvers->add_resolver(discount_curve);
		if (resolvers->capacity())
			flow->var = resolver_discount.add_var(in->payment_date());
	}
	return cashflow;
}

// Count variables
static int count_vars_in_fra_cashflow(const CFFra *in)
{
	int n = 0;
	const CFFloatingPeriod &p = in->floating_period();
	// If a floating period has two different tenors then
	// we need to allow for 2 extra vars
	if (p.index2() != IsdaIndex::ISDA_INDEX_UNSPECIFIED && p.tenor2() != Tenor::TENOR_UNSPECIFIED) {
		// value and maturity date for second tenor
		n += 2;
	}
	// value and maturity date for first tenor
	n += 2;
	// payment date
	n += 1;
	return n;
}

void CashflowStream::add(Cashflow *flow)
{
	if (!first)
		first = flow;
	if (last)
		last->next = flow;
	last = flow;
}

void Cashflows::add(CashflowStream *stream)
{
	if (!first_)
		first_ = stream;
	if (last_)
		last_->next = stream;
	last_ = stream;
	n_streams_++;
}

/**
 * New version
 * Converts the CFCollection to internal cashflow format
 * Our aim here to allocate the memory in such a way that it is
 * easier to keep data together - thereby improvng cache locality
 * and (hopefully) performance. We also do not allocate the
 * Resolvers if derivaties are not requested - this should further
 * optimise the memory layout.
 */
Cashflows *construct_cashflows(RegionAllocator *A, const CFCollection *cfcollection, const ValuationContext &ctx,
			       const CurveMapper *curve_mapper)
{
	Cashflows *flows = new (*A) Cashflows();

	for (int i = 0; i < cfcollection->streams_size(); i++) {
		CashflowStream *leg = new (*A) CashflowStream();
		flows->add(leg);
		auto &stream = cfcollection->streams(i);
		leg->factor = stream.factor() != 0.0 ? stream.factor() : 1.0;
		for (int j = 0; j < stream.cashflows_size(); j++) {
			auto &cf = stream.cashflows(j);
			// Below we only setup vars if dervivatives are needed
			// However we still need VariableResolvers populated to
			// hold the logical curves used for mapping to real
			// curves
			int nvars = 0;
			VariableResolvers *resolvers = nullptr;
			if (cf.has_floating()) {
				if (ctx.derivative_order() > 0) {
					nvars = count_vars_in_floating_cashflows(&cf.floating());
				}
				resolvers = new (*A) VariableResolvers(A, nvars);
				auto flow = create_floating_cashflow(A, ctx, &cf.floating(), resolvers, curve_mapper);
				if (!flow)
					return nullptr;
				leg->add(flow);
			} else if (cf.has_ois()) {
				if (ctx.derivative_order() > 0) {
					nvars = 3;
				}
				resolvers = new (*A) VariableResolvers(A, nvars);
				auto flow = create_ois_cashflow(A, ctx, &cf.ois(), resolvers, curve_mapper);
				if (!flow) {
					return nullptr;
				}
				leg->add(flow);
			} else if (cf.has_simple()) {
				if (ctx.derivative_order() > 0) {
					nvars = 1;
				}
				resolvers = new (*A) VariableResolvers(A, nvars);
				auto flow = create_simple_cashflow(A, ctx, &cf.simple(), resolvers, curve_mapper);
				if (!flow) {
					return nullptr;
				}
				leg->add(flow);
			} else if (cf.has_fra()) {
				if (ctx.derivative_order() > 0) {
					nvars = count_vars_in_fra_cashflow(&cf.fra());
				}
				resolvers = new (*A) VariableResolvers(A, nvars);
				auto flow = create_fra_cashflow(A, ctx, &cf.fra(), resolvers, curve_mapper);
				if (!flow)
					return nullptr;
				leg->add(flow);
			} else {
				return nullptr;
			}
		}
	}
	return flows;
}

} // namespace redukti
