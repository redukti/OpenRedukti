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
#ifndef _REDUKTI_CASHFLOW_INTERNAL_H_
#define _REDUKTI_CASHFLOW_INTERNAL_H_

#include <cashflow.h>

#include <algorithm>
#include <array>
#include <map>
#include <memory>

namespace redukti
{

class ValuationContextImpl : public ValuationContext
{
	private:
	// Business date
	Date evaluation_date_;
	// required derivative order
	// allowed values are 0, 1 or 2
	int derivative_order_;
	// If a coupon is paying before this date then
	// assume it is already paid therefore value of coupon is 0
	Date payment_cutoff_date_;
	// handle to fixings data
	FixingDataService *fixing_data_service_;
	// include today's fixing (e.g. eod)
	// If false then the curve will be used to determine the
	// rate. The bootstrapper requires this to be false so we set
	// the default value to false
	bool include_todays_fixing_;

	public:
	ValuationContextImpl(Date edate, int order = 0, FixingDataService *fixing_service = nullptr,
			     Date payment_cutoff = 0, bool include_todays_fixing = false)
	    : evaluation_date_(edate), derivative_order_(order),
	      payment_cutoff_date_(payment_cutoff == 0 ? edate : payment_cutoff), fixing_data_service_(fixing_service),
	      include_todays_fixing_(include_todays_fixing)
	{
	}
	Date evaluation_date() const override final { return evaluation_date_; }
	Date payment_cutoff_date() const override final { return payment_cutoff_date_; }
	int derivative_order() const override final { return derivative_order_; }
	bool include_todays_fixing() const override final { return include_todays_fixing_; }
	void set_derivative_order(int order) { derivative_order_ = order; }
	bool get_fixing(IndexId fixing_key, Date fixingDt, double &fixing, StatusCode &status) const override final;
};

class SimpleCurveMapper : public CurveMapper
{
	private:
	std::map<PricingCurve, PricingCurve> mappings_;

	public:
	PricingCurve map_index_tenor(PricingCurveType curve_type, Currency currency, IndexFamily family,
				     Tenor tenor) const override final;
	void add_mapping(PricingCurveType curve_type, Currency currency, IndexFamily family, Tenor tenor,
			 PricingCurve curve);
	void add_mapping(PricingCurve from_curve, PricingCurve to_curve);
	void clear_mappings() { mappings_.clear(); }
};

struct VariableResolvers;
struct VariableIterator {
	VariableResolvers *resolvers_;
	int32_t curve_;
	int32_t pos_;
	Date date_;

	public:
	VariableIterator(VariableResolvers *parent, int32_t curve, int32_t pos, Date d)
	    : resolvers_(parent), curve_(curve), pos_(pos), date_(d)
	{
	}
	VariableIterator &operator++();
	Date date() const { return date_; }
	int var() const;
};

inline bool operator==(const VariableIterator &a, const VariableIterator &b)
{
	return a.resolvers_ == b.resolvers_ && a.curve_ == b.curve_ && a.pos_ == b.pos_;
}

inline bool operator!=(const VariableIterator &a, const VariableIterator &b) { return !(a == b); }

struct VariableResolver {
	private:
	VariableResolvers *resolvers_;
	// Position of curve in the VariableResolvers.curves_
	// If -1 then not valid
	int32_t curve_;

	public:
	VariableResolver(VariableResolvers *owner, int32_t id) : resolvers_(owner), curve_(id) {}
	PricingCurve curve() const;
	// Add or get the variable assignment for date
	int add_var(Date date);
	// Retrieve the variable assignment for a date, -1 if no assignment
	// found
	int get_var(Date date) const;

	VariableIterator begin() const;
	VariableIterator end() const;
};

struct DateVar {
	unsigned date : 24, curveid : 8;

	bool is_null() const { return date == 0; }
};

struct VariableResolvers {
	enum { MAX_CURVES_PER_RESOLVER = 3 };

	private:
	// Maximum of 3 pricing curves per resolver
	PricingCurve curves_[MAX_CURVES_PER_RESOLVER];
	// Size of the data array below
	uint32_t datasize_;
	// Data is stored as an array
	// Each entry is made up of:
	// curve id - 8 bits
	// date value - 32 bits
	DateVar *data_;

	public:
	VariableResolvers(Allocator *A, uint32_t datasize) : datasize_(datasize)
	{
		data_ = new (*A) DateVar[datasize_];
		reset();
	}

	// Search and add a curve to the resolvers set
	// If successful then the returned resolver will
	// have give curve's location
	// A value of -1 indicates failure
	VariableResolver add_resolver(PricingCurve curve);

	// search for a var resolver by curve
	// If not found returns resolver with curve = -1
	VariableResolver get_resolver(PricingCurve curve);

	// Gets resolver at position
	// If invalid then returns resolver with curve = -1
	VariableResolver get_resolver_at(int32_t pos);

	// How many resolvers, i.e. curves
	size_t num_resolvers() const;

	// how many variables do we have?
	int num_vars() const;

	// Obtain a registered curve of specified type
	// Optional index specifies which relative curve to get
	// E.g. whether to get the first (index = 0) forward curve
	// or second (index = 1) - note that first or second do not
	// have any meaning other than that may return different
	// curves
	PricingCurve get_curve_of_type(PricingCurveType type, int index = 0) const;

	// Initialize
	void reset();

	// Adds a var and returns the var id
	// if var id is -1 then error occurred - there was
	// no more space available
	int add_var(int curve_id, Date date);

	// Searches for a var by curve id and date
	// if found returns 1, and *var is set
	// If error (no more space) return -1, and *var is not set
	// If not found returns 0 and *var is set to the next available slot
	// This interface  allows this to be reused in add_var()
	int search_var(int curve_id, Date date, PricingCurve *curve, int *var) const;

	// Searches for next date/var for given curve
	// search starts at var - set to 0 for first call
	// returns next var which can be used to iterate
	// by passig var+1 in the next call
	// This method is primarily used by the iterator implementation
	int next_var(int curve_id, Date *date, int var) const;

	// Gets the curve at specified position
	PricingCurve get_curve(int32_t id) const;

	int capacity() const { return datasize_; }

	private:
	VariableResolvers(const VariableResolvers &) = delete;
	VariableResolvers &operator=(const VariableResolvers &) = delete;
};

inline int VariableIterator::var() const { return pos_ >= resolvers_->capacity() ? -1 : pos_; }

// SimpleCashflow just has an amount and a payment date
// The payment date is a variable on the discount curve
struct NSimpleCashflow {
	NSimpleCashflow() : amount(0), var(0), paydate(0) {}

	// Scalar amount
	// If fixed rate then notional * fixed rate * day fraction
	// If known amount then amount
	// We do not retain the details (saves memory, details can be worked
	// out)
	double amount;
	// variable number - used to calculate derivatives
	int var;
	// payment date
	Date paydate;
	// Is the party paying or receiving
	bool pay;

	private:
	NSimpleCashflow(const NSimpleCashflow &) = delete;
	NSimpleCashflow &operator=(const NSimpleCashflow &) = delete;
};

struct FloatingCalculation {
	// index - this also identifies the tenor etc.
	const InterestRateIndex *index;
	// day count as per index
	double index_tau;
	// fixing date
	Date fixing_date;
	// value date
	Date value_date;
	// variable for value date
	int value_date_var;
	// maturity date
	Date maturity_date;
	// variable for maturity date
	int maturity_date_var;

	FloatingCalculation();

	private:
	FloatingCalculation(const FloatingCalculation &) = delete;
	FloatingCalculation &operator=(const FloatingCalculation &) = delete;
};

// Holds data for a single floating rate
// calculation period
struct FloatingPeriod {
	FloatingPeriod *next;
	double notional;
	// accrual period year fraction
	double period_tau;
	// index period year fraction
	// spread if any
	double spread;
	// weight & distance for interpolating
	double weight;
	double distance;
	// First floating calculation (required)
	FloatingCalculation float1;
	// Second floating calculation (optional - used in stubs)
	FloatingCalculation float2;

	FloatingPeriod();

	private:
	FloatingPeriod(const FloatingPeriod &) = delete;
	FloatingPeriod &operator=(const FloatingPeriod &) = delete;
};

struct NFloatingCashflow {
	// First floating period in this cash flow
	FloatingPeriod *head;
	// Last floating period in this cash flow
	FloatingPeriod *tail;
	Date paydate;
	// Variable on discount curve for payment date
	int var;
	// F for Flat or S for Straight or 0
	CompoundingMethod compounding_method;

	NFloatingCashflow();
	void push_back(FloatingPeriod *fp);
	void set_paydate(Date dt, int variable);

	private:
	NFloatingCashflow(const NFloatingCashflow &) = delete;
	NFloatingCashflow &operator=(const NFloatingCashflow &) = delete;
};

struct NOISCashflow {
	// notional amount
	double notional;
	// index - this also identifies the tenor etc.
	const InterestRateIndex *index;
	// day count as per index
	double period_tau;
	// telecopic period
	double tele_tau;
	// fixing date
	Date fixing_date;
	// value date
	Date value_date;
	// variable for value date
	int value_date_var;
	// maturity date
	Date maturity_date;
	// variable for maturity date
	int maturity_date_var;
	// compounded interest for past fixings
	double compound_factor;
	// payment date
	Date paydate;
	// Variable on discount curve for payment date
	int var;

	// methods
	void set_paydate(Date dt, int variable);
};

struct NFraCashflow {
	double fixed_rate;
	FloatingPeriod *floating_period;
	// payment date
	Date paydate;
	// Variable on discount curve for payment date
	int var;
	// For AUD and NZD curencies
	// This is
	double double_discounting_factor;

	NFraCashflow() : fixed_rate(0), floating_period(nullptr), paydate(0), var(-1), double_discounting_factor(0) {}
};

enum class CashflowType { CF_NONE, CF_SIMPLE, CF_FLOATING, CF_OIS, CF_FRA };

// the Cashflow contains a union so that the PV templates can
// identify the correct PV function to execute
struct Cashflow {
	CashflowType type;

	union {
		NSimpleCashflow *simple;
		NFloatingCashflow *floating;
		NOISCashflow *ois;
		NFraCashflow *fra;
	};

	VariableResolvers *resolvers;
	Cashflow *next;

	Cashflow(VariableResolvers *r) : type(CashflowType::CF_NONE), simple(nullptr), resolvers(r), next(nullptr) {}
	void set(NSimpleCashflow *s)
	{
		simple = s;
		type = CashflowType::CF_SIMPLE;
	}
	void set(NFloatingCashflow *f)
	{
		floating = f;
		type = CashflowType::CF_FLOATING;
	}
	void set(NOISCashflow *f)
	{
		ois = f;
		type = CashflowType::CF_OIS;
	}
	void set(NFraCashflow *f)
	{
		fra = f;
		type = CashflowType::CF_FRA;
	}

	private:
	Cashflow(const Cashflow &) = delete;
	Cashflow &operator=(const Cashflow &) = delete;
};

// Leg is a collection of cashflow objects
struct CashflowStream {
	double factor;
	Cashflow *first;
	Cashflow *last;
	CashflowStream *next;

	CashflowStream() : factor(1.0), first(nullptr), last(nullptr), next(nullptr) {}
	void add(Cashflow *flow);

	private:
	CashflowStream(const CashflowStream &) = delete;
	CashflowStream &operator=(const CashflowStream &) = delete;
};

// Cashflows is a collection of 2 legs
class Cashflows
{
	public:
	int n_streams_;
	CashflowStream *first_;
	CashflowStream *last_;

	Cashflows() : n_streams_(0), first_(nullptr), last_(nullptr) {}
	void add(CashflowStream *stream);
	int n_streams() const { return n_streams_; }

	private:
	Cashflows(const Cashflows &) = delete;
	Cashflows &operator=(const Cashflows &) = delete;
};
class CFSimple;
class CFOis;
class CFFloating;
class CFFra;
extern Cashflow *create_simple_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFSimple *cf,
					VariableResolvers *resolvers, const CurveMapper *curve_mapper);

extern Cashflow *create_ois_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFOis *ois,
				     VariableResolvers *resolvers, const CurveMapper *curve_mapper);

extern Cashflow *create_floating_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFFloating *in,
					  VariableResolvers *resolvers, const CurveMapper *curve_mapper);

extern Cashflow *create_fra_cashflow(RegionAllocator *A, const ValuationContext &ctx, const CFFra *in,
				     VariableResolvers *resolvers, const CurveMapper *curve_mapper);

} // namespace redukti

#endif
