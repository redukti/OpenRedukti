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

#include <bootstrap.pb.h>
#include <curve.pb.h>
#include <enums.pb.h>

#include <bootstrap.h>
#include <cashflow.h>
#include <cashflow_pricing.h>
#include <converters.h>
#include <curve.h>
#include <logger.h>
#include <matrix.h>
#include <status.h>

#include <internal/cashflow_internal.h>
#include <internal/cashflow_pricing_internal.h>
#include <internal/raviapi_internal.h>

#include <cminpack.h>

#include <array>
#include <map>
#include <memory>
#include <string>

#include <cmath>

// TODO
// Long term - look at using autodiff

// NOTE
// Following error occurs when the data is uninitialized or invalid and it is
// submitted to the solver
// Causes a crash
// DEBUG tid(480) C:\d\redukti\Redukti-Risk\src\bootstrap.cpp:1372
// (handle_bootstrap_request) Starting curve build
// ** On entry to DGECON parameter number  4 had an illegal value

using namespace google::protobuf;

namespace redukti
{

// Given the CurveDefinitionId locate the definition
class CurveDefinitionProvider
{
	public:
	virtual ~CurveDefinitionProvider() = default;
	virtual const IRCurveDefinition *get_definition_by_id(int id) const = 0;
};

// Simple implementation of id to definition mapping
class CurveDefinitionProviderImpl : public CurveDefinitionProvider
{
	private:
	std::map<int, const IRCurveDefinition *> mapping;

	public:
	~CurveDefinitionProviderImpl() = default;
	void add(const IRCurveDefinition *def) { mapping.insert({(int)def->id(), def}); }
	const IRCurveDefinition *get_definition_by_id(int id) const final
	{
		auto &&iter = mapping.find(id);
		if (iter != mapping.end())
			return iter->second;
		return nullptr;
	}
};

// When bootstrapping we fix the curves for each pricing call
// Based on instrument definition - the curves are referenced by
// their respective definition ids
// Note objects of this class are copyable
class BootstrapCurveMapper : public CurveMapper
{
	private:
	PricingCurve forward_curve_;
	PricingCurve discount_curve_;
	// The following reference the curves by position
	int forward_curve_index_;
	int discount_curve_index_;

	public:
	BootstrapCurveMapper() : forward_curve_index_(-1), discount_curve_index_(-1) {}
	BootstrapCurveMapper(PricingCurve fc, int fc_idx, PricingCurve dc, int dc_idx)
	    : forward_curve_(fc), discount_curve_(dc), forward_curve_index_(fc_idx), discount_curve_index_(dc_idx)
	{
	}
	PricingCurve map_index_tenor(PricingCurveType curve_type, Currency currency, IndexFamily family,
				     Tenor tenor) const override final
	{
		assert(forward_curve_index_ != -1 && discount_curve_index_ != -1);
		return curve_type == PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT ? discount_curve_ : forward_curve_;
	}
	// Returns the forward curve index
	int forward_curve() const
	{
		assert(forward_curve_index_ != -1);
		return forward_curve_index_;
	}
	// Returns the discount curve index
	int discount_curve() const
	{
		assert(forward_curve_index_ != -1);
		return discount_curve_index_;
	}
};

// For bootstrapping we fix the curves to be used, as all instruments are
// vanilla and do not have stubs etc. we only need two curves
class BootstrapCurveProvider : public CurveProvider
{
	private:
	CurveReference *discount_curve_;
	CurveReference *forward_curve_;

	public:
	BootstrapCurveProvider(CurveReference *dc, CurveReference *fc) : discount_curve_(dc), forward_curve_(fc) {}
	const CurveReference *get_curve(PricingCurve curve) const
	{
		if (curve.curve_type() == PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT)
			return discount_curve_;
		else
			return forward_curve_;
	}
};

// CurveHolder contains the curve being built during
// bootstrapping. In addition to the base curve it is also
// used to bump pillar points on the curve so that the Jacobian
// can be computed numerically. There is a current curve
// pointer which can be switched from base curve to bumped curve;
// this allows the pricing functions to be unaware of which
// curve is being used
struct CurveHolder : public CurveReference {
	CurveType curve_type_;
	uint32_t n_maturities_;
	Date as_of_date_;
	const DayFraction *fraction_;
	IRRateType interpolated_on_;
	Date *maturities_;
	// Base zero rates - computed by bootstrap process
	double *rates_;
	// Temp curve populated by bumping a pillar
	double *rates_plus_h_;
	// Temp curve populated by bumping a pillar
	double *rates_minus_h_;
	// Copy of the bootstrapped rates for
	// used during construction of the PAR delta
	// Jacobian
	double *saved_rates_;

	// Curves may be interpolated on discount factors rather than rates
	// values hold the data that is interpolated, so they can be either
	// zero rates or discount factors
	double *values_;
	double *values_plus_h_;
	double *values_minus_h_;
	double *saved_values_;

	std::unique_ptr<YieldCurve, Deleter<YieldCurve>> base_curve_;
	std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve_plus_h_;
	std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve_minus_h_;

	// The curve used for pricing - could be one of the three above
	YieldCurve *pricing_curve_;

	CurveHolder(CurveType curve_type, uint32_t n_maturities);
	~CurveHolder();
	YieldCurve *get() const noexcept final { return pricing_curve_; }

	void dump(const char *heading)
	{
		printf("%s\n", heading);
		for (int i = 0; i < n_maturities_; i++) {
			printf("Curve[%d] = {%d, %.12f}\n", i, maturities_[i], values_[i]);
		}
	}

	void init(Date asOfDate, Allocator *alloc, CurveId curveId, InterpolatorType interp, IRRateType interpolated_on)
	{
		assert(!base_curve_);
		as_of_date_ = asOfDate;
		interpolated_on_ = interpolated_on;
		std::copy(rates_, rates_ + n_maturities_, rates_plus_h_);
		std::copy(rates_, rates_ + n_maturities_, rates_minus_h_);
		if (curve_type_ == CurveType::CURVE_TYPE_INTERPOLATED &&
		    interpolated_on == IRRateType::DISCOUNT_FACTOR) {
			for (int i = 0; i < n_maturities_; i++) {
				// Convert to discount factor
				double t = fraction_->year_fraction(asOfDate, maturities_[i]);
				values_[i] = values_plus_h_[i] = values_minus_h_[i] = std::exp(-rates_[i] * t);
			}
		} else {
			std::copy(rates_, rates_ + n_maturities_, values_);
			std::copy(rates_, rates_ + n_maturities_, values_plus_h_);
			std::copy(rates_, rates_ + n_maturities_, values_minus_h_);
		}
		switch (curve_type_) {
		default:
		case CurveType::CURVE_TYPE_INTERPOLATED:
			base_curve_ = make_curve(alloc, curveId, asOfDate, maturities_, values_, n_maturities_, interp,
						 interpolated_on_);
			curve_plus_h_ = make_curve(alloc, curveId, asOfDate, maturities_, values_plus_h_, n_maturities_,
						   interp, interpolated_on_);
			curve_minus_h_ = make_curve(alloc, curveId, asOfDate, maturities_, values_minus_h_,
						    n_maturities_, interp, interpolated_on_);
			break;
		case CurveType::CURVE_TYPE_SVENSSON_PARAMETRIC: {
			// A this point the have the same initial guess parameters
			// TODO check that input has 6 parameters
			assert(6 == n_maturities_);
			base_curve_ = make_svensson_curve(alloc, curveId, asOfDate, values_, n_maturities_);
			curve_plus_h_ = make_svensson_curve(alloc, curveId, asOfDate, values_, n_maturities_);
			curve_minus_h_ = make_svensson_curve(alloc, curveId, asOfDate, values_, n_maturities_);
			break;
		}
		}
		pricing_curve_ = base_curve_.get();
	}

	// Bumps a particular rate (identified by pillar) up and down
	// We always bump rates, never discount factors, when interpolating
	// in discount factor space we derive the discount factors from the
	// bumped rates
	void bump(uint32_t pillar, const double h)
	{
		std::copy(rates_, rates_ + n_maturities_, rates_plus_h_);
		rates_plus_h_[pillar] = rates_plus_h_[pillar] + h;
		std::copy(rates_, rates_ + n_maturities_, rates_minus_h_);
		rates_minus_h_[pillar] = rates_minus_h_[pillar] - h;
		for (int i = 0; i < n_maturities_; i++) {
			if (curve_type_ == CurveType::CURVE_TYPE_INTERPOLATED &&
			    interpolated_on_ == IRRateType::DISCOUNT_FACTOR) {
				// Convert to discount factor
				double t = fraction_->year_fraction(as_of_date_, maturities_[i]);
				values_plus_h_[i] = std::exp(-rates_plus_h_[i] * t);
				values_minus_h_[i] = std::exp(-rates_minus_h_[i] * t);
			} else {
				values_plus_h_[i] = rates_plus_h_[i];
				values_minus_h_[i] = rates_minus_h_[i];
			}
		}
		curve_plus_h_->update_rates(values_plus_h_, n_maturities_);
		curve_minus_h_->update_rates(values_minus_h_, n_maturities_);
	}

	void save()
	{
		std::copy(rates_, rates_ + n_maturities_, saved_rates_);
		std::copy(values_, values_ + n_maturities_, saved_values_);
	}
	void restore()
	{
		std::copy(saved_rates_, saved_rates_ + n_maturities_, rates_);
		std::copy(saved_values_, saved_values_ + n_maturities_, values_);
	}

	void set_vector(const double *x, unsigned n, bool is_correction = false)
	{
		assert(n == n_maturities_);
		for (int pillar = 0; pillar < n_maturities_; pillar++) {
			if (is_correction)
				rates_[pillar] = rates_[pillar] - x[pillar];
			else
				rates_[pillar] = x[pillar];
			if (curve_type_ == CurveType::CURVE_TYPE_INTERPOLATED &&
			    interpolated_on_ == IRRateType::DISCOUNT_FACTOR) {
				values_[pillar] = std::exp(-rates_[pillar] *
							   fraction_->year_fraction(as_of_date_, maturities_[pillar]));
			} else {
				values_[pillar] = rates_[pillar];
			}
		}
		if (is_trace_enabled()) {
			dump("Updated curve values");
		}
		base_curve_->update_rates(values_, n_maturities_);
	}

	void get_vector(double *x, unsigned n)
	{
		assert(n == n_maturities_);
		for (uint32_t pillar = 0; pillar < n_maturities_; pillar++) {
			x[pillar] = rates_[pillar];
		}
	}
};

CurveHolder::CurveHolder(CurveType curve_type, uint32_t n)
    : curve_type_(curve_type), n_maturities_(n), as_of_date_(0), pricing_curve_(nullptr)
{
	// FIXME use allocator
	fraction_ = get_day_fraction(DayCountFraction::ACT_365_FIXED);
	maturities_ = new Date[n];
	rates_ = new double[n];
	rates_plus_h_ = new double[n];
	rates_minus_h_ = new double[n];
	saved_rates_ = new double[n];
	values_ = new double[n];
	values_plus_h_ = new double[n];
	values_minus_h_ = new double[n];
	saved_values_ = new double[n];
}

CurveHolder::~CurveHolder()
{
	// FIXME use allocator
	delete[] maturities_;
	delete[] rates_;
	delete[] rates_plus_h_;
	delete[] rates_minus_h_;
	delete[] saved_rates_;
	delete[] values_;
	delete[] values_plus_h_;
	delete[] values_minus_h_;
	delete[] saved_values_;
}

// A wrapper for each instrument
// Holds the cashflow structure obtained from Lua script
// as well as the processed one
struct CashflowInstrument {
	private:
	std::unique_ptr<CFCollection> cfcollection_;
	Cashflows *cashflows_;
	Date maturity_;
	// As the instrument is arranged in sorted order by maturity
	// its position may be different from the input data
	// So we track the original position where this instrument was located
	// This can be used to get the instrument attributes from the input
	// curve
	int original_instrument_position_;
	// A constraint does not contribute a pillar
	bool is_constraint_;

	BootstrapCurveMapper mapper_;
	Date business_date_;

	std::unique_ptr<CFCollection> cfcollection_backup_;
	Cashflows *cashflows_backup_;

	public:
	CashflowInstrument() : maturity_(0), original_instrument_position_(-1), is_constraint_(false), business_date_(0)
	{
	}
	~CashflowInstrument() {}
	void backup_cashflows()
	{
		assert(!cfcollection_backup_);
		cfcollection_backup_ = std::move(cfcollection_);
		cashflows_backup_ = cashflows_;
	}
	void restore_cashflows()
	{
		assert(cfcollection_backup_);
		cfcollection_ = std::move(cfcollection_backup_);
		cashflows_ = cashflows_backup_;
	}
	void copy_cashflows(std::unique_ptr<CashflowInstrument> other)
	{
		cfcollection_ = std::move(other->cfcollection_);
	}
	void set_cfcollection(std::unique_ptr<CFCollection> collection)
	{
		assert(collection);
		cfcollection_ = std::move(collection);
	}
	void set_maturity(Date dt)
	{
		assert(is_valid_date(dt));
		maturity_ = dt;
	}
	void set_original_instrument_position(int pos) { original_instrument_position_ = pos; }
	void set_is_constraint(bool value) { is_constraint_ = value; }
	Date maturity() const { return maturity_; }
	void set_business_date(Date business_date) { business_date_ = business_date; }
	void set_curve_mapper(const BootstrapCurveMapper &mapper) { mapper_ = mapper; }
	int original_instrument_position() const { return original_instrument_position_; }
	bool build_cashflows(RegionAllocator *A)
	{
		ValuationContextImpl ctx(business_date_, 0);
		cashflows_ = construct_cashflows(A, cfcollection_.get(), ctx, &mapper_);
		return !!cashflows_;
	}
	double evaluate(FixedRegionAllocator *A, const ValuationContext &ctx, CurveProvider *curve_provider,
			Sensitivities &sens, StatusCode &status)
	{
		return compute_present_value(A, ctx, cashflows_, curve_provider, sens, status);
	}
	double evaluate(FixedRegionAllocator *A, const ValuationContext &ctx, CurveReference *forward_curve,
			CurveReference *discount_curve, Sensitivities &sens, StatusCode &status)
	{
		return compute_present_value(A, cashflows_, discount_curve, forward_curve, forward_curve, ctx, sens,
					     status);
	}
	int forward_curve() { return mapper_.forward_curve(); }
	int discount_curve() { return mapper_.discount_curve(); }

	void get_name(char *buffer, size_t n)
	{
		snprintf(buffer, n, "CashflowInstrument{pos=%d,maturity=%d,constraint=%d}",
			 original_instrument_position_, maturity_, (int)is_constraint_);
	}
};

struct CashflowInstrumentByCurve {
	std::map<Date, std::unique_ptr<CashflowInstrument>> instruments_by_maturity;
};

class LuaWrapper
{
	private:
	lua_State *L_;

	public:
	LuaWrapper()
	{
		L_ = luaL_newstate(); /* create Lua state */
		if (L_ == nullptr) {
			die("Failed to create Lua State\n");
		}
		luaL_checkversion(L_); /* check that interpreter has correct version */
		luaL_openlibs(L_);     /* open standard libraries */
	}
	~LuaWrapper()
	{
		if (L_ != nullptr)
			lua_close(L_);
	}
	bool load_luascript(std::string lua_script)
	{
		if (luaL_loadfile(L_, lua_script.c_str()) != LUA_OK) {
			error("Failed to load Lua script %s: %s\n", lua_script.c_str(), lua_tostring(L_, -1));
			lua_pop(L_, 1); // pop the error message
			return false;
		}
		if (lua_pcall(L_, 0, 0, 0) != LUA_OK) {
			error("Failed to run Lua script %s: %s\n", lua_script.c_str(), lua_tostring(L_, -1));
			lua_pop(L_, 1); // pop the error message
			return false;
		}
		return true;
	}
	lua_State *L() const { return L_; }
};

typedef double (*PenaltyFunction)(const double *x, const double *lower_bound, const double *upper_bound, unsigned n);
double default_penalty_function(const double *x, const double *lower_bound, const double *upper_bound, unsigned n)
{
	double penalty = 0.0;
	for (int i = 0; i < n; i++) {
		if ((lower_bound[i] != -HUGE_VAL && x[i] < lower_bound[i]) ||
		    (upper_bound[i] != HUGE_VAL && x[i] > upper_bound[i])) {
			double y = x[i];
			if (fabs(y) < 1.0)
				y = fabs(y) + 1.0;
			penalty += (y * y);
		}
	}
	return penalty;
}

class CurveBuilderOptions
{
	public:
	SolverType solver_type;
	int max_iterations;
	bool generate_par_sensitivities;
	PenaltyFunction penalty_function;

	public:
	CurveBuilderOptions()
	    : solver_type(SolverType::SOLVER_TYPE_LEVENBERG_MARQUARDT), max_iterations(10),
	      generate_par_sensitivities(false), penalty_function(default_penalty_function)
	{
	}
};

class ParSensitivities
{
	private:
	int m_;			     // rows
	int n_;			     // columns
	std::vector<double> values_; // row major data

	public:
	ParSensitivities(int m, int n) : m_(m), n_(n), values_(m * n) {}
	void set(int row, int col, double v)
	{
		int pos = row * n_ + col;
		values_[pos] = v;
	}
	void dump(FILE *fp)
	{
		if (!is_trace_enabled())
			return;
		fprintf(fp, "Par Sensitivities dump\n");
		for (int row = 0; row < m_; row++) {
			fprintf(fp, "row[%d]", row);
			for (int col = 0; col < n_; col++) {
				double v = values_[row * n_ + col];
				fprintf(fp, "\t%f", v);
			}
			fprintf(fp, "\n");
		}
	}
	int m() const { return m_; }
	int n() const { return n_; }
	double value(int row, int col) const
	{
		assert(row < m_);
		assert(col < n_);
		int pos = row * n_ + col;
		return values_[pos];
	}
};

class SolverFunction;
class CurveBuilder
{
	private:
	const CurveDefinitionProvider *definition_provider_;
	// The input curves containing raw instruments data
	const ParCurveSet *input_curves_;
	// The processed instruments stored by each input curve; curves are in
	// same order as input_curves_ but instruments are sorted by maturity
	// date
	std::vector<std::unique_ptr<CashflowInstrumentByCurve>> instruments_by_curve_;
	// The global list of instruments, sorted by curve index and maturity
	// date as the composite key (uint64)
	std::map<uint64_t, CashflowInstrument *> global_instruments_list_;
	// The work in progress curves where the processing results are stored
	std::vector<std::unique_ptr<CurveHolder>> bootstrap_curves_;
	// A lookup map to allow mapping of curve definition id to its index,
	// i.e. position in the vectors above
	std::map<int, int> curve_by_name_;
	std::vector<std::unique_ptr<ParSensitivities>> par_sensitivities_by_curve_;
	// The as of date
	Date business_date_;
	LuaWrapper &lua_;
	StatusCode last_error_code_;
	char last_error_message_[1024];
	DynamicRegionAllocator cashflow_allocator_;

	public:
	CurveBuilder(LuaWrapper &lua, Date business_date, const CurveDefinitionProvider *definition_provider);
	~CurveBuilder();
	bool add_input_curves(const ParCurveSet *curves);
	const ParCurve &get_input_curve(int i) const
	{
		assert(i >= 0 && i < input_curves_->par_curves_size());
		return input_curves_->par_curves(i);
	}
	const char *last_error_message() const { return last_error_message_; }
	StatusCode last_error_code() const { return last_error_code_; }
	// Generate an instrument definition
	std::unique_ptr<CashflowInstrument> build_instrument(const ParCurve &curve, int instrument_position,
							     double par_rate);
	// Adds instruments for the curve into the global list
	// The instruments are added such that the global instrument list
	// is sorted by curve_index and instrument maturity date
	// The 'curve_index' here means the position of the par curve
	// in the input par curve set.
	// curve_index is also the position of the curve in the input_curves_
	// and instruments_by_curve_ vectors
	void add_to_global_instrument_list(int curve_index, const CashflowInstrumentByCurve *curve);
	// Construct the initial curve based on the input data
	// If maturities are derived from instruments then this curve will
	// have as many pillars as there are input instruments, else
	// the maturities will be based on fixed set of tenors
	bool create_bootstrap_curve(const ParCurve &input_curve, const CashflowInstrumentByCurve *instruments_by_curve);
	bool prepare_instruments();
	int num_instruments() const { return global_instruments_list_.size(); }
	int num_pillars() const
	{
		int n = 0;
		for (auto &&curve : bootstrap_curves_) {
			n += curve->n_maturities_;
		}
		return n;
	}
	int num_curves() const { return bootstrap_curves_.size(); }
	std::map<uint64_t, CashflowInstrument *> &sorted_instruments() { return global_instruments_list_; }
	CurveHolder *get_curve_by_index(int id) const
	{
		assert(id >= 0 && id < bootstrap_curves_.size());
		return bootstrap_curves_[id].get();
	}
	CurveHolder *get_curve_by_definition_id(int id) const
	{
		auto iter = curve_by_name_.find(id);
		if (iter == std::end(curve_by_name_))
			return nullptr;
		return get_curve_by_index(iter->second);
	}
	ParSensitivities *get_sensitivities_by_index(int id) const
	{
		assert(id >= 0);
		if (id >= 0 && id < par_sensitivities_by_curve_.size())
			return par_sensitivities_by_curve_[id].get();
		else
			return nullptr;
	}
	ParSensitivities *get_sensitivities_by_curve_definition_id(int id) const
	{
		auto iter = curve_by_name_.find(id);
		if (iter == std::end(curve_by_name_))
			return nullptr;
		return get_sensitivities_by_index(iter->second);
	}
	Date business_date() const { return business_date_; }
	// The list of curves being bootstrapped
	// In the same order as the input curves
	std::vector<std::unique_ptr<CurveHolder>> &curves() { return bootstrap_curves_; }
	StatusCode build_curves(const CurveBuilderOptions &options);
	std::unique_ptr<SolverFunction> make_solver(const CurveBuilderOptions &options);
};

CurveBuilder::CurveBuilder(LuaWrapper &lua, Date business_date, const CurveDefinitionProvider *definition_provider)
    : definition_provider_(definition_provider), business_date_(business_date), lua_(lua), last_error_code_(kOk),
      cashflow_allocator_(64 * 1024)
{
	last_error_message_[0] = 0;
}

CurveBuilder::~CurveBuilder() {}

/*
** Message handler used to run all chunks
*/
static int msghandler(lua_State *L)
{
	const char *msg = lua_tostring(L, 1);
	if (msg == NULL) {				 /* is error object not a string? */
		if (luaL_callmeta(L, 1, "__tostring") && /* does it have a metamethod */
		    lua_type(L, -1) == LUA_TSTRING)      /* that produces a string? */
			return 1;			 /* that is the message */
		else
			msg = lua_pushfstring(L, "(error object is a %s value)", luaL_typename(L, 1));
	}
	luaL_traceback(L, L, msg, 1); /* append a standard traceback */
	return 1;		      /* return the traceback */
}

/*
** Interface to 'lua_pcall', which sets appropriate message function
** and C-signal handler. Used to run all chunks.
*/
static int docall(lua_State *L, int narg, int nres)
{
	int status;
	int base = lua_gettop(L) - narg;  /* function index */
	lua_pushcfunction(L, msghandler); /* push message handler */
	lua_insert(L, base);		  /* put it under function and args */
	// globalL = L;  /* to be available to 'laction' */
	// signal(SIGINT, laction);  /* set C-signal handler */
	status = lua_pcall(L, narg, nres, base);
	// signal(SIGINT, SIG_DFL); /* reset C-signal handler */
	lua_remove(L, base); /* remove message handler from the stack */
	return status;
}

// Generates CashflowInstrument for the specified
// instrument in the parcurve; Lua script is called to obtain
// cashflow definition
std::unique_ptr<CashflowInstrument> CurveBuilder::build_instrument(const ParCurve &par_curve, int instrument_position,
								   double par_rate)
{
	std::unique_ptr<CashflowInstrument> cfinst;
	const IRCurveDefinition *curve_defn =
	    definition_provider_->get_definition_by_id(par_curve.curve_definition_id());
	if (curve_defn == nullptr) {
		last_error_code_ = StatusCode::kBTS_CurveDefinitionNotFound;
		snprintf(last_error_message_, sizeof last_error_message_, "%s: ID %d", error_message(last_error_code_),
			 par_curve.curve_definition_id());
		error("%s\n", last_error_message_);
		return cfinst;
	}
	auto const &inst = par_curve.instruments(instrument_position);
	// auto par_rate = par_curve.par_rates().values(instrument_position);
	auto lua_function_name = inst.instrument_type();
	lua_State *L_ = lua_.L();
	int top = lua_gettop(L_);
	if (lua_getglobal(L_, lua_function_name.c_str()) == LUA_TFUNCTION) {
		lua_pushinteger(L_, business_date_);
		lua_pushstring(L_, get_default_converter()->currency_to_string(curve_defn->currency()));
		lua_pushstring(L_, get_default_converter()->index_family_to_string(curve_defn->index_family()));
		lua_pushstring(L_, get_default_converter()->tenor_to_string(curve_defn->tenor()).c_str());
		lua_pushstring(L_, inst.instrument_key().c_str());
		lua_pushnumber(L_, par_rate);
		int nargs = 6;
		if (inst.floating_tenor() != Tenor::TENOR_UNSPECIFIED) {
			lua_pushstring(L_, get_default_converter()->tenor_to_string(inst.floating_tenor()).c_str());
			nargs++;
		}
		if (docall(L_, nargs, 2) == LUA_OK) {
			if (lua_type(L_, -2) == LUA_TNUMBER && lua_type(L_, -1) == LUA_TUSERDATA) {
				Date maturity_date = (Date)lua_tointeger(L_, -2);
				auto cfcollection = raviapi_get_CFCollection(L_, -1);
				if (is_valid_date(maturity_date) && cfcollection) {
					cfinst = std::unique_ptr<CashflowInstrument>(new CashflowInstrument());
					cfinst->set_cfcollection(std::move(cfcollection));
					cfinst->set_maturity(maturity_date);
					cfinst->set_original_instrument_position(instrument_position);
					// cfinst->set_is_constraint(inst.is_constraint());
					cfinst->set_is_constraint(false);
					cfinst->set_business_date(business_date_);
				} else {
					last_error_code_ = StatusCode::kBTS_CashflowDefinitionFailed;
					snprintf(last_error_message_, sizeof last_error_message_, "%s: \n%s",
						 error_message(last_error_code_), inst.DebugString().c_str());
					error("%s\n", last_error_message_);
				}
			} else {
				last_error_code_ = StatusCode::kBTS_CashflowDefinitionFailed;
				snprintf(last_error_message_, sizeof last_error_message_, "%s: \n%s",
					 error_message(last_error_code_), inst.DebugString().c_str());
				error("%s\n", last_error_message_);
			}
			// Pop the return values
			lua_pop(L_, 2);
		} else {
			last_error_code_ = StatusCode::kBTS_LuaScriptFailure;
			snprintf(last_error_message_, sizeof last_error_message_, "%s: %s \n%s",
				 error_message(last_error_code_), lua_tostring(L_, -1), inst.DebugString().c_str());
			error("%s\n", last_error_message_);
			// Pop the error object
			lua_pop(L_, 1);
		}
	} else {
		error("Could not find a Lua pricing routine for %s\n", lua_function_name.c_str());
		/* pop the value pushed by lua_getglobal() */
		last_error_code_ = StatusCode::kBTS_LuaPricingRoutineNotFound;
		snprintf(last_error_message_, sizeof last_error_message_, "%s: %s", error_message(last_error_code_),
			 lua_function_name.c_str());
		error("%s\n", last_error_message_);
		lua_pop(L_, 1);
	}
	assert(top == lua_gettop(L_));
	return cfinst;
}

// Adds instruments for the curve into the global list
// The instruments are added such that the global instrument list
// is sorted by curve_index and instrument maturity date
// The 'curve_index' here means the position of the par curve
// in the input par curve set.
void CurveBuilder::add_to_global_instrument_list(int curve_index, const CashflowInstrumentByCurve *curve)
{
	for (auto &&iter = std::begin(curve->instruments_by_maturity); iter != std::end(curve->instruments_by_maturity);
	     iter++) {
		uint64_t id = ((uint64_t)curve_index) << 32 | (uint64_t)(iter->first);
		global_instruments_list_.insert({id, iter->second.get()});
	}
}

// Construct the initial curve based on the input data
// If maturities are derived from instruments then this curve will
// have as many pillars as there are input instruments, else
// the maturities will be based on fixed set of tenors
// The curve will be added to the list of curves
bool CurveBuilder::create_bootstrap_curve(const ParCurve &input_curve,
					  const CashflowInstrumentByCurve *instruments_by_curve)
{
	const IRCurveDefinition *defn = definition_provider_->get_definition_by_id(input_curve.curve_definition_id());
	if (defn->curve_type() == CurveType::CURVE_TYPE_SVENSSON_PARAMETRIC) {
		if (defn->maturity_generation_rule() != MaturityGenerationRule::MATURITY_GENERATION_RULE_SVENSSON ||
		    defn->interpolated_on() != IRRateType::ZERO_RATE) {
			last_error_code_ = StatusCode::kBTS_BadInterpolatedOn;
			snprintf(last_error_message_, sizeof last_error_message_,
				 "%s: Svensson curves must be on zero rates, and maturity "
				 "generation rule must be set to MATURITY_GENERATION_RULE_SVENSSON",
				 error_message(last_error_code_));
			error("%s\n", last_error_message_);
			return false;
		}
	}
	assert(defn != nullptr);
	std::unique_ptr<CurveHolder> ycurve;
	switch (defn->maturity_generation_rule()) {
	case MaturityGenerationRule::MATURITY_GENERATION_RULE_FIXED_TENORS: {
		// input curve has defined set of maturities
		// FIXME we need to implement some logic to base the initial
		// guess on par rates as otherwise the solver may not
		// converge. For now for discount factors we set guess to 1.0,
		// for rates to 5%.
		ycurve = std::make_unique<CurveHolder>(defn->curve_type(), defn->tenors_size());
		auto converter = get_default_converter();
		for (int i = 0; i < defn->tenors_size(); i++) {
			int offset = converter->tenor_to_days(defn->tenors(i));
			Date maturity = business_date_ + offset;
			ycurve->maturities_[i] = maturity;
			ycurve->rates_[i] = 0.05;
		}
		break;
	}
	case MaturityGenerationRule::MATURITY_GENERATION_RULE_SVENSSON: {
		ycurve = std::make_unique<CurveHolder>(defn->curve_type(),
						       6); // We have six parameters
		for (int i = 0; i < 6; i++) {
			ycurve->maturities_[i] = 0; // Dummy value as we don't have maturities
			ycurve->rates_[i] = 0.5;    // Initial guess
		}
		break;
	}
	default:
	case MaturityGenerationRule::MATURITY_GENERATION_RULE_DERIVE_FROM_INSTRUMENTS: {
		// We will use the input instruments to define the set of
		// maturities
		uint32_t n = input_curve.instruments_size();
		ycurve = std::make_unique<CurveHolder>(defn->curve_type(), n);
		// Iterate over the instruments by maturity date
		int i = 0;
		for (auto &&iter = std::begin(instruments_by_curve->instruments_by_maturity);
		     iter != std::end(instruments_by_curve->instruments_by_maturity); iter++) {
			auto par_rate = input_curve.par_rates().values(iter->second->original_instrument_position());
			if (par_rate < -0.1 || par_rate > 0.3) {
				// Hmmm bad input
				last_error_code_ = StatusCode::kBTS_BadInputCurves;
				snprintf(last_error_message_, sizeof last_error_message_,
					 "%s: index %d par rate %.12f, expected range (-0.1,0.3)",
					 error_message(last_error_code_), i, par_rate);
				error("%s\n", last_error_message_);
				return false;
			}
			ycurve->maturities_[i] = iter->first;
			ycurve->rates_[i] = par_rate;
			if (ycurve->maturities_[i] <= business_date_ ||
			    (i > 0 && ycurve->maturities_[i - 1] >= ycurve->maturities_[i])) {
				last_error_code_ = StatusCode::kBTS_BadMaturitiesSortOrder;
				snprintf(last_error_message_, sizeof last_error_message_, "%s: index %d",
					 error_message(last_error_code_), i);
				error("%s\n", last_error_message_);
				return false;
			}
			i++;
		}
		break;
	}
	}
	ycurve->init(business_date_, get_default_allocator(),
		     make_curve_id(PricingCurveType::PRICING_CURVE_TYPE_UNSPECIFIED, defn->currency(),
				   defn->index_family(), defn->tenor(), business_date_),
		     defn->interpolator_type(), defn->interpolated_on());
	if (is_trace_enabled()) {
		ycurve->dump("Initial Curve Values");
	}
	this->bootstrap_curves_.push_back(std::move(ycurve));
	return true;
}

bool CurveBuilder::add_input_curves(const ParCurveSet *curve)
{
	input_curves_ = curve;
	bool success = true;
	for (int k = 0; k < input_curves_->par_curves_size(); k++) {
		// TODO validate input curve
		// TODO If specs supplied then do the maturities make sense /
		// are sorted / > business date etc.
		auto &underlying = input_curves_->par_curves(k);
		if (curve_by_name_.find(underlying.curve_definition_id()) != curve_by_name_.end()) {
			last_error_code_ = StatusCode::kBTS_DuplicateCurve;
			snprintf(last_error_message_, sizeof last_error_message_, "%s: ID %d",
				 error_message(last_error_code_), underlying.curve_definition_id());
			error("%s\n", last_error_message_);
			return false;
		}
		std::unique_ptr<CashflowInstrumentByCurve> inst_by_curve =
		    std::make_unique<CashflowInstrumentByCurve>();
		// For each instrument in the curve call a Lua script which
		// returns the instruments maturity date and the cashflow
		// representation
		for (int i = 0; i < underlying.instruments_size(); i++) {
			auto par_rate = underlying.par_rates().values(i);
			std::unique_ptr<CashflowInstrument> cfinst = build_instrument(underlying, i, par_rate);
			if (cfinst) {
				inst_by_curve->instruments_by_maturity.insert(
				    std::make_pair(cfinst->maturity(), std::move(cfinst)));
			} else {
				success = false;
				break;
			}
		}
		if (success) {
			auto curve_ptr = inst_by_curve.get();
			instruments_by_curve_.push_back(std::move(inst_by_curve));
			// Add instruments in this curve to the global
			// insruments list
			// k is the curve index
			add_to_global_instrument_list(k, curve_ptr);
			if (!create_bootstrap_curve(underlying, curve_ptr)) {
				success = false;
				break;
			}
			curve_by_name_.insert({underlying.curve_definition_id(), k});
		} else {
			break;
		}
	}
	return success;
}

bool CurveBuilder::prepare_instruments()
{
	// Resolve references to forward / discount curves in instruments
	ValuationContextImpl ctx(business_date_, 0);
	for (auto &&iter : global_instruments_list_) {
		int curve_index = (int)(iter.first >> 32);
		assert(curve_index >= 0 && curve_index < input_curves_->par_curves_size());
		auto &curve = input_curves_->par_curves(curve_index);
		CashflowInstrument *cfinst = iter.second;
		assert(cfinst->original_instrument_position() < curve.instruments_size());
		auto &input_instrument = curve.instruments(cfinst->original_instrument_position());
		// The forward and discount curves need to be mapped to the ones
		// specified in the instrument defintion - the curves are
		// identified by name in the input so we need to convert from
		// name to index, and create a mapper that forces the mapping
		// output to resolve to the predefined curves
		int discount_curve = -1;
		int forward_curve = -1;
		auto iter2 = curve_by_name_.find(input_instrument.forward_curve_definition_id());
		if (iter2 != curve_by_name_.end())
			forward_curve = iter2->second;
		iter2 = curve_by_name_.find(input_instrument.discount_curve_definition_id());
		if (iter2 != curve_by_name_.end())
			discount_curve = iter2->second;
		// It is an error if the curves were not found
		if (forward_curve < 0) {
			last_error_code_ = kBTS_ForwardCurveReferenceNotFound;
			snprintf(last_error_message_, sizeof last_error_message_, "%s: ID %d",
				 error_message(last_error_code_), input_instrument.forward_curve_definition_id());
			error("%s\n", last_error_message_);
			return false;
		}
		if (forward_curve < 0 || discount_curve < 0) {
			last_error_code_ = StatusCode::kBTS_DiscountCurveReferenceNotFound;
			snprintf(last_error_message_, sizeof last_error_message_, "%s: ID %d",
				 error_message(last_error_code_), input_instrument.discount_curve_definition_id());
			error("%s\n", last_error_message_);
			return false;
		}
		PricingCurve fc(
		    PricingCurveType::PRICING_CURVE_TYPE_FORWARD,
		    definition_provider_
			->get_definition_by_id(input_curves_->par_curves(forward_curve).curve_definition_id())
			->currency());
		PricingCurve dc(
		    PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT,
		    definition_provider_
			->get_definition_by_id(input_curves_->par_curves(discount_curve).curve_definition_id())
			->currency());
		BootstrapCurveMapper mapper(fc, forward_curve, dc, discount_curve);
		cfinst->set_curve_mapper(mapper);
		if (!cfinst->build_cashflows(&cashflow_allocator_)) {
			last_error_code_ = kBTS_CashflowGenerationFailed;
			snprintf(last_error_message_, sizeof last_error_message_, "%s: %s",
				 error_message(last_error_code_), input_instrument.DebugString().c_str());
			error("%s\n", last_error_message_);
			return false;
		}
	}
	return true;
}

int minpack_lmder_function(void *p, int m, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag);

class SolverFunction
{
	protected:
	CurveBuilder *curve_builder_;
	std::map<std::string, std::shared_ptr<YieldCurve>> output_curves_;
	int max_iterations_;
	Allocator *alloc_;
	FILE *debug_output_;
	double *A_;  // sized m x n
	double *bx_; // sized m
	double *upper_bounds_;
	double *lower_bounds_;
	Sensitivities dummy_sens_;
	PenaltyFunction penalty_function_;
	FixedRegionAllocator *pricing_allocator_;

	public:
	SolverFunction(CurveBuilder *builder, int max_iterations, Allocator *alloc, FILE *debug_output)
	    : curve_builder_(builder), max_iterations_(max_iterations), alloc_(alloc), debug_output_(debug_output),
	      dummy_sens_(alloc), penalty_function_(default_penalty_function)
	{
		A_ = new double[curve_builder_->num_instruments() * curve_builder_->num_pillars()];
		// Note that bx is sized to be as big as the m() = num
		// instruments although we are only interested in n (pillars)
		// because the LAPACK routines require this. But when applying
		// corrections and checking the solution we must only look at
		// n() elements in bx
		bx_ = new double[curve_builder_->num_instruments()];
		upper_bounds_ = new double[curve_builder_->num_pillars()];
		std::fill_n(upper_bounds_, curve_builder_->num_pillars(), HUGE_VAL);
		lower_bounds_ = new double[curve_builder_->num_pillars()];
		std::fill_n(lower_bounds_, curve_builder_->num_pillars(), -HUGE_VAL);
		pricing_allocator_ = get_threadspecific_allocators()->tempspace_allocator;
	}

	virtual ~SolverFunction()
	{
		delete[] bx_;
		delete[] upper_bounds_;
		delete[] lower_bounds_;
		delete[] A_;
	}
	// number of functions
	// each instrument is a row in the Jacobian
	int32_t m() { return (int32_t)curve_builder_->num_instruments(); }

	// number of variables
	// consists of all pillars of all curves arranged
	// by curve
	// each pillar is a column in the Jacobian
	int32_t n() { return (int32_t)curve_builder_->num_pillars(); }

	double *solution_vector() { return bx_; }

	double *jacobian() { return A_; }

	// Reset any cached data
	void reset() {}

	void set_penalty_function(PenaltyFunction f) { this->penalty_function_ = f; }

	double evaluate_instrument(CashflowInstrument *cfinst, StatusCode &status, double *x, unsigned n)
	{
		CurveReference *forward_curve = curve_builder_->get_curve_by_index(cfinst->forward_curve());
		CurveReference *discount_curve = curve_builder_->get_curve_by_index(cfinst->discount_curve());
		ValuationContextImpl ctx(curve_builder_->business_date());
		BootstrapCurveProvider provider(discount_curve, forward_curve);
		return cfinst->evaluate(pricing_allocator_, ctx, &provider, dummy_sens_, status) +
		       (penalty_function_ != nullptr ? penalty_function_(x, lower_bounds_, upper_bounds_, n) : 0.0);
	}

	void dump_solution_vector(const char *description, bool force = false)
	{
		if (!is_trace_enabled() && !force)
			return;
		unsigned int instrument_num = 0;
		char instrument_name[128];
		printf("%s\n", description);
		for (auto &item : curve_builder_->sorted_instruments()) {
			item.second->get_name(instrument_name, sizeof instrument_name);
			printf("fx[%d] = %.12f, %s\n", instrument_num, bx_[instrument_num], instrument_name);
			instrument_num++;
		}
		printf("Curve values\n");
		for (auto &ch : curve_builder_->curves()) {
			ch->dump(ch->base_curve_->name().c_str());
		}
	}

	// Evaluate the system of equations and populate the vector bx
	bool evaluate()
	{
		reset();
		StatusCode status = StatusCode::kOk;
		double *x = (double *)alloca(sizeof(double) * n());
		get_vector(x, n());
		// and populate vector
		unsigned int instrument_num = 0;
		for (auto &item : curve_builder_->sorted_instruments()) {
			bx_[instrument_num] = evaluate_instrument(item.second, status, x, n());
			if (status != StatusCode::kOk) {
				char instrument_name[128];
				item.second->get_name(instrument_name, sizeof instrument_name);
				error("Failed to evaluate instrument %s\n", instrument_name);
				return false;
			}
			instrument_num++;
		}
		dump_solution_vector("Computed Solution Vector");
		return true;
	}

	// calculate the Jacobian numerically
	bool compute_jacobian()
	{
		// Get the current solution vector
		// we can use this as input to penalty function
		// We ignore the bumped values here
		double *x = (double *)alloca(sizeof(double) * n());
		get_vector(x, n());

		double epsmch = __cminpack_func__(dpmpar)(1);
		double eps = sqrt(epsmch);

		redukti_matrix_t X = {m(), n(), A_};
		// To populate the Jacobian , we process by curve pillar by
		// instrument so that we can setup the curve once for each
		// pillar rate_num tracks where we are in terms of the
		// curve/pillar combination - i.e. rate_num is 0 for first
		// curve/pillar

		int rate_num = 0;
		int curve_num = 0;
		StatusCode status = StatusCode::kOk;
		for (auto &ch : curve_builder_->curves()) {
			// for each pillar on the curve we
			// compute the partial derivative for each instrument's
			// PV with respect to the pillar - we use a numerical
			// approach to ensure we can handle unknown
			// interpolation methods
			for (uint32_t pillar = 0; pillar < ch->n_maturities_; pillar++) {
				const double h =
				    ch->curve_type_ == CurveType::CURVE_TYPE_SVENSSON_PARAMETRIC
					? (ch->rates_[pillar] == 0.0 ? eps : eps * fabs(ch->rates_[pillar]))
					: 0.0001; // use 1 basis point shifts for regular curves.
				ch->bump(pillar, h);

				bool forward_diff = false;
				// FIXME this check should use the lower_bounds array
				if (ch->curve_type_ == CurveType::CURVE_TYPE_SVENSSON_PARAMETRIC &&
				    (pillar == 0 || pillar == 4 || pillar == 5)) {
					forward_diff = true;
				}

				size_t instrument_num = 0;
				for (auto &item : curve_builder_->sorted_instruments()) {
					double partial = 0.0;

					if (!forward_diff) {
						// Central diff
						// First calculate PV using -h
						// bump
						ch->pricing_curve_ = ch->curve_minus_h_.get();
						double minus_h_pv = evaluate_instrument(item.second, status, x, n());
						if (status != StatusCode::kOk)
							return false;

						// Now calculate PV using +h
						// bump
						ch->pricing_curve_ = ch->curve_plus_h_.get();
						double plus_h_pv = evaluate_instrument(item.second, status, x, n());
						if (status != StatusCode::kOk)
							return false;

						// Calculate partial derivative
						if (plus_h_pv != minus_h_pv) {
							partial = (plus_h_pv - minus_h_pv) / (2 * h);
						}
					} else {
						// FIXME as the base curve doesn't change we don't need
						// to evaluate every time
						// Can we just use bx here?
						ch->pricing_curve_ = ch->base_curve_.get();
						double base_pv = evaluate_instrument(item.second, status, x, n());
						if (status != StatusCode::kOk)
							return false;

						// Now calculate PV using +h
						// bump
						ch->pricing_curve_ = ch->curve_plus_h_.get();
						double plus_h_pv = evaluate_instrument(item.second, status, x, n());
						if (status != StatusCode::kOk)
							return false;

						// Calculate partial derivative
						if (plus_h_pv != base_pv) {
							partial = (plus_h_pv - base_pv) / h;
						}
					}
					// reset pricing curve
					ch->pricing_curve_ = ch->base_curve_.get();

					trace("Setting matrix A for curve %d, curve pillar %d, "
					      "global pillar %d instrument %d at [%d]\n",
					      curve_num, (int)pillar, rate_num, instrument_num,
					      rate_num * m() + instrument_num);

					// Set value in matrix
					// rate_num is the column
					// instrument_num is the row
					// m() is the number of instruments (rows)
					// A[rate_num * m() + instrument_num] = partial;
					redukti_matrix_set(&X, instrument_num, rate_num, partial);
					trace("A[%d, %d] = %f\n", (int)instrument_num, rate_num, partial);

					instrument_num++;
				}
				rate_num++;
			}
			curve_num++;
		}
		// dump("Jacobian:", A, debug_output_);
		// dump("Vector:", bx, debug_output_);
		return true;
	}

	void set_vector(const double *x, size_t n)
	{
		assert(n == this->n());
		int rate_num = 0;
		for (auto &ch : curve_builder_->curves()) {
			// for each pillar on the curve we
			// compute the partial derivative for each instrument
			unsigned n = ch->n_maturities_;
			ch->set_vector(x + rate_num, n);
			rate_num += n;
		}
	}

	void get_vector(double *x, size_t len)
	{
		assert(len == n());
		int rate_num = 0;
		for (auto &ch : curve_builder_->curves()) {
			// for each pillar on the curve we
			// compute the partial derivative for each instrument
			unsigned n = ch->n_maturities_;
			ch->get_vector(x + rate_num, n);
			rate_num += n;
		}
	}

	void compute_bounds()
	{
		int rate_num = 0;
		for (auto &ch : curve_builder_->curves()) {
			unsigned n = ch->n_maturities_;
			if (ch->curve_type_ == CurveType::CURVE_TYPE_SVENSSON_PARAMETRIC) {
				double *ub = upper_bounds_ + rate_num;
				double *lb = lower_bounds_ + rate_num;
				// TODO These bounds need to be reviewed
				ub[0] = 25.0; // TODO Beta0 should be close to the overnight rate according to some
					      // docs TBC
				lb[0] = 0.0;  // Beta0 must be positive
				ub[1] = HUGE_VAL;
				lb[1] = -HUGE_VAL;
				ub[2] = HUGE_VAL;
				lb[2] = -HUGE_VAL;
				ub[3] = HUGE_VAL;
				lb[3] = -HUGE_VAL;
				lb[4] = 0.0; // tau1 must be positive
				ub[4] = HUGE_VAL;
				lb[5] = 0.0; // tau2 must be positive
				ub[5] = HUGE_VAL;
			}
			rate_num += n;
		}
		apply_bounds();
	}

	virtual void apply_bounds() {}
	virtual bool solve() = 0;
};

class LinearLeastSquaresSolver : public SolverFunction
{
	public:
	LinearLeastSquaresSolver(CurveBuilder *builder, int max_iterations, Allocator *alloc, FILE *debug_output)
	    : SolverFunction(builder, max_iterations, alloc, debug_output)
	{
	}

	bool solution()
	{
		// Now solve.
		// We use a linear square solver - the -1.0 argument says that
		// the solver should estimate the reciprocal condition number of
		// the matrix
		redukti_matrix_t X = {m(), n(), A_};
		redukti_matrix_t c = {m(), 1, bx_}; // NOTE bx must be sized to m() not n()
		double rcond = -1.0;

		if (!redukti_matrix_estimate_rcond(&X, &rcond)) {
			redukti_matrix_dump_to(&X, "Failed to estimate rcond for Matrix", stderr);
			error("Failed to solve the system of equations\n");
			return false;
		}
		if (!redukti_matrix_lsq_solve(&X, &c, rcond /* 0.0000001*/, true)) {
			redukti_matrix_dump_to(&X, "Failed to solve equation: Matrix", stderr);
			redukti_matrix_dump_to(&c, "Failed to solve equation: vector", stderr);
			error("Failed to solve the system of equations\n");
			return false;
		}
		return true;
	}

	void apply_corrections()
	{
		// Now we apply corrections that are in bx to the
		// zero rates
		// dump("Corrections vector:", bx, debug_output_);
		int rate_num = 0;
		for (auto &ch : curve_builder_->curves()) {
			// for each pillar on the curve we
			// compute the partial derivative for each instrument
			unsigned n = ch->n_maturities_;
			ch->set_vector(bx_ + rate_num, n, true);
			rate_num += n;
		}
	}

	bool is_solved(int iter)
	{
		// Check if corrections were good enough that we can stop
		// iterating
		bool iterate = false;
		// NOTE we are only interested in n() elements of bx as
		// they correspond to pillars
		for (size_t i = 0; i < n(); i++) {
			if (std::fabs(bx_[i]) > 1e-10) {
				iterate = true;
				break;
			}
		}
		return iterate;
	}

	bool solve() override final
	{
		for (int iter = 0; iter < max_iterations_; iter++) {
			// fprintf(stderr, "Pricing instruments: iteration
			// %d\n", iter);
			if (!evaluate())
				return false;
			if (!compute_jacobian())
				return false;
			if (!solution())
				return false;
			apply_corrections();
			// Check if corrections were good enough that we can
			// stop iterating
			bool iterate = is_solved(iter);
			if (!iterate)
				break;
		}
		fprintf(stdout, "+");
		return true;
	}
};

class LMDER_Solver : public SolverFunction
{
	public:
	LMDER_Solver(CurveBuilder *builder, int max_iterations, Allocator *alloc, FILE *debug_output)
	    : SolverFunction(builder, max_iterations, alloc, debug_output)
	{
	}

	bool solve() override final
	{
#if 0
		size_t lwa = 5 * this->n() + this->m();
		std::unique_ptr<double[]> wa = std::unique_ptr<double[]>(new double[lwa]);
#endif

		std::unique_ptr<int[]> ipvt = std::unique_ptr<int[]>(new int[this->n()]);
		std::unique_ptr<double[]> wa1 = std::unique_ptr<double[]>(new double[n()]);
		std::unique_ptr<double[]> wa2 = std::unique_ptr<double[]>(new double[n()]);
		std::unique_ptr<double[]> wa3 = std::unique_ptr<double[]>(new double[n()]);
		std::unique_ptr<double[]> wa4 = std::unique_ptr<double[]>(new double[m()]);
		std::unique_ptr<double[]> qtf = std::unique_ptr<double[]>(new double[n()]);
		std::unique_ptr<double[]> diag = std::unique_ptr<double[]>(new double[n()]);

		double tol = std::sqrt(dpmpar(1));
		int info = 0;

		std::unique_ptr<double[]> x = std::unique_ptr<double[]>(new double[this->n()]);
		get_vector(x.get(), this->n());

#if 0
		// we switched to the full api so that we can set max evaluations
		info = lmder1(minpack_lmder_function, this, this->m(), this->n(), x.get(), this->solution_vector(),
			      this->jacobian(), this->m(), tol, ipvt.get(), wa.get(), lwa);
#endif

		double ftol = tol;
		double xtol = tol;
		double gtol = 0.0;
		int maxfev = 10000;
		int mode = 1; // autoscale
		double factor = 100.;
		int nprint = 0;
		int nfev = 0;
		int njev = 0;
		info = lmder(minpack_lmder_function, reinterpret_cast<void *>(this), this->m(), this->n(), x.get(),
			     this->solution_vector(), this->jacobian(), this->m(), ftol, xtol, gtol, maxfev, diag.get(),
			     mode, factor, nprint, &nfev, &njev, ipvt.get(), qtf.get(), wa1.get(), wa2.get(), wa3.get(),
			     wa4.get());

		debug("info = %d, nfev = %d, njev = %d\n", info, nfev, njev);
		set_vector(x.get(), this->n());

		// debug
		// unsigned int instrument_num = 0;
		// for (auto &item : curve_builder_->sorted_instruments()) {
		//  double pv = evaluate_instrument(item.second);
		//  printf("bx[%d] = %f\n", instrument_num, pv);
		//  instrument_num++;
		//}
		fprintf(stdout, "%d", info);
		return info == 1 || info == 2 || info == 3;
	}
};

/* for lmder1 and lmder */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate lmder1/lmder */
int minpack_lmder_function(void *p, int m, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag)
{
	LMDER_Solver *fn = reinterpret_cast<LMDER_Solver *>(p);
	assert((size_t)n == fn->n());
	assert((size_t)m == fn->m());
	assert((size_t)ldfjac == fn->m());

	fn->set_vector(x, n);
	if (iflag == 1) {
		if (!fn->evaluate())
			return -1;
		if (fvec != fn->solution_vector()) {
			std::copy_n(fn->solution_vector(), m, fvec);
		}
	} else if (iflag == 2) {
		if (!fn->compute_jacobian())
			return -1;
		if (fjac != fn->jacobian()) {
			std::copy_n(fn->jacobian(), m * n, fjac);
		}
	}
	return 0;
}

std::unique_ptr<SolverFunction> CurveBuilder::make_solver(const CurveBuilderOptions &options)
{
	switch (options.solver_type) {
	default:
	case SolverType::SOLVER_TYPE_LEVENBERG_MARQUARDT:
		return std::make_unique<LMDER_Solver>(this, options.max_iterations, get_default_allocator(), stderr);
	case SolverType::SOLVER_TYPE_LINEAR_LEAST_SQUARE:
		return std::make_unique<LinearLeastSquaresSolver>(this, options.max_iterations, get_default_allocator(),
								  stderr);
	}
}

// This function builds the base curves and optionally
// also builds the par sensitivities matrix (Jacobian) for
// each curve
StatusCode CurveBuilder::build_curves(const CurveBuilderOptions &options)
{
	auto solver = make_solver(options);
	solver->compute_bounds();
	solver->set_penalty_function(options.penalty_function);
	bool result = solver->solve();
	if (!result) {
		solver->dump_solution_vector("Final solution vector", true);
		return StatusCode::kBTS_SolverFailure;
	}
	if (!options.generate_par_sensitivities)
		return StatusCode::kOk;
	auto &curves = this->curves();
	// First save all the bootstrapped rates as we can use these
	// as our initial guess curve.
	for (auto &ch : curves) {
		ch->save();
	}
	par_sensitivities_by_curve_.clear();
	for (int curve_index = 0; curve_index < instruments_by_curve_.size(); curve_index++) {
		if (curves[curve_index]->curve_type_ == CurveType::CURVE_TYPE_SVENSSON_PARAMETRIC)
			// PAR sensitivities don't make sense when the curve is made up
			// of parameters
			continue;
		// m is the number of instruments
		// and also the rows of the par sensitivities matrix
		int m = instruments_by_curve_[curve_index]->instruments_by_maturity.size();
		// n is the number of pillars in the curve
		// and also the number of columns in the par sensitivities
		// matrix And n therefore matches the number of zero
		// sensitivities that we get
		int n = curves[curve_index]->n_maturities_;
		auto par_sens = std::make_unique<ParSensitivities>(m, n);
		auto par_sens_ptr = par_sens.get();
		par_sensitivities_by_curve_.push_back(std::move(par_sens));

		std::unique_ptr<double[]> up_rates = std::make_unique<double[]>(n);
		std::unique_ptr<double[]> down_rates = std::make_unique<double[]>(n);

		auto &curve_holder = curves[curve_index];

		DynamicRegionAllocator temp_allocator(32 * 1024);
		int row = 0; // tracks instrument number - goes up to m
		for (auto &entry : instruments_by_curve_[curve_index]->instruments_by_maturity) {
			auto &cfinst = entry.second;

			// Backup instrument cashflows as we will
			// use bumped cashflows below, and after that we
			// need to restore the original cashflows
			cfinst->backup_cashflows();

			// Get the par rate that was originally provided
			auto input_curve = get_input_curve((int)curve_index);
			auto par_rate = input_curve.par_rates().values(cfinst->original_instrument_position());

			const double h = 0.00001;
			// Bump par rate up
			double par_rate_up = par_rate + h;

			// Rebuild cashflows using the bumped par rate
			// We do this by just building a temporary instrument
			// instance
			auto tmp_cfinst =
			    build_instrument(input_curve, cfinst->original_instrument_position(), par_rate_up);
			// Copy the cashflows from the temp instance
			cfinst->copy_cashflows(std::move(tmp_cfinst));
			cfinst->build_cashflows(&temp_allocator);

			// By default we have baseline rates as initial guess
			// Run solver
			bool result = solver->solve();
			if (!result) {
				return StatusCode::kBTS_SolverFailure;
			}
			// Copy the up rates
			std::copy(curve_holder->rates_, curve_holder->rates_ + n, up_rates.get());

			// Restore the baseline rates
			for (auto &ch : curves) {
				ch->restore();
			}
			// Bump par rate down
			double par_rate_down = par_rate - h;

			// Rebuild cashflows using the bumped par rate
			// We do this by just building a temporary instrument
			// instance
			tmp_cfinst =
			    build_instrument(input_curve, cfinst->original_instrument_position(), par_rate_down);
			cfinst->copy_cashflows(std::move(tmp_cfinst));
			cfinst->build_cashflows(&temp_allocator);

			// Rebuild curves
			result = solver->solve();
			if (!result) {
				return StatusCode::kBTS_SolverFailure;
			}

			// Copy the down rates
			std::copy(curve_holder->rates_, curve_holder->rates_ + n, down_rates.get());

			// Restore the baseline curves
			for (auto &ch : curves) {
				ch->restore();
			}

			// col represents curve pillar - goes up to n
			for (int col = 0; col < n; col++) {
				double difference = (up_rates[col] - down_rates[col]) / (2 * h);
				par_sens_ptr->set(row, col, difference);
			}

			// Restore the original cashlows as
			// we are done with this instrument
			cfinst->restore_cashflows();
			row++;

			temp_allocator.release();
		}
		par_sens_ptr->dump(stdout);
	}
	return StatusCode::kOk;
}

class BootstrapperImpl : public CurveBuilderService
{
	private:
	std::unique_ptr<LuaWrapper> lua_;
	bool luaok_;

	public:
	BootstrapperImpl(std::string script);
	~BootstrapperImpl() {}
	BootstrapCurvesReply *handle_bootstrap_request(const BootstrapCurvesRequest *request,
						       BootstrapCurvesReply *reply) override final;
};

BootstrapperImpl::BootstrapperImpl(std::string script)
{
	lua_ = std::unique_ptr<LuaWrapper>(new LuaWrapper());
	luaok_ = lua_->load_luascript(script);
}

BootstrapCurvesReply *BootstrapperImpl::handle_bootstrap_request(const BootstrapCurvesRequest *request,
								 BootstrapCurvesReply *reply)
{
	auto header = reply->mutable_header();
	header->set_response_code(StandardResponseCode::SRC_ERROR);
	header->set_response_message(error_message(StatusCode::kInternalError));
	// debug("Request %s\n", request->DebugString().c_str());
	if (!luaok_) {
		header->set_response_message(error_message(StatusCode::kBTS_LuaInitFailure));
		header->set_response_sub_code(StatusCode::kBTS_LuaInitFailure);
		return reply;
	}
	if (!is_valid_date(request->business_date())) {
		header->set_response_message(error_message(StatusCode::kBTS_BadBusinessDate));
		header->set_response_sub_code(StatusCode::kBTS_BadBusinessDate);
		return reply;
	}
	if (request->curve_definitions_size() == 0) {
		header->set_response_message(error_message(StatusCode::kBTS_MissingCurveDefinitions));
		header->set_response_sub_code(StatusCode::kBTS_MissingCurveDefinitions);
		return reply;
	}
	if (!request->has_par_curve_set()) {
		header->set_response_message(error_message(StatusCode::kBTS_MissingParCurves));
		header->set_response_sub_code(StatusCode::kBTS_MissingParCurves);
		return reply;
	}
	if (request->curve_definitions_size() != request->par_curve_set().par_curves_size()) {
		header->set_response_message(error_message(StatusCode::kBTS_DefinitonsAndCurvesMismatch));
		header->set_response_sub_code(StatusCode::kBTS_DefinitonsAndCurvesMismatch);
		return reply;
	}
	CurveDefinitionProviderImpl defn_provider;
	for (int i = 0; i < request->curve_definitions_size(); i++) {
		if (request->curve_definitions(i).curve_type() == CurveType::CURVE_TYPE_SVENSSON_PARAMETRIC) {
			if (request->curve_definitions(i).maturity_generation_rule() !=
			    MaturityGenerationRule::MATURITY_GENERATION_RULE_SVENSSON) {
				header->set_response_message(error_message(StatusCode::kBTS_BadMaturityGenerationRule));
				header->set_response_sub_code(StatusCode::kBTS_BadInterpolatedOn);
				return reply;
			} else if (request->curve_definitions(i).interpolated_on() != IRRateType::ZERO_RATE) {
				header->set_response_message(error_message(StatusCode::kBTS_BadInterpolatedOn));
				header->set_response_sub_code(StatusCode::kBTS_BadInterpolatedOn);
				return reply;
			}
		}
		if (request->curve_definitions(i).maturity_generation_rule() ==
			MaturityGenerationRule::MATURITY_GENERATION_RULE_FIXED_TENORS &&
		    request->curve_definitions(i).tenors_size() == 0) {
			header->set_response_message(error_message(StatusCode::kBTS_MissingFixedTenors));
			header->set_response_sub_code(StatusCode::kBTS_MissingFixedTenors);
			return reply;
		}
		if (request->curve_definitions(i).interpolated_on() == IRRateType::ZERO_RATE ||
		    request->curve_definitions(i).interpolated_on() == IRRateType::DISCOUNT_FACTOR) {
			defn_provider.add(&request->curve_definitions(i));
		} else {
			header->set_response_message(error_message(StatusCode::kBTS_BadInterpolatedOn));
			header->set_response_sub_code(StatusCode::kBTS_BadInterpolatedOn);
			return reply;
		}
	}
	CurveBuilder curve_builder(*lua_, request->business_date(), &defn_provider);
	if (!curve_builder.add_input_curves(&request->par_curve_set())) {
		header->set_response_message(curve_builder.last_error_message());
		header->set_response_sub_code(curve_builder.last_error_code());
		return reply;
	}
	if (!curve_builder.prepare_instruments()) {
		header->set_response_message(curve_builder.last_error_message());
		header->set_response_sub_code(curve_builder.last_error_code());
		return reply;
	}
	CurveBuilderOptions options;
	options.generate_par_sensitivities = request->generate_par_sensitivities();
	options.solver_type = request->solver_type();
	options.max_iterations = request->max_solver_iterations() != 0 ? request->max_solver_iterations() : 20;
	debug("Starting curve build: Levenberg_marquardt solver: %d, max_iterations: "
	      "%d\n",
	      options.solver_type == SolverType::SOLVER_TYPE_LEVENBERG_MARQUARDT, options.max_iterations);
	auto status = curve_builder.build_curves(options);
	debug("Curve build completed %s\n", status != StatusCode::kOk ? "with ERROR" : "successfully");
	if (status != StatusCode::kOk) {
		header->set_response_message(error_message(status));
		header->set_response_sub_code(status);
		return reply;
	}
	for (int i = 0; i < request->curve_definitions_size(); i++) {
		auto &def = request->curve_definitions(i);
		auto curveholder = curve_builder.get_curve_by_definition_id(def.id());
		ParSensitivities *parsens = options.generate_par_sensitivities
						? curve_builder.get_sensitivities_by_curve_definition_id(def.id())
						: nullptr;
		if (curveholder == nullptr || (options.generate_par_sensitivities && parsens == nullptr &&
					       def.curve_type() == CurveType::CURVE_TYPE_INTERPOLATED)) {
			error("PAR sensitivities requested but failed to generate those");
			return reply;
		}
		ZeroCurve *zc = reply->add_curves();
		for (int i = 0; i < curveholder->n_maturities_; i++) {
			zc->add_maturities(curveholder->maturities_[i]);
			zc->add_values(curveholder->rates_[i]);
		}
		zc->set_curve_definition_id(def.id());
		if (parsens) {
			ZeroCurveParSensitivities *zc_par_sens = reply->add_par_sensitivities();
			zc_par_sens->set_num_instruments(parsens->m());
			zc_par_sens->set_num_maturities(parsens->n());
			zc_par_sens->set_curve_definition_id(def.id());
			for (int row = 0; row < parsens->m(); row++) {
				for (int col = 1; col < parsens->n(); col++) {
					auto v = parsens->value(row, col);
					if (std::fabs(v) > 1e-10) {
						uint32_t return_row = row;
						uint32_t return_col = col;
						uint32_t return_index = (uint32_t)((return_col << 16) | return_row);
						zc_par_sens->mutable_values()->insert(
						    google::protobuf::MapPair<google::protobuf::uint32, double>(
							return_index, v));
					}
				}
			}
		}
	}
	header->set_response_code(StandardResponseCode::SRC_OK);
	header->set_response_message("");
	return reply;
}

std::unique_ptr<CurveBuilderService> get_curve_builder_service(std::string script)
{
	return std::make_unique<BootstrapperImpl>(script);
}

int test_bootstrap()
{
	// CODE not released
	printf("Bootstrap Tests skipped\n");
	return 0;
}

} // namespace redukti
