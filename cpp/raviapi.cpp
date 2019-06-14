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
#include <cashflow.pb.h>
#include <curve.pb.h>
#include <enums.pb.h>
#include <schedule.pb.h>

#include <autodiff.h>
#include <bootstrap.h>
#include <buffer.h>
#include <calendars.h>
#include <cashflow.h>
#include <cashflow_pricing.h>
#include <converters.h>
#include <curve.h>
#include <datasource.h>
#include <date.h>
#include <dayfractions.h>
#include <fixings.h>
#include <hashtable.h>
#include <index.h>
#include <internal/cashflow_internal.h>
#include <internal/cashflow_pricing_internal.h>
#include <interpolators.h>
#include <logger.h>
#include <schedule.h>
#include <statistics.h>

#include <map>
#include <sstream>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif
#include <lauxlib.h>
#include <lua.h>
#include <lualib.h>
#ifdef __cplusplus
}
#endif

#include <assert.h>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <ctime>

extern "C" {

using namespace redukti;
using namespace google::protobuf;

static const char *Type_CFCollection = "Redukti_CFCollection";
static const char *Type_Index = "Redukti_Index";
static const char *Type_DayFraction = "Redukti_DayFraction";
static const char *Type_Calendar = "Redukti_Calendar";
static const char *Type_ADouble = "Redukti_ADouble";
static const char *Type_Interpolator = "Redukti_Interpolator";
static const char *Type_IRCurve = "Redukti_IRCurve";
static const char *Type_CurveMapper = "Redukti_CurveMapper";
static const char *Type_ValuationContext = "Redukti_ValuationContext";
static const char *Type_PricingCashflows = "Redukti_PricingCashflows";
static const char *Type_CurveProvider = "Redukti_CurveProvider";
static const char *Type_PricingResult = "Redukti_PricingResult";
static const char *Type_ZeroCurveSet = "Redukti_ZeroCurveSet";
static const char *Type_FixingService = "Redukti_FixingService";

#define test_Type_CFCollection(L, idx) ((CFCollectionHolder *)luaL_testudata(L, idx, Type_CFCollection))
#define check_Type_CFCollection(L, idx) ((CFCollectionHolder *)luaL_checkudata(L, idx, Type_CFCollection))
#define test_Type_Index(L, idx) ((IndexHolder *)luaL_testudata(L, idx, Type_Index))
#define check_Type_Index(L, idx) ((IndexHolder *)luaL_checkudata(L, idx, Type_Index))
#define test_Type_DayFraction(L, idx) ((DayFractionHolder *)luaL_testudata(L, idx, Type_DayFraction))
#define check_Type_DayFraction(L, idx) ((DayFractionHolder *)luaL_checkudata(L, idx, Type_DayFraction))
#define test_Type_Calendar(L, idx) ((CalendarHolder *)luaL_testudata(L, idx, Type_Calendar))
#define check_Type_Calendar(L, idx) ((CalendarHolder *)luaL_checkudata(L, idx, Type_Calendar))
#define test_Type_ADouble(L, idx) ((struct redukti_adouble_t *)luaL_testudata(L, idx, Type_ADouble))
#define check_Type_ADouble(L, idx) ((struct redukti_adouble_t *)luaL_checkudata(L, idx, Type_ADouble))
#define test_Type_Interpolator(L, idx) ((InterpolatorHolder *)luaL_testudata(L, idx, Type_Interpolator))
#define check_Type_Interpolator(L, idx) ((InterpolatorHolder *)luaL_checkudata(L, idx, Type_Interpolator))
#define test_Type_IRCurve(L, idx) ((IRCurveHolder *)luaL_testudata(L, idx, Type_IRCurve))
#define check_Type_IRCurve(L, idx) ((IRCurveHolder *)luaL_checkudata(L, idx, Type_IRCurve))
#define test_Type_CurveMapper(L, idx) ((CurveMapperHolder *)luaL_testudata(L, idx, Type_CurveMapper))
#define check_Type_CurveMapper(L, idx) ((CurveMapperHolder *)luaL_checkudata(L, idx, Type_CurveMapper))
#define test_Type_ValuationContext(L, idx) ((ValuationContextHolder *)luaL_testudata(L, idx, Type_ValuationContext))
#define check_Type_ValuationContext(L, idx) ((ValuationContextHolder *)luaL_checkudata(L, idx, Type_ValuationContext))
#define test_Type_PricingCashflows(L, idx) ((PricingCashflowsHolder *)luaL_testudata(L, idx, Type_PricingCashflows))
#define check_Type_PricingCashflows(L, idx) ((PricingCashflowsHolder *)luaL_checkudata(L, idx, Type_PricingCashflows))
#define test_Type_CurveProvider(L, idx) ((CurveProviderHolder *)luaL_testudata(L, idx, Type_CurveProvider))
#define check_Type_CurveProvider(L, idx) ((CurveProviderHolder *)luaL_checkudata(L, idx, Type_CurveProvider))
#define test_Type_PricingResult(L, idx) ((PricingResult *)luaL_testudata(L, idx, Type_PricingResult))
#define check_Type_PricingResult(L, idx) ((PricingResult *)luaL_checkudata(L, idx, Type_PricingResult))
#define test_Type_ZeroCurveSet(L, idx) ((ZeroCurveSetHolder *)luaL_testudata(L, idx, Type_ZeroCurveSet))
#define check_Type_ZeroCurveSet(L, idx) ((ZeroCurveSetHolder *)luaL_checkudata(L, idx, Type_ZeroCurveSet))
#define test_Type_FixingService(L, idx) ((FixingServiceHolder *)luaL_testudata(L, idx, Type_FixingService))
#define check_Type_FixingService(L, idx) ((FixingServiceHolder *)luaL_checkudata(L, idx, Type_FixingService))

/* garbage collected */
struct CFCollectionHolder {
	CFCollection *cashflows;
};

/* garbage collected */
struct InterpolatorHolder {
	std::unique_ptr<double> x_values;
	std::unique_ptr<double> y_values;
	std::unique_ptr<Interpolator, Deleter<Interpolator>> interp;

	InterpolatorHolder(std::unique_ptr<double> x, std::unique_ptr<double> y,
			   std::unique_ptr<Interpolator, Deleter<Interpolator>> i)
	    : x_values(std::move(x)), y_values(std::move(y)), interp(std::move(i))
	{
	}
};

/* garbage collected */
struct IRCurveHolder : public CurveReference {
	std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve;

	IRCurveHolder(std::unique_ptr<YieldCurve, Deleter<YieldCurve>> c) : curve(std::move(c)) {}
	virtual ~IRCurveHolder() {}
	virtual YieldCurve *get() const noexcept override final { return curve.get(); }
};

struct IndexHolder {
	InterestRateIndex *index;
};

struct DayFractionHolder {
	const DayFraction *fraction;
};

struct CalendarHolder {
	const Calendar *calendar;
};

struct CurveMapperHolder {
	SimpleCurveMapper mapper;
};

struct CurveProviderHolder {
	SimpleCurveProvider provider;
};

struct ValuationContextHolder {
	ValuationContextImpl context;

	ValuationContextHolder(Date business_date, int order, FixingDataService *fixing_service)
	    : context(business_date, order, fixing_service)
	{
	}
};

struct PricingCashflowsHolder {
	DynamicRegionAllocator fallback_allocator;
	StackRegionWithFallbackAllocator<2048> stackregion_allocator;

	Cashflows *cashflows;

	PricingCashflowsHolder() : stackregion_allocator(&fallback_allocator), cashflows(nullptr) {}
};

struct PricingResult {
	char buf[1 * 1024 * 1024];
	DynamicRegionAllocator sensitivities_allocator;
	FixedRegionAllocator pricing_allocator;

	double value;
	Sensitivities sensitivities;
	StatusCode status;

	PricingResult()
	    : sensitivities_allocator(), pricing_allocator(buf, buf + sizeof buf), value(0.0),
	      sensitivities(&sensitivities_allocator), status(StatusCode::kNotAvailable)
	{
	}
};

struct ZeroCurveSetHolder {
	Arena arena;
	BootstrapCurvesRequest *request;
	ZeroCurveSet *zero_curve_set;

	ZeroCurveSetHolder() : request(nullptr), zero_curve_set(nullptr) {}
};

struct FixingServiceHolder {
	FixingDataService fixing_service;
};

// Utility to extract a number field from a table
// Table must be on top of the stack
// If successful, *result will be set
// Else *result will be set to supplied default
// Return value indicates whether number was obtained (true) or
// default value used (false)
static bool table_get_number(lua_State *L, int idx, const char *key, lua_Number *result, lua_Number default_value)
{
	bool rc = false;
	idx = lua_absindex(L, idx);
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	if (!lua_isnumber(L, -1))
		*result = default_value;
	else {
		*result = lua_tonumber(L, -1);
		rc = true;
	}
	lua_pop(L, 1); /* remove number */
	return rc;
}

// Retrieves an integer value from a Lua table / Ravi array
// from specified index. If successful *val is initialized
// to retrieved value else to 0.
// Returns true for success else false
static bool array_get_integer(lua_State *L, int table_stack_index, int index, lua_Integer *val)
{
	bool rc = false;
	lua_rawgeti(L, table_stack_index, index);
	if (lua_isnumber(L, -1)) {
		*val = lua_tointeger(L, -1);
		rc = true;
	} else {
		debug("Error: expected integer at %d\n", index);
		*val = 0;
	}
	lua_pop(L, 1);
	return rc;
}

// Utility to extract a integer field from a table
// Table must be on top of the stack
// If successful *result is set and true returned
static bool table_get_integer(lua_State *L, int idx, const char *key, int *result)
{
	bool rc = false;
	idx = lua_absindex(L, idx);
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	if (!lua_isnumber(L, -1)) {
		*result = 0;
		goto done;
	}
	*result = (int)lua_tointeger(L, -1);
	rc = true;
done:
	lua_pop(L, 1); /* remove number */
	return rc;
}

static bool table_get_lua_Integer(lua_State *L, int idx, const char *key, lua_Integer *result)
{
	bool rc = false;
	idx = lua_absindex(L, idx);
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	if (!lua_isnumber(L, -1)) {
		*result = 0;
		goto done;
	}
	*result = lua_tointeger(L, -1);
	rc = true;
done:
	lua_pop(L, 1); /* remove number */
	return rc;
}

// Utility to extract a boolean field from a table
// Table must be on top of the stack
// If successful *result is set and true returned
static bool table_get_boolean(lua_State *L, int idx, const char *key, bool *result)
{
	bool rc = false;
	idx = lua_absindex(L, idx);
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	if (!lua_isboolean(L, -1)) {
		*result = false;
		goto done;
	}
	*result = (int)lua_toboolean(L, -1);
	rc = true;
done:
	lua_pop(L, 1); /* remove value */
	return rc;
}

// Utility to extract a string field from a table
// Table must be on top of the stack
// The string is copied to the suppled buffer
//
static bool table_get_string(lua_State *L, int idx, const char *key, char *buffer, size_t buffer_length)
{
	bool rc = false;
	idx = lua_absindex(L, idx);
	if (buffer_length == 0)
		goto done;
	*buffer = 0;
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	if (!lua_isstring(L, -1))
		goto done;
	snprintf(buffer, buffer_length, "%s", lua_tostring(L, -1));
	rc = true;
done:
	lua_pop(L, 1); /* remove string */
	return rc;
}

// Utility for getting a named number array member
// from a table
// Table must be at the top of the stack
// Caller must supply a values buffer and
// the length of the buffer
// Returns actual length loaded or 0 on error
static size_t table_get_doublearray(lua_State *L, int idx, const char *key, double *values, const size_t len)
{
	size_t rawlen = 0;
	idx = lua_absindex(L, idx);
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	bool errored = false;
	if (!lua_istable(L, -1))
		errored = true;
	if (!errored) {
		rawlen = lua_rawlen(L, -1);
		if (rawlen > len) {
			rawlen = len;
		}
	}
	if (!errored) {
		for (size_t j = 0; j < rawlen; j++) {
			lua_rawgeti(L, -1, j + 1);
			if (lua_isnumber(L, -1))
				values[j] = lua_tonumber(L, -1);
			else
				errored = true;
			lua_pop(L, 1);
			if (errored)
				break;
		}
	}
	lua_pop(L, 1); /* remove table[key] */
	return errored ? 0 : rawlen;
}

// Utility for getting a named integer array member
// from a table
// Table must be at the top of the stack
// Caller must supply a values buffer and
// the length of the buffer
// Returns actual length loaded or 0 on error
static size_t table_get_intarray(lua_State *L, int idx, const char *key, int32_t *values, const size_t len)
{
	size_t rawlen = 0;
	idx = lua_absindex(L, idx);
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	bool errored = false;
	if (!lua_istable(L, -1))
		errored = true;
	if (!errored) {
		rawlen = lua_rawlen(L, -1);
		if (rawlen > len) {
			rawlen = len;
		}
	}
	if (!errored) {
		for (size_t j = 0; j < rawlen; j++) {
			lua_rawgeti(L, -1, j + 1);
			if (lua_isinteger(L, -1))
				values[j] = (int32_t)lua_tointeger(L, -1);
			else
				errored = true;
			lua_pop(L, 1);
			if (errored)
				break;
		}
	}
	lua_pop(L, 1); /* remove table[key] */
	return errored ? 0 : rawlen;
}

// Utility for getting a named integer array member
// from a table
static size_t table_get_arraylength(lua_State *L, int idx, const char *key)
{
	size_t len = 0;
	idx = lua_absindex(L, idx);
	lua_pushstring(L, key);
	lua_gettable(L, idx); /* get table[key] */
	if (!lua_istable(L, -1)) {
		len = 0;
	} else {
		len = lua_rawlen(L, -1);
	}
	lua_pop(L, 1); /* remove table[key] */
	return len;
}

// reduki.date(day,month,year) - takes day, month, year and returns a date
// reduki.date(str) - converts string to date
// reduki.date{ d,m,y } - takes table and converts to date
// reduki.date() returns today's date
static int make_date(lua_State *L)
{
	Date date;
	lua_Integer d, m, y;
	if (lua_gettop(L) == 3) {
		// d,m,y
		d = luaL_checkinteger(L, 1);
		m = luaL_checkinteger(L, 2);
		y = luaL_checkinteger(L, 3);
	} else if (lua_istable(L, 1)) {
		if (!array_get_integer(L, 1, 1, &d) || !array_get_integer(L, 1, 2, &m) ||
		    !array_get_integer(L, 1, 3, &y)) {
			luaL_error(L, "Expect 3 int elements\n");
			return 0;
		}
	} else if (lua_isstring(L, 1)) {
		const char *s = lua_tostring(L, 1);
		if (!parse_date(s, &date)) {
			luaL_error(L, "Invalid argument, failed to parse date %s", s);
			return 0;
		}
		lua_pushinteger(L, date);
		return 1;
	} else {
		time_t now = time(nullptr);
		struct tm *datetime = localtime(&now);
		y = datetime->tm_year + 1900;
		m = datetime->tm_mon + 1;
		d = datetime->tm_mday;
	}
	date = make_date(d, m, y);
	lua_pushinteger(L, date);
	return 1;
}

// redukti.date_parts - takes a date and returns day, month, year
static int date_parts(lua_State *L)
{
	lua_Integer date = luaL_checkinteger(L, 1);
	luaL_argcheck(L, is_valid_date(date), 1, "invalid date");
	Date dt = (Date)date;
	auto ymd = date_components(dt);
	lua_pushinteger(L, ymd.d);
	lua_pushinteger(L, ymd.m);
	lua_pushinteger(L, ymd.y);
	return 3;
}

// add period to date
// redukti.addperiod(date, period) -> date
// period is a string of the form [-]n[D|M|W|Y|T]
static int add_period(lua_State *L)
{
	Date d = luaL_checkinteger(L, 1);
	luaL_argcheck(L, is_valid_date(d), 1, "invalid date");
	const char *period = luaL_checkstring(L, 2);
	Period p;
	luaL_argcheck(L, get_default_converter()->period_from_string(period, &p), 2, "invalid period");
	d = add(d, p);
	lua_pushinteger(L, d);
	return 1;
}

static void parse_calendars(const char *startp, const char *endp, std::vector<BusinessCenter> &centers)
{
	buffer_ptr<char> buf(80);
	std::vector<const char *> tokens;
	parse_delimited(startp, endp, tokens, buf);
	for (int i = 0; i < tokens.size(); i++) {
		BusinessCenter center = get_default_converter()->business_center_from_string(tokens[i]);
		if (center == BusinessCenter::BUSINESS_CENTER_UNSPECIFIED)
			continue;
		centers.push_back(center);
	}
}

// dy.make_schedule - takes a table as input
// returns 4 tables - adjusted start dates, adjusted end dates, payment flags,
// payment date
static int build_schedule(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1, "table expected");

	Date start = 0, end = 0, first_reg_start = 0, last_reg_end = 0, first_pay = 0, last_reg_pay = 0;
	char roll_convention[80] = {0}, term[80] = {0}, pay_freq[80] = {0}, calc_freq[80] = {0};
	char calc_conv[80] = {0}, pay_conv[80] = {0};
	char calc_calendars[80] = {0}, pay_calendars[80] = {0};
	char stub_type[80] = {0};

	table_get_integer(L, 1, "effective_date", &start);
	table_get_integer(L, 1, "termination_date", &end);
	table_get_integer(L, 1, "first_regular_period_start_date", &first_reg_start);
	table_get_integer(L, 1, "last_regular_period_end_date", &last_reg_end);
	table_get_integer(L, 1, "first_payment_date", &first_pay);
	table_get_integer(L, 1, "last_regular_payment_date", &last_reg_pay);
	table_get_string(L, 1, "roll_convention", roll_convention, sizeof roll_convention);
	table_get_string(L, 1, "term", term, sizeof term);
	table_get_string(L, 1, "payment_frequency", pay_freq, sizeof pay_freq);
	table_get_string(L, 1, "calculation_frequency", calc_freq, sizeof calc_freq);
	table_get_string(L, 1, "calculation_day_convention", calc_conv, sizeof calc_conv);
	table_get_string(L, 1, "payment_day_convention", pay_conv, sizeof pay_conv);
	table_get_string(L, 1, "calculation_business_centers", calc_calendars, sizeof calc_calendars);
	table_get_string(L, 1, "payment_business_centers", pay_calendars, sizeof pay_calendars);
	table_get_string(L, 1, "stub_type", stub_type, sizeof stub_type);

	luaL_argcheck(L, is_valid_date(start), 1, "invalid effective date");

	std::vector<BusinessCenter> payment_centers;
	std::vector<BusinessCenter> calc_centers;
	parse_calendars(std::begin(pay_calendars), std::end(pay_calendars), payment_centers);
	parse_calendars(std::begin(calc_calendars), std::end(calc_calendars), calc_centers);

	luaL_argcheck(L, payment_centers.size() > 0, 1, "payment calendars invalid");

	Period payFreq;
	luaL_argcheck(L, get_default_converter()->period_from_string(pay_freq, &payFreq), 1,
		      "payment frequency invalid");

	StatusCode status_code = StatusCode::kOk;
	for (;;) {
		ScheduleParameters parms;

		parms.set_effective_date(start);

		for (int i = 0; i < payment_centers.size(); i++)
			parms.add_payment_calendars(payment_centers[i]);
		parms.set_payment_frequency(
		    get_default_converter()->tenor_from_period_unit_and_len(payFreq.units(), payFreq.length()));
		parms.set_calculation_frequency(parms.payment_frequency());

		if (is_valid_date(end)) {
			parms.set_termination_date(end);
		} else if (term[0]) {
			Period tradeterm;
			if (get_default_converter()->period_from_string(term, &tradeterm))
				parms.set_term(get_default_converter()->tenor_from_period_unit_and_len(
				    tradeterm.units(), tradeterm.length()));
		}

		if (is_valid_date(first_reg_start)) {
			parms.set_first_regular_period_start_date(first_reg_start);
		}

		if (is_valid_date(last_reg_end)) {
			parms.set_last_regular_period_end_date(last_reg_end);
		}

		if (is_valid_date(first_pay)) {
			parms.set_first_payment_date(first_pay);
		}

		if (is_valid_date(last_reg_pay)) {
			parms.set_last_regular_payment_date(last_reg_pay);
		}

		if (roll_convention[0]) {
			RollConvention rc = get_default_converter()->roll_convention_from_string(roll_convention);
			parms.set_roll_convention(rc);
		}

		if (calc_freq[0]) {
			Period calcFreq;
			if (get_default_converter()->period_from_string(calc_freq, &calcFreq))
				parms.set_calculation_frequency(get_default_converter()->tenor_from_period_unit_and_len(
				    calcFreq.units(), calcFreq.length()));
		}

		if (calc_centers.size() > 0) {
			for (int i = 0; i < calc_centers.size(); i++)
				parms.add_period_calendars(calc_centers[i]);
		} else {
			for (int i = 0; i < payment_centers.size(); i++)
				parms.add_period_calendars(payment_centers[i]);
		}
		if (pay_conv[0]) {
			BusinessDayConvention bdc =
			    get_default_converter()->business_day_convention_from_string(pay_conv);
			parms.set_payment_convention(bdc);
			parms.set_period_convention(bdc);
		}
		if (calc_conv[0]) {
			BusinessDayConvention bdc =
			    get_default_converter()->business_day_convention_from_string(calc_conv);
			parms.set_period_convention(bdc);
		}
		if (stub_type[0]) {
			if (stub_type[0] == 'f' || stub_type[0] == 'F')
				parms.set_stub_location(StubLocation::SHORT_FRONT_STUB);
			else if (stub_type[0] == 'b' || stub_type[0] == 'B')
				parms.set_stub_location(StubLocation::SHORT_BACK_STUB);
		}
		size_t len = 0;
		Schedule schedule;
		status_code = build_schedule(parms, schedule);
		if (status_code != ResponseSubCode::kOk)
			break;
		ravi_create_integer_array(L, schedule.adjusted_start_dates_size(), 0);
		lua_Integer *startarray = ravi_get_integer_array_rawdata(L, lua_gettop(L), &len);
		assert(len == schedule.adjusted_start_dates_size() + 1);
		ravi_create_integer_array(L, schedule.adjusted_start_dates_size(), 0);
		lua_Integer *endarray = ravi_get_integer_array_rawdata(L, lua_gettop(L), &len);
		ravi_create_integer_array(L, schedule.adjusted_start_dates_size(), 0);
		lua_Integer *payarray = ravi_get_integer_array_rawdata(L, lua_gettop(L), &len);
		for (unsigned int i = 0; i < schedule.adjusted_start_dates_size(); i++) {
			startarray[i + 1] = schedule.adjusted_start_dates(i);
			endarray[i + 1] = schedule.adjusted_end_dates(i);
			payarray[i + 1] = schedule.adjusted_payment_dates(i);
		}
		const char *has_stubs = "n";
		if (schedule.has_front_stub() && schedule.has_back_stub())
			has_stubs = "a";
		else if (schedule.has_front_stub())
			has_stubs = "f";
		else if (schedule.has_back_stub())
			has_stubs = "b";
		lua_pushstring(L, has_stubs);
		return 4;
	}
	luaL_error(L, "failed to build schedule: %s", error_message(status_code));
	return 0;
}

static bool build_cashflow_floatingperiod(lua_State *L, CFFloatingPeriod *period)
{
	double notional = 0.0;
	Date accrual_start_date, accrual_end_date;
	char index[80] = {0};
	char tenor[10] = {0};
	double spread = 0.0;
	IsdaIndex index1 = ISDA_INDEX_UNSPECIFIED;
	IsdaIndex index2 = ISDA_INDEX_UNSPECIFIED;
	Tenor tenor1 = TENOR_UNSPECIFIED;
	Tenor tenor2 = TENOR_UNSPECIFIED;
	table_get_number(L, -1, "notional", &notional, 0.0);
	if (notional <= 0.0)
		return false;
	table_get_number(L, -1, "spread", &spread, 0.0);
	if (table_get_string(L, -1, "index", index, sizeof index))
		index1 = get_default_converter()->isda_index_from_string(index);
	if (index1 == ISDA_INDEX_UNSPECIFIED)
		return false;
	if (table_get_string(L, -1, "index2", index, sizeof index))
		index2 = get_default_converter()->isda_index_from_string(index);
	Period p;
	if (table_get_string(L, -1, "tenor", tenor, sizeof tenor)) {
		if (get_default_converter()->period_from_string(tenor, &p))
			tenor1 = get_default_converter()->tenor_from_period_unit_and_len(p.units(), p.length());
	}
	if (tenor1 == TENOR_UNSPECIFIED)
		return false;
	if (table_get_string(L, -1, "tenor2", tenor, sizeof tenor)) {
		if (get_default_converter()->period_from_string(tenor, &p))
			tenor2 = get_default_converter()->tenor_from_period_unit_and_len(p.units(), p.length());
	}
	if (tenor2 != TENOR_UNSPECIFIED && index2 == IsdaIndex::ISDA_INDEX_UNSPECIFIED)
		return false;
	table_get_integer(L, -1, "accrual_start_date", &accrual_start_date);
	table_get_integer(L, -1, "accrual_end_date", &accrual_end_date);
	if (is_valid_date(accrual_start_date))
		period->set_accrual_start_date(accrual_start_date);
	else
		return false;
	if (is_valid_date(accrual_end_date))
		period->set_accrual_end_date(accrual_end_date);
	else
		return false;
	period->set_notional(notional);
	period->set_spread(spread);
	period->set_index(index1);
	period->set_tenor(tenor1);
	period->set_index2(index2);
	period->set_tenor2(tenor2);
	return true;
}

static bool build_cashflow_fra(lua_State *L, CFCollection *cashflows, CFStream *stream, CFFra *fra,
			       CFFloatingPeriod *period)
{
	char currency_str[10] = {0};
	char day_count_fraction_str[40] = {0};
	Currency currency = Currency::CURRENCY_UNSPECIFIED;
	DayCountFraction day_count_fraction = DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED;
	Date payment_date;
	double fixed_rate = 0.0;
	table_get_integer(L, -1, "payment_date", &payment_date);
	if (!is_valid_date(payment_date))
		return false;
	if (table_get_string(L, -1, "currency", currency_str, sizeof currency_str))
		currency = get_default_converter()->currency_from_string(currency_str);
	if (currency == Currency::CURRENCY_UNSPECIFIED)
		return false;
	if (table_get_string(L, -1, "day_count_fraction", day_count_fraction_str, sizeof day_count_fraction_str))
		day_count_fraction = get_default_converter()->day_count_fraction_from_string(day_count_fraction_str);
	if (day_count_fraction == DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED)
		return false;
	table_get_number(L, -1, "fixed_rate", &fixed_rate, 0.0);
	if (fixed_rate <= 0)
		return false;
	fra->set_currency(currency);
	fra->set_day_count_fraction(day_count_fraction);
	fra->set_payment_date(payment_date);
	fra->set_fixed_rate(fixed_rate);
	return build_cashflow_floatingperiod(L, period);
}

static bool build_cashflow_floating(lua_State *L, CFCollection *cashflows, CFStream *stream, CFFloating *floating)
{
	int top = lua_gettop(L);
	char compounding_method_str[80] = {0};
	char currency_str[10] = {0};
	char day_count_fraction_str[40] = {0};
	CompoundingMethod compounding_method = CompoundingMethod::COMPOUNDING_METHOD_NONE;
	Currency currency = Currency::CURRENCY_UNSPECIFIED;
	DayCountFraction day_count_fraction = DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED;
	if (table_get_string(L, -1, "compounding_method", compounding_method_str, sizeof compounding_method_str))
		compounding_method = get_default_converter()->compounding_method_from_string(compounding_method_str);
	if (table_get_string(L, -1, "currency", currency_str, sizeof currency_str))
		currency = get_default_converter()->currency_from_string(currency_str);
	if (currency == Currency::CURRENCY_UNSPECIFIED)
		return false;
	if (table_get_string(L, -1, "day_count_fraction", day_count_fraction_str, sizeof day_count_fraction_str))
		day_count_fraction = get_default_converter()->day_count_fraction_from_string(day_count_fraction_str);
	if (day_count_fraction == DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED)
		return false;
	Date payment_date;
	table_get_integer(L, -1, "payment_date", &payment_date);
	if (!is_valid_date(payment_date))
		return false;
	floating->set_currency(currency);
	floating->set_compounding_method(compounding_method);
	floating->set_day_count_fraction(day_count_fraction);
	floating->set_payment_date(payment_date);
	bool result = true;
	int count = 0;
	lua_pushstring(L, "periods");
	lua_gettable(L, -2);
	if (lua_type(L, -1) == LUA_TTABLE) {
		for (int i = 1; i <= lua_rawlen(L, -1); i++) {
			lua_rawgeti(L, -1, i);
			if (lua_type(L, -1) == LUA_TTABLE) {
				auto period = floating->add_floating_periods();
				if (!build_cashflow_floatingperiod(L, period))
					result = false;
				else
					count++;
			}
			lua_pop(L, 1);
		}
	} else {
		result = false;
	}
	lua_pop(L, 1);
	assert(top == lua_gettop(L));
	return result && count > 0;
}

static bool build_cashflow_ois(lua_State *L, CFCollection *cashflows, CFStream *stream, CFOis *ois)
{
	char day_count_fraction_str[40] = {0};
	DayCountFraction day_count_fraction = DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED;
	if (table_get_string(L, -1, "day_count_fraction", day_count_fraction_str, sizeof day_count_fraction_str))
		day_count_fraction = get_default_converter()->day_count_fraction_from_string(day_count_fraction_str);
	if (day_count_fraction == DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED)
		return false;
	Date payment_date;
	table_get_integer(L, -1, "payment_date", &payment_date);
	if (!is_valid_date(payment_date))
		return false;
	double notional = 0.0;
	Date accrual_start_date, accrual_end_date;
	char index[80] = {0};
	IsdaIndex index1 = ISDA_INDEX_UNSPECIFIED;
	table_get_number(L, -1, "notional", &notional, 0.0);
	if (notional <= 0.0)
		return false;
	if (table_get_string(L, -1, "index", index, sizeof index))
		index1 = get_default_converter()->isda_index_from_string(index);
	if (index1 == IsdaIndex::ISDA_INDEX_UNSPECIFIED)
		return false;
	table_get_integer(L, -1, "accrual_start_date", &accrual_start_date);
	table_get_integer(L, -1, "accrual_end_date", &accrual_end_date);
	if (is_valid_date(accrual_start_date))
		ois->set_accrual_start_date(accrual_start_date);
	else
		return false;
	if (is_valid_date(accrual_end_date))
		ois->set_accrual_end_date(accrual_end_date);
	else
		return false;
	ois->set_notional(notional);
	ois->set_index(index1);
	ois->set_payment_date(payment_date);
	ois->set_day_count_fraction(day_count_fraction);
	return true;
}

static bool build_cashflow_simple(lua_State *L, CFCollection *cashflows, CFStream *stream, CFSimple *simple)
{
	char currency_str[10] = {0};
	Currency currency = Currency::CURRENCY_UNSPECIFIED;
	if (table_get_string(L, -1, "currency", currency_str, sizeof currency_str))
		currency = get_default_converter()->currency_from_string(currency_str);
	if (currency == Currency::CURRENCY_UNSPECIFIED)
		return false;
	Date payment_date;
	table_get_integer(L, -1, "payment_date", &payment_date);
	if (!is_valid_date(payment_date))
		return false;
	double amount = 0.0;
	table_get_number(L, -1, "amount", &amount, 0.0);
	simple->set_currency(currency);
	simple->set_amount(amount);
	simple->set_payment_date(payment_date);
	return true;
}

static bool build_cashflow_stream(lua_State *L, CFCollection *cashflows, CFStream *stream)
{
	bool result = true;
	int count = 0;
	int top = lua_gettop(L);
	size_t n = lua_rawlen(L, -1);
	for (int i = 1; i <= n; i++) {
		lua_rawgeti(L, -1, i);
		if (lua_type(L, -1) == LUA_TTABLE) {
			char cashflow_type[80] = {0};
			if (table_get_string(L, -1, "type", cashflow_type, sizeof cashflow_type)) {
				if (strcmp(cashflow_type, "simple") == 0) {
					auto cashflow = stream->add_cashflows();
					auto cf = new CFSimple();
					cashflow->set_allocated_simple(cf);
					if (!build_cashflow_simple(L, cashflows, stream, cf))
						result = false;
					else
						count++;
				} else if (strcmp(cashflow_type, "ois") == 0) {
					auto cashflow = stream->add_cashflows();
					auto cf = new CFOis();
					cashflow->set_allocated_ois(cf);
					if (!build_cashflow_ois(L, cashflows, stream, cf))
						result = false;
					else
						count++;
				} else if (strcmp(cashflow_type, "floating") == 0) {
					auto cashflow = stream->add_cashflows();
					auto cf = new CFFloating();
					cashflow->set_allocated_floating(cf);
					if (!build_cashflow_floating(L, cashflows, stream, cf))
						result = false;
					else
						count++;
				} else if (strcmp(cashflow_type, "fra") == 0) {
					auto cashflow = stream->add_cashflows();
					auto cf = new CFFra();
					auto fp = new CFFloatingPeriod();
					cf->set_allocated_floating_period(fp);
					cashflow->set_allocated_fra(cf);
					if (!build_cashflow_fra(L, cashflows, stream, cf, fp))
						result = false;
					else
						count++;
				}
			}
		} else {
			result = false;
		}
		lua_pop(L, 1);
	}
	assert(top == lua_gettop(L));
	lua_Number factor = 1.0;
	if (table_get_number(L, -1, "factor", &factor, 1.0))
		stream->set_factor(factor);
	return result && count > 0;
}

// streams table on top
// can't throw errors
static bool build_cashflow_streams(lua_State *L, CFCollection *cashflows)
{
	bool result = true;
	int top = lua_gettop(L);
	size_t n = lua_rawlen(L, -1);
	for (int i = 1; i <= n; i++) {
		lua_rawgeti(L, -1, i);
		if (lua_type(L, -1) == LUA_TTABLE) {
			auto stream = cashflows->add_streams();
			if (!build_cashflow_stream(L, cashflows, stream))
				result = false;
		}
		lua_pop(L, 1);
	}
	assert(top == lua_gettop(L));
	return result;
}

static int build_cfcollection(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1, "table expected");
	std::unique_ptr<CFCollection> cashflows = std::unique_ptr<CFCollection>(new CFCollection());
	if (build_cashflow_streams(L, cashflows.get())) {
		CFCollectionHolder *holder = (CFCollectionHolder *)lua_newuserdata(L, sizeof(CFCollectionHolder));
		luaL_getmetatable(L, Type_CFCollection);
		lua_setmetatable(L, -2);
		holder->cashflows = cashflows.release();
		return 1;
	} else {
		luaL_error(L, "failed to generate cashflows due to one or more errors");
	}
	return 0;
}

/* __gc function for Type_CFCollection */
static int collect_Type_CFCollection(lua_State *L)
{
	CFCollectionHolder *holder = check_Type_CFCollection(L, 1);
	if (holder->cashflows)
		delete holder->cashflows;
	holder->cashflows = nullptr;
	return 0;
}

static int cfcollection_tostring(lua_State *L)
{
	CFCollectionHolder *holder = check_Type_CFCollection(L, 1);
	if (holder->cashflows != nullptr) {
		auto string = holder->cashflows->ShortDebugString();
		lua_pushstring(L, string.c_str());
	} else {
		lua_pushstring(L, "<deleted>");
	}
	return 1;
}

static int get_index(lua_State *L)
{
	const char *isda_index_str = nullptr;
	const char *currency_str = nullptr;
	const char *index_str = nullptr;
	const char *tenor_str = nullptr;
	if (lua_gettop(L) == 3) {
		currency_str = luaL_checkstring(L, 1);
		index_str = luaL_checkstring(L, 2);
		tenor_str = luaL_checkstring(L, 3);
	} else {
		isda_index_str = luaL_checkstring(L, 1);
		tenor_str = luaL_checkstring(L, 2);
	}
	Currency currency = Currency::CURRENCY_UNSPECIFIED;
	IndexFamily index_family = IndexFamily::INDEX_FAMILY_UNSPECIFIED;
	IsdaIndex isda_index = IsdaIndex::ISDA_INDEX_UNSPECIFIED;
	if (currency_str) {
		currency = get_default_converter()->currency_from_string(currency_str);
		luaL_argcheck(L, currency != Currency::CURRENCY_UNSPECIFIED, 1, "Invalid currency");
	}
	if (index_str) {
		index_family = get_default_converter()->index_family_from_string(index_str);
		luaL_argcheck(L, index_family != IndexFamily::INDEX_FAMILY_UNSPECIFIED, 2, "Invalid index");
	}
	if (isda_index_str) {
		isda_index = get_default_converter()->isda_index_from_string(isda_index_str);
		luaL_argcheck(L, isda_index != IsdaIndex::ISDA_INDEX_UNSPECIFIED, 1, "Invalid index");
	}
	Period p;
	bool result = get_default_converter()->period_from_string(tenor_str, &p);
	luaL_argcheck(L, result, isda_index_str ? 2 : 3, "Invalid tenor");
	Tenor t = get_default_converter()->tenor_from_period_unit_and_len(p.units(), p.length());
	luaL_argcheck(L, t != TENOR_UNSPECIFIED, 3, "Invalid tenor");
	InterestRateIndex *index = nullptr;
	if (isda_index != IsdaIndex::ISDA_INDEX_UNSPECIFIED)
		index = get_default_index_service()->get_index(isda_index, t);
	else
		index = get_default_index_service()->get_index(currency, index_family, t);
	if (index == nullptr)
		luaL_error(L, "Index not found");
	IndexHolder *holder = (IndexHolder *)lua_newuserdata(L, sizeof(IndexHolder));
	luaL_getmetatable(L, Type_Index);
	lua_setmetatable(L, -2);
	holder->index = index;
	return 1;
}

static int index_tostring(lua_State *L)
{
	IndexHolder *idx = check_Type_Index(L, 1);
	char name[80];
	auto converter = get_default_converter();
	snprintf(name, sizeof name, "%s %s-%s", Type_Index, converter->isda_index_to_string(idx->index->isda_index()),
		 converter->tenor_to_string(idx->index->tenor()).c_str());
	lua_pushstring(L, name);
	return 1;
}

static int get_day_count_fraction(lua_State *L)
{
	const char *fraction_str = luaL_checkstring(L, 1);
	DayCountFraction fraction = get_default_converter()->day_count_fraction_from_string(fraction_str);
	luaL_argcheck(L, fraction != DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED, 2, "Invalid day count fraction");
	DayFractionHolder *holder = (DayFractionHolder *)lua_newuserdata(L, sizeof(DayFractionHolder));
	luaL_getmetatable(L, Type_DayFraction);
	lua_setmetatable(L, -2);
	holder->fraction = get_day_fraction(fraction);
	return 1;
}

// calendar(calendars [, join_rule])
//   calendars (string) can contain comma separated values
//   join rule defaults to JOIN_HOLIDAYS
static int get_calendar(lua_State *L)
{
	const char *calendar_str = luaL_checkstring(L, 1);
	std::vector<BusinessCenter> centers;
	parse_calendars(calendar_str, calendar_str + strlen(calendar_str), centers);
	luaL_argcheck(L, centers.size() > 0, 1, "Invalid calendar");
	JointCalendarRule joint_calendar_rule = JointCalendarRule::JOIN_HOLIDAYS;
	if (lua_gettop(L) == 2) {
		const char *joinmethod = luaL_checkstring(L, 2);
		joint_calendar_rule = get_default_converter()->joint_calendar_rule_from_string(joinmethod);
	}
	CalendarHolder *holder = (CalendarHolder *)lua_newuserdata(L, sizeof(CalendarHolder));
	luaL_getmetatable(L, Type_Calendar);
	lua_setmetatable(L, -2);
	holder->calendar = build_calendar(get_calendar_factory(), centers, joint_calendar_rule);
	return 1;
}

// dy.index_maturity - computes a maturity date for given index, start date
static int index_maturity(lua_State *L)
{
	IndexHolder *idx = check_Type_Index(L, 1);
	int start = luaL_checkinteger(L, 2);
	luaL_argcheck(L, is_valid_date(start), 2, "invalid date");
	auto md = idx->index->maturity_date(start);
	lua_pushinteger(L, md);
	return 1;
}

// dy.index_adjust_date - takes index, date and number of days as parameters
// and returns an adjusted date using
// fixing calendar of the index
static int index_adjust_date(lua_State *L)
{
	IndexHolder *idx = check_Type_Index(L, 1);
	lua_Integer start = luaL_checkinteger(L, 2);
	lua_Integer days = luaL_checkinteger(L, 3);
	luaL_argcheck(L, is_valid_date(start), 2, "invalid date");
	auto calendar = idx->index->fixing_calendar();
	auto dt = calendar->advance(start, Period(days, PeriodUnit::DAYS), idx->index->day_convention());
	assert(is_valid_date(dt));
	lua_pushinteger(L, dt);
	return 1;
}

// dy.index_date_components - takes index, and start date as parameters
// and returns fixing date, value date and maturity date
static int index_date_components(lua_State *L)
{
	IndexHolder *idx = check_Type_Index(L, 1);
	lua_Integer start = luaL_checkinteger(L, 2);
	luaL_argcheck(L, is_valid_date(start), 2, "invalid date");
	auto fixing_dt = idx->index->fixing_date(start);
	auto value_dt = idx->index->value_date(fixing_dt);
	auto maturity_dt = idx->index->maturity_date(value_dt);
	// TODO check for errors
	lua_pushinteger(L, fixing_dt);
	lua_pushinteger(L, value_dt);
	lua_pushinteger(L, maturity_dt);
	return 3;
}

static int day_fraction_evaluate(lua_State *L)
{
	DayFractionHolder *fraction = check_Type_DayFraction(L, 1);
	lua_Integer start = luaL_checkinteger(L, 2);
	luaL_argcheck(L, is_valid_date(start), 2, "invalid date");
	lua_Integer end = luaL_checkinteger(L, 3);
	luaL_argcheck(L, is_valid_date(end), 3, "invalid date");
	if (fraction->fraction->id() == DayCountFraction::THIRTYE_360_ISDA) {
		int isfinal = (int)luaL_checkinteger(L, 4);
		// TODO check for errors
		lua_pushnumber(L, fraction->fraction->year_fraction(start, end, isfinal != 0));
	} else if (fraction->fraction->id() == DayCountFraction::ACT_ACT_ISMA) {
		lua_Integer refStart = start;
		lua_Integer refEnd = end;
		if (lua_gettop(L) == 4) {
			refStart = luaL_checkinteger(L, 4);
			luaL_argcheck(L, is_valid_date(refStart), 4, "invalid date");
		}
		if (lua_gettop(L) == 5) {
			refEnd = luaL_checkinteger(L, 5);
			luaL_argcheck(L, is_valid_date(refEnd), 5, "invalid date");
		}
		lua_pushnumber(L, fraction->fraction->year_fraction(start, end, refStart, refEnd));
	} else {
		// TODO check for errors
		lua_pushnumber(L, fraction->fraction->year_fraction(start, end));
	}
	return 1;
}

static int day_fraction_tostring(lua_State *L)
{
	DayFractionHolder *fraction = check_Type_DayFraction(L, 1);
	char name[80];
	auto converter = get_default_converter();
	snprintf(name, sizeof name, "%s %s", Type_DayFraction,
		 converter->day_count_fraction_to_string(fraction->fraction->id()));
	lua_pushstring(L, name);
	return 1;
}

static int calendar_advance(lua_State *L)
{
	CalendarHolder *calendar = check_Type_Calendar(L, 1);
	Date d = luaL_checkinteger(L, 2);
	luaL_argcheck(L, is_valid_date(d), 2, "invalid date");
	BusinessDayConvention bdc = BusinessDayConvention::FOLLOWING;
	if (lua_gettop(L) == 4 && lua_isstring(L, 4)) {
		auto temp = get_default_converter()->business_day_convention_from_string(lua_tostring(L, 4));
		luaL_argcheck(L,
			      temp != BusinessDayConvention::BUSINESS_DAY_CONVENTION_UNSPECIFIED &&
				  BusinessDayConvention_IsValid(temp),
			      4, "invalid business day convention");
	}
	if (lua_type(L, 3) == LUA_TSTRING) {
		const char *period = luaL_checkstring(L, 3);
		Period p;
		luaL_argcheck(L, get_default_converter()->period_from_string(period, &p), 3, "invalid period");
		d = calendar->calendar->advance(d, p, bdc);
	} else {
		lua_Integer n = luaL_checkinteger(L, 3);
		d = calendar->calendar->advance(d, (int)n, PeriodUnit::DAYS, bdc);
	}
	lua_pushinteger(L, d);
	return 1;
}

static int calendar_tostring(lua_State *L)
{
	CalendarHolder *calendar = check_Type_Calendar(L, 1);
	char name[80];
	std::array<BusinessCenter, 4> ids;
	calendar->calendar->get_ids(ids);
	auto converter = get_default_converter();
	name[0] = 0;
	snprintf(name, sizeof name, "%s ", Type_Calendar);
	for (int i = 0; i < ids.size(); i++) {
		if (ids[i] != BusinessCenter::BUSINESS_CENTER_UNSPECIFIED) {
			const char *s = converter->business_center_to_string(ids[i]);
			auto n = strlen(name);
			snprintf(name + n, sizeof name - n, "%s%s", (i > 0 ? ":" : ""), s);
		}
	}
	lua_pushstring(L, name);
	return 1;
}

static bool parse_future_expiry(lua_State *L, const char *future_expiry, Date *d)
{
	auto len = strlen(future_expiry);
	bool ok = false;
	if (len == 4) {
		char buf[5];
		strcpy(buf, future_expiry);
		// FIXME
		int year = 2000 + atoi(buf + 2);
		*(buf + 2) = 0;
		int month = atoi(buf);
		*d = make_date(1, month, year);
		ok = true;
	} else if (len == 5) {
		// FIXME
		int year = 2000 + atoi(future_expiry + 3);
		int month = 0;
		switch (*future_expiry) {
		case 'J':
		case 'j':
			if (future_expiry[1] == 'a' || future_expiry[1] == 'A')
				month = 1;
			else if (future_expiry[2] == 'l' || future_expiry[2] == 'L')
				month = 7;
			else
				month = 6;
			break;
		case 'F':
		case 'f':
			month = 2;
			break;
		case 'A':
		case 'a':
			if (future_expiry[1] == 'p' || future_expiry[1] == 'P')
				month = 4;
			else
				month = 8;
			break;
		case 'm':
		case 'M':
			if (future_expiry[2] == 'y' || future_expiry[2] == 'Y')
				month = 5;
			else
				month = 3;
			break;
		case 'S':
		case 's':
			month = 9;
			break;
		case 'O':
		case 'o':
			month = 10;
			break;
		case 'N':
		case 'n':
			month = 11;
			break;
		case 'D':
		case 'd':
			month = 12;
			break;
		default:
			ravi_writestringerror(L, "bad month in future expiry '%s'\n", future_expiry);
			return false;
		}
		*d = make_date(1, month, year);
		ok = true;
	} else {
		ravi_writestringerror(L, "failed to parse '%s'\n", future_expiry);
	}
	return ok;
}

// dy.future_imm_dates - takes a currency and future expiry (mmyy or MMMyy)
// and returns that start and end IMM dates of the future instrument
static int future_imm_dates(lua_State *L)
{
	const char *currency = luaL_checkstring(L, -2);
	const char *future_expiry = luaL_checkstring(L, -1);

	RollConvention rc = RollConvention::IMM;
	if (strcmp(currency, "AUD") == 0)
		rc = RollConvention::IMMAUD;
	else if (strcmp(currency, "CAD") == 0)
		rc = RollConvention::IMMCAD;
	else if (strcmp(currency, "NZD") == 0)
		rc = RollConvention::IMMNZD;

	Date expiry;
	if (!parse_future_expiry(L, future_expiry, &expiry))
		luaL_error(L,
			   "Invalid argument, failed to determine imm future "
			   "date for %s",
			   future_expiry);
	else {
		Date start = sub(expiry, Period(3, PeriodUnit::MONTHS));
		expiry = adjust_date(expiry, rc, 0);
		start = adjust_date(start, rc, 0);

		lua_pushinteger(L, start);
		lua_pushinteger(L, expiry);
		return 2;
	}
	return 0;
}

// dy.fra_tenors - takes a string of format nnXnn and returns two numbers
static int fra_tenors(lua_State *L)
{
	const char *fra = luaL_checkstring(L, -1);

	const char *cp = fra;
	while (*cp && !isdigit(*cp)) {
		cp++;
	}
	int p1 = atoi(cp);
	while (*cp && !(*cp == 'X' || *cp == 'x'))
		cp++;
	if (*cp == 'X' || *cp == 'x')
		cp++;
	int p2 = atoi(cp);

	lua_pushinteger(L, p1);
	lua_pushinteger(L, p2);

	return 2;
}

// Reads a CSV into a table
// Columns are converted as specified
// Supported conversions are string, number, integer, date
static int read_csv(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1, "table expected with attributes 'file', 'conversion'");
	char filename[1024] = {0};
	char instructions[100] = {0};

	luaL_argcheck(L, table_get_string(L, 1, "file", filename, sizeof filename), 1, "file attribute expected");
	table_get_string(L, 1, "conversion", instructions, sizeof instructions);

	bool has_headings = false;
	table_get_boolean(L, 1, "heading", &has_headings);

	lua_Integer missing_integer = 0;
	table_get_lua_Integer(L, 1, "default_integer", &missing_integer);
	lua_Number missing_number = 0;
	table_get_number(L, 1, "default_number", &missing_number, 0.0);
	Date missing_date = 0;
	table_get_integer(L, 1, "default_date", &missing_date);
	bool use_headings = false;
	if (has_headings)
		table_get_boolean(L, 1, "fields", &use_headings);

	int ncols = strlen(instructions);

	lua_newtable(L);
	int lno = 0;
	CSVDataSource ds(std::string(filename), has_headings);
	if (!ds.is_valid()) {
		return 0;
	}
	auto line = std::make_unique<Line>(1024);
	while (ds.next(line.get())) {
		Date d = 0;
		lua_newtable(L);
		for (int f = 0; f < line->fields.size(); f++) {
			if (use_headings && f >= line->headings.size())
				break;
			char t = f < ncols ? instructions[f] : 's';
			switch (t) {
			case 'n':
				if (line->fields[f] && line->fields[f][0])
					lua_pushnumber(L, strtod(line->fields[f], nullptr));
				else
					lua_pushnumber(L, missing_number);
				break;
			case 'i':
				if (line->fields[f] && line->fields[f][0])
					lua_pushinteger(L, atoi(line->fields[f]));
				else
					lua_pushinteger(L, missing_integer);
				break;
			case 'd':
				if (parse_date(line->fields[f], &d))
					lua_pushinteger(L, d);
				else
					lua_pushinteger(L, missing_date);
				break;
			case '-':
				lua_pushnil(L);
				break;
			case 's':
			default:
				lua_pushstring(L, line->fields[f]);
				break;
			}
			if (use_headings) {
				const char *h = line->headings[f];
				lua_setfield(L, -2, h);
			} else
				lua_rawseti(L, -2, f + 1);
		}
		lua_rawseti(L, -2, lno + 1);
		lno++;
	}
	return 1;
}

// Reads a matrix from a CSV or tab delimited file into a table
// All rows and columns are assumed to contain numeric data
static int read_matrix(lua_State *L)
{
	luaL_argcheck(L, lua_isstring(L, 1), 1, "filename expected");
	char filename[1024] = {0};
	snprintf(filename, sizeof filename, "%s", lua_tostring(L, 1));
	lua_newtable(L);
	int lno = 0;
	CSVDataSource ds(std::string(filename), false);
	auto line = std::make_unique<Line>(64 * 1024);
	while (ds.next(line.get())) {
		ravi_create_number_array(L, line->fields.size(), 0.0);
		for (int f = 0; f < line->fields.size(); f++) {
			lua_pushnumber(L, strtod(line->fields[f], nullptr));
			lua_rawseti(L, -2, f + 1);
		}
		lua_rawseti(L, -2, lno + 1);
		lno++;
	}
	return 1;
}

static void replace_all(std::string &str, const std::string &text, const std::string &value)
{
	if (text == value)
		return;
	for (size_t pos = 0; pos != std::string::npos;) {
		pos = str.find(text, pos);
		if (pos == std::string::npos)
			break;
		str.replace(pos, text.length(), value);
	}
}

/**
strrepl(str, table)
For each key appearing in table, this function searches
str replaces all occurences with the associated value.
Faster than doing this in Lua
*/
static int string_replace(lua_State *L)
{
	const char *s = luaL_checkstring(L, 1);
	luaL_argcheck(L, lua_istable(L, 2), 2, "table expected");
	std::string str(s);
	lua_pushnil(L); // push first key
	while (lua_next(L, 2)) {
		if (lua_isstring(L, -2) && lua_isstring(L, -1)) {
			std::string k(lua_tostring(L, -2));
			std::string v(lua_tostring(L, -1));
			replace_all(str, k, v);
		}
		lua_pop(L, 1);
	}
	lua_pushstring(L, str.c_str());
	return 1;
}

static int stats_mean(lua_State *L)
{
	luaL_argcheck(L, ravi_is_number_array(L, 1), 1, "number[] expected");
	size_t len = 0;
	double *vec = ravi_get_number_array_rawdata(L, 1, &len);
	vec++, len--;
	double m = mean(vec, vec + len);
	lua_pushnumber(L, m);
	return 1;
}

static int stats_variance(lua_State *L)
{
	luaL_argcheck(L, ravi_is_number_array(L, 1), 1, "number[] expected");
	size_t len = 0;
	double *vec = ravi_get_number_array_rawdata(L, 1, &len);
	vec++, len--;
	double m = variance(vec, vec + len);
	lua_pushnumber(L, m);
	return 1;
}

static int stats_stddev(lua_State *L)
{
	luaL_argcheck(L, ravi_is_number_array(L, 1), 1, "number[] expected");
	size_t len = 0;
	double *vec = ravi_get_number_array_rawdata(L, 1, &len);
	vec++, len--;
	double m = standard_deviation(vec, vec + len);
	lua_pushnumber(L, m);
	return 1;
}

static int stats_lagged_difference(lua_State *L)
{
	luaL_argcheck(L, ravi_is_number_array(L, 1), 1, "number[] expected");
	luaL_argcheck(L, lua_isinteger(L, 2), 2, "integer lag expected");
	DiffType diff_type = DiffType::DIFFERENCE_ABSOLUTE;
	if (lua_gettop(L) == 3) {
		diff_type = DiffType::DIFFERENCE_RELATIVE;
	}
	int lag = lua_tointeger(L, 2);
	size_t len = 0;
	double *vec = ravi_get_number_array_rawdata(L, 1, &len);
	vec++, len--;
	luaL_argcheck(L, lag >= 0 && lag < len, 2, "lag must be >=0 and < length of number array");
	ravi_create_number_array(L, len - lag, 0.0);
	size_t len2 = 0;
	double *vec2 = ravi_get_number_array_rawdata(L, -1, &len2);
	vec2++;
	lagged_difference(vec2, vec, len, lag, diff_type);
	return 1;
}

static int stats_volatility(lua_State *L)
{
	luaL_argcheck(L, ravi_is_number_array(L, 1), 1, "number[] expected");
	luaL_argcheck(L, lua_isnumber(L, 2), 2, "number parameter (lambda) expected");
	luaL_argcheck(L, lua_isnumber(L, 3), 3, "number parameter (volatility seed) expected");
	VolatilityEstimator disp = VolatilityEstimator::STANDARD_DEVIATION;
	double lambda = lua_tonumber(L, 2);
	luaL_argcheck(L, lambda >= 0 && lambda <= 1.0, 2, "lambda must be >=0 and <= 1.0");
	double seed = lua_tonumber(L, 3);
	size_t len = 0;
	double *vec = ravi_get_number_array_rawdata(L, 1, &len);
	vec++, len--;
	ravi_create_number_array(L, len, 0.0);
	size_t len2 = 0;
	double *vec2 = ravi_get_number_array_rawdata(L, -1, &len2);
	vec2++;
	ewma_volatility(vec2, vec, len, lambda, seed, disp);
	return 1;
}

static struct redukti_adouble_t *alloc_adouble(lua_State *L, int nvars, int order, int var, double v)
{
	auto adsize = redukti_adouble_alloc_size(nvars, order);
	redukti_adouble_t *adnum = (redukti_adouble_t *)lua_newuserdata(L, adsize);
	redukti_adouble_init(adnum, nvars, order, var, v);
	luaL_getmetatable(L, Type_ADouble);
	lua_setmetatable(L, -2);
	return adnum;
}

static int create_adouble_internal(lua_State *L, int order)
{
	int nvars = 0;
	bool table_input = false;
	if (lua_gettop(L) == 1 && lua_istable(L, 1)) {
		nvars = lua_rawlen(L, 1);
		table_input = true;
	} else {
		nvars = lua_gettop(L);
	}
	luaL_argcheck(L, nvars > 0, 1, "missing argument");
	luaL_argcheck(L, nvars <= 100, 1, "number of variables restricted to 100");
	int n = 0;
	for (int i = 0; i < nvars; i++) {
		double v = 0.0;
		if (table_input) {
			lua_rawgeti(L, 1, i + 1);
			v = lua_tonumber(L, -1);
			lua_pop(L, 1);
		} else {
			v = lua_tonumber(L, i + 1);
		}
		alloc_adouble(L, nvars, order, i, v);
		n++;
	}
	return n;
}

// Create an adouble of order 1
static int create_adouble_order1(lua_State *L) { return create_adouble_internal(L, 1); }

// Create an adouble of order 2
static int create_adouble_order2(lua_State *L) { return create_adouble_internal(L, 2); }

static int adouble_tostring(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	luaL_Buffer b;
	luaL_buffinit(L, &b);
	luaL_addstring(&b, "{\n");
	luaL_addstring(&b, "  value=");
	lua_pushnumber(L, redukti_adouble_get_value(ad));
	luaL_addvalue(&b);
	if (ad->order_ > 0)
		luaL_addstring(&b, ",");
	luaL_addstring(&b, "\n");
	if (ad->order_ > 0) {
		luaL_addstring(&b, "  firstorder = {\n");
		int count = 0;
		for (int i = 0; i < ad->vars_; i++) {
			double v = redukti_adouble_get_derivative1(ad, i);
			if (v != 0.0) {
				if (count) {
					luaL_addstring(&b, ",\n");
				}
				luaL_addstring(&b, "    [");
				lua_pushinteger(L, i + 1);
				luaL_addvalue(&b);
				luaL_addstring(&b, "] = ");
				lua_pushnumber(L, v);
				luaL_addvalue(&b);
				count++;
			}
		}
		luaL_addstring(&b, "\n  }\n");
		if (ad->order_ == 2) {
			count = 0;
			luaL_addstring(&b, "  secondorder = {\n");
			for (int i = 0; i < ad->vars_; i++) {
				for (int j = 0; j < ad->vars_; j++) {
					double v = redukti_adouble_get_derivative2(ad, i, j);
					if (v != 0.0) {
						if (count) {
							luaL_addstring(&b, ",\n");
						}
						luaL_addstring(&b, "    [");
						lua_pushinteger(L, i + 1);
						luaL_addvalue(&b);
						luaL_addstring(&b, ", ");
						lua_pushinteger(L, j + 1);
						luaL_addvalue(&b);
						luaL_addstring(&b, "] = ");
						lua_pushnumber(L, v);
						luaL_addvalue(&b);
						count++;
					}
				}
			}
			if (count)
				luaL_addstring(&b, "\n");
			luaL_addstring(&b, "  }\n");
		}
	}
	luaL_addstring(&b, "}");
	luaL_pushresult(&b);
	return 1;
}

static int adouble_get(lua_State *L)
{
	int nargs = lua_gettop(L);
	if (nargs == 0) {
		return 0;
	}
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	switch (nargs) {
	case 1: {
		lua_pushnumber(L, redukti_adouble_get_value(ad));
		return 1;
	}
	case 2: {
		int first_order_deriv = luaL_checkinteger(L, 2);
		luaL_argcheck(L, ad->order_ >= 1 && first_order_deriv >= 1 && first_order_deriv <= ad->vars_, 2,
			      "variable reference out of range");
		lua_pushnumber(L, redukti_adouble_get_derivative1(ad, first_order_deriv - 1));
		return 1;
	}
	case 3: {
		int second_order_deriv1 = luaL_checkinteger(L, 2);
		int second_order_deriv2 = luaL_checkinteger(L, 3);
		luaL_argcheck(L, ad->order_ == 2 && second_order_deriv1 >= 1 && second_order_deriv1 <= ad->vars_, 2,
			      "variable reference out of range");
		luaL_argcheck(L, ad->order_ == 2 && second_order_deriv2 >= 1 && second_order_deriv2 <= ad->vars_, 2,
			      "variable reference out of range");
		lua_pushnumber(L,
			       redukti_adouble_get_derivative2(ad, second_order_deriv1 - 1, second_order_deriv2 - 1));
		return 1;
	}
	}
	luaL_argerror(L, 1, "Invalid arguments");
	return 0;
}

static int adouble_get_gradient(lua_State *L)
{
	int nargs = lua_gettop(L);
	if (nargs == 0) {
		return 0;
	}
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	if (ad->order_ < 1) {
		return 0;
	}
	ravi_create_number_array(L, ad->vars_, 0.0);
	size_t len = 0;
	double *v = ravi_get_number_array_rawdata(L, -1, &len);
	for (int i = 1; i <= ad->vars_; i++)
		v[i] = redukti_adouble_get_derivative1(ad, i - 1);
	return 1;
}

static int adouble_get_hessian(lua_State *L)
{
	int nargs = lua_gettop(L);
	if (nargs == 0) {
		return 0;
	}
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	if (ad->order_ < 2) {
		return 0;
	}
	lua_createtable(L, ad->vars_, 0);
	for (int j = 1; j <= ad->vars_; j++) {
		ravi_create_number_array(L, ad->vars_, 0.0);
		size_t len = 0;
		double *v = ravi_get_number_array_rawdata(L, -1, &len);
		for (int i = 1; i <= ad->vars_; i++)
			v[i] = redukti_adouble_get_derivative2(ad, j - 1, i - 1);
		lua_rawseti(L, -2, j);
	}
	return 1;
}

static int adouble_add(lua_State *L)
{
	if (lua_isnumber(L, 1)) {
		lua_Number x = lua_tonumber(L, 1);
		redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, x);
		redukti_adouble_add(result, ad, 1.0);
		return 1;
	}
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	if (lua_isnumber(L, 2)) {
		lua_Number x = lua_tonumber(L, 2);
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, x);
		redukti_adouble_add(result, ad, 1.0);
		return 1;
	}
	redukti_adouble_t *ad2 = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
	luaL_argcheck(L, ad->vars_ == ad2->vars_, 2, "vars mismatch");
	luaL_argcheck(L, ad->order_ == ad2->order_, 2, "order mismatch");
	redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
	redukti_adouble_assign(result, ad);
	redukti_adouble_add(result, ad2, 1.0);
	return 1;
}

static int adouble_sub(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	if (lua_isnumber(L, 2)) {
		lua_Number x = lua_tonumber(L, 2);
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		result->data_[0] -= x;
		return 1;
	}
	redukti_adouble_t *ad2 = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
	luaL_argcheck(L, ad->vars_ == ad2->vars_, 2, "vars mismatch");
	luaL_argcheck(L, ad->order_ == ad2->order_, 2, "order mismatch");
	redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
	redukti_adouble_assign(result, ad);
	redukti_adouble_add(result, ad2, -1.0);
	return 1;
}

static int adouble_mul(lua_State *L)
{
	if (lua_isnumber(L, 1)) {
		// number * adouble
		lua_Number x = lua_tonumber(L, 1);
		redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		redukti_adouble_scalar_multiply(result, x);
		return 1;
	}
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	if (lua_isnumber(L, 2)) {
		// adouble * number
		lua_Number x = lua_tonumber(L, 2);
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		redukti_adouble_scalar_multiply(result, x);
		return 1;
	}
	// adouble * adouble
	redukti_adouble_t *ad2 = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
	luaL_argcheck(L, ad->vars_ == ad2->vars_, 2, "vars mismatch");
	luaL_argcheck(L, ad->order_ == ad2->order_, 2, "order mismatch");
	redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
	redukti_adouble_assign(result, ad);
	auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
	redukti_adouble_t *tempbuf = (redukti_adouble_t *)alloca(adsize);
	redukti_adouble_multiply(result, ad2, tempbuf);
	return 1;
}

static int adouble_div(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		if (lua_isnumber(L, 2)) {
			// adouble * 1.0 / n
			lua_Number x = lua_tonumber(L, 2);
			luaL_argcheck(L, x != 0.0, 2, "divide by zero");
			redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
			redukti_adouble_assign(result, ad);
			redukti_adouble_scalar_multiply(result, 1.0 / x);
			return 1;
		}
		// adouble / adouble
		redukti_adouble_t *ad2 = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
		luaL_argcheck(L, ad->vars_ == ad2->vars_, 2, "vars mismatch");
		luaL_argcheck(L, ad->order_ == ad2->order_, 2, "order mismatch");
		luaL_argcheck(L, ad2->data_[0] != 0.0, 2, "divide by zero");
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_t *tempbuf2 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_divide(result, ad2, tempbuf1, tempbuf2);
		return 1;
	} else {
		// n * 1.0 / adouble
		lua_Number x = luaL_checknumber(L, 1);
		ad = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_power(result, -1.0, tempbuf1);
		redukti_adouble_scalar_multiply(result, x);
		return 1;
	}
}

static int adouble_exp(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_exp(result, tempbuf1);
	} else {
		lua_pushnumber(L, std::exp(luaL_checknumber(L, 1)));
	}
	return 1;
}

static int adouble_power(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	lua_Number n = luaL_checknumber(L, 2);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_power(result, n, tempbuf1);
	} else {
		lua_pushnumber(L, std::pow(luaL_checknumber(L, 1), n));
	}
	return 1;
}

static int adouble_sqrt(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_power(result, 0.5, tempbuf1);
	} else {
		lua_pushnumber(L, std::sqrt(luaL_checknumber(L, 1)));
	}
	return 1;
}

static int adouble_log(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_log(result, tempbuf1);
	} else {
		lua_pushnumber(L, std::log(luaL_checknumber(L, 1)));
	}
	return 1;
}

static int adouble_sin(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_sin(result, tempbuf1);
	} else {
		lua_pushnumber(L, std::sin(luaL_checknumber(L, 1)));
	}
	return 1;
}

static int adouble_cos(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_cos(result, tempbuf1);
	} else {
		lua_pushnumber(L, std::cos(luaL_checknumber(L, 1)));
	}
	return 1;
}

static int adouble_tan(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		auto adsize = redukti_adouble_alloc_size(ad->vars_, ad->order_);
		redukti_adouble_t *tempbuf1 = (redukti_adouble_t *)alloca(adsize);
		redukti_adouble_tan(result, tempbuf1);
	} else {
		lua_pushnumber(L, std::tan(luaL_checknumber(L, 1)));
	}
	return 1;
}

// absolute
static int adouble_abs(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad) {
		redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
		redukti_adouble_assign(result, ad);
		if (result->data_[0] < 0.0) {
			redukti_adouble_scalar_multiply(result, -1.0);
		}
	} else {
		lua_pushnumber(L, std::fabs(luaL_checknumber(L, 1)));
	}
	return 1;
}

// unary minus
static int adouble_unm(lua_State *L)
{
	redukti_adouble_t *ad = (redukti_adouble_t *)luaL_checkudata(L, 1, Type_ADouble);
	redukti_adouble_t *result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
	redukti_adouble_assign(result, ad);
	redukti_adouble_scalar_multiply(result, -1.0);
	return 1;
}

// Compare two adouble objects and return
// 0.0 if equal
// < 0.0 if a < b
// > 0.0 if a > b
static double adouble_compare(lua_State *L)
{
	redukti_adouble_t *ad1 = (redukti_adouble_t *)luaL_testudata(L, 1, Type_ADouble);
	if (ad1) {
		if (lua_isnumber(L, 2)) {
			lua_Number x = lua_tonumber(L, 2);
			return ad1->data_[0] - x;
		} else {
			redukti_adouble_t *ad2 = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
			luaL_argcheck(L, ad1->vars_ == ad2->vars_, 2, "vars mismatch");
			luaL_argcheck(L, ad1->order_ == ad2->order_, 2, "order mismatch");
			return ad1->data_[0] - ad2->data_[0];
		}
	} else if (lua_isnumber(L, 1)) {
		lua_Number x = lua_tonumber(L, 2);
		redukti_adouble_t *ad2 = (redukti_adouble_t *)luaL_checkudata(L, 2, Type_ADouble);
		return x - ad2->data_[0];
	} else {
		luaL_argerror(L, 1, "wrong type");
		return -1;
	}
}

static int adouble_eq(lua_State *L)
{
	lua_pushboolean(L, adouble_compare(L) == 0.0);
	return 1;
}

static int adouble_lt(lua_State *L)
{
	lua_pushboolean(L, adouble_compare(L) < 0.0);
	return 1;
}

static int adouble_le(lua_State *L)
{
	lua_pushboolean(L, adouble_compare(L) <= 0.0);
	return 1;
}

/* __gc function for Type_Interpolator */
static int collect_Type_Interpolator(lua_State *L)
{
	InterpolatorHolder *holder = check_Type_Interpolator(L, 1);
	holder->~InterpolatorHolder();
	return 0;
}

static int create_interpolator(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1,
		      "table expected"); // will do longjmp

	// We have to be very careful below in our error handling
	// as Lua error handling causes a longjmp which will lead to
	// leaks with C++ objects
	// This code is structured so that we do not
	// raise Lua errors while we have C++ objects in scope

	InterpolatorType it = InterpolatorType::LINEAR; // default
	int order = 0;
	InterpolationOptions options;
	char str[80] = {0};

	if (table_get_integer(L, 1, "order", &order) && order >= 1 && order <= 2) {
		options.differentiation_order = order;
	}
	if (table_get_string(L, 1, "interpolator", str, sizeof str)) {
		it = get_default_converter()->interpolator_type_from_string(str);
	}

	size_t xlen = table_get_arraylength(L, 1, "x"), ylen = table_get_arraylength(L, 1, "y");
	const char *errmsg = nullptr;
	if (xlen < 2) {
		errmsg = "at least two elements must be present in x and y values";
	} else if (xlen != ylen) {
		errmsg = "number of elements in x and y must agree";
	} else if (xlen > 50 && options.differentiation_order > 1) {
		errmsg = "second order derivatives limit number of elements to 50";
	} else if (xlen > 100) {
		errmsg = "length of x and y must be <= 100";
	}
	while (!errmsg) {
		std::unique_ptr<double> xvalues = std::unique_ptr<double>(new double[xlen]);
		std::unique_ptr<double> yvalues = std::unique_ptr<double>(new double[ylen]);
		auto xn = table_get_doublearray(L, 1, "x", xvalues.get(), xlen);
		auto yn = table_get_doublearray(L, 1, "y", yvalues.get(), ylen);

		if (xn != xlen || yn != ylen) {
			errmsg = "failed to get the number of values into array as "
				 "expected";
			break;
		}
		double *x = xvalues.get();
		for (int i = 1; i < xlen; i++) {
			if (x[i] <= x[i - 1]) {
				errmsg = "x values must be sorted in ascending "
					 "order";
				break;
			}
		}
		if (errmsg)
			break;
		std::unique_ptr<Interpolator, Deleter<Interpolator>> interp =
		    make_interpolator(it, xvalues.get(), yvalues.get(), xlen, &GlobalAllocator, options);
		InterpolatorHolder *iptr = (InterpolatorHolder *)lua_newuserdata(L, sizeof(InterpolatorHolder));
		new (iptr) InterpolatorHolder(std::move(xvalues), std::move(yvalues), std::move(interp));
		luaL_getmetatable(L, Type_Interpolator);
		lua_setmetatable(L, -2);
		break;
	}
	if (errmsg) {
		luaL_argerror(L, 1, errmsg);
	}
	return errmsg ? 0 : 1;
}

static int interpolator_get(lua_State *L)
{
	InterpolatorHolder *iptr = check_Type_Interpolator(L, 1);
	lua_Number x = luaL_checknumber(L, 2);
	int rc = 0;
	if (lua_gettop(L) == 2) {
		lua_pushnumber(L, iptr->interp->interpolate(x));
		rc = 1;
	} else {
		int i = (int)lua_tointeger(L, 2);

		auto allocSet = get_threadspecific_allocators();
		auto A = allocSet->tempspace_allocator;
		FixedRegionAllocatorGuard g(A);
		if (i == 1) {
			auto ad = iptr->interp->interpolate_with_sensitivities(x, A);
			if (ad) {
				auto result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
				redukti_adouble_assign(result, ad.get());
				rc = 1;
			}
		} else {
			auto ad = iptr->interp->interpolate_with_numeric_sensitivities(x, A);
			if (ad) {
				auto result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
				redukti_adouble_assign(result, ad.get());
				rc = 1;
			}
		}
	}
	return rc;
}

/* __gc for Type_ValuationContext */
static int collect_Type_ValuationContext(lua_State *L)
{
	ValuationContextHolder *holder = check_Type_ValuationContext(L, 1);
	holder->~ValuationContextHolder();
	return 0;
}

static int make_valuation_context(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1, "table expected");
	FixingServiceHolder *fixing_service = test_Type_FixingService(L, 2);

	int as_of_date = 0;
	int temp = 0;
	int order = 0;
	if (table_get_integer(L, 1, "reference_date", &temp) && is_valid_date(temp)) {
		as_of_date = temp;
	}
	if (table_get_integer(L, 1, "order", &temp) && temp >= 1 && temp <= 2) {
		order = temp;
	}
	luaL_argcheck(L, is_valid_date(as_of_date), 1, "invalid reference date");
	ValuationContextHolder *iptr = (ValuationContextHolder *)lua_newuserdata(L, sizeof(ValuationContextHolder));
	new (iptr)
	    ValuationContextHolder(as_of_date, order, fixing_service ? &fixing_service->fixing_service : nullptr);
	iptr->context.set_derivative_order(order);
	luaL_getmetatable(L, Type_ValuationContext);
	lua_setmetatable(L, -2);
	if (fixing_service) {
		// We need to keep a reference to the fixing service in our uservalue
		// so that it doesn't get garbage collected while this valuation context
		// is alive
		lua_newtable(L);
		assert(test_Type_FixingService(L, -3));
		lua_pushvalue(L, -3);
		lua_setfield(L, -2, "fixing_service");
		assert(test_Type_ValuationContext(L, -2));
		lua_setuservalue(L, -2);
	}
	return 1;
}

/* __gc for Type_CurveMapper */
static int collect_Type_CurveMapper(lua_State *L)
{
	CurveMapperHolder *holder = check_Type_CurveMapper(L, 1);
	holder->~CurveMapperHolder();
	return 0;
}

static int make_curve_mapper(lua_State *L)
{
	CurveMapperHolder *iptr = (CurveMapperHolder *)lua_newuserdata(L, sizeof(CurveMapperHolder));
	new (iptr) CurveMapperHolder();
	luaL_getmetatable(L, Type_CurveMapper);
	lua_setmetatable(L, -2);
	return 1;
}

static int add_curve_mapping(lua_State *L)
{
	CurveMapperHolder *holder = check_Type_CurveMapper(L, 1);
	lua_Integer from_id = luaL_checkinteger(L, 2);
	lua_Integer to_id = luaL_checkinteger(L, 3);

	PricingCurve c1 = PricingCurve((uint32_t)from_id);
	PricingCurve c2 = PricingCurve((uint32_t)to_id);

	luaL_argcheck(L, c1.is_valid(), 2, "invalid pricing curve id");
	luaL_argcheck(L, c2.is_valid(), 3, "invalid pricing curve id");

	holder->mapper.add_mapping(c1, c2);
	return 0;
}

static int get_curve_mapping(lua_State *L)
{
	CurveMapperHolder *holder = check_Type_CurveMapper(L, 1);
	lua_Integer from_id = luaL_checkinteger(L, 2);

	PricingCurve c1 = PricingCurve((uint32_t)from_id);

	luaL_argcheck(L, c1.is_valid(), 2, "invalid pricing curve id");

	auto mapped = holder->mapper.map_index_tenor(c1.curve_type(), c1.currency(), c1.index_family(), c1.tenor());

	lua_pushinteger(L, mapped.id());
	return 1;
}

/* __gc for Type_CurveMapper */
static int collect_Type_CurveProvider(lua_State *L)
{
	CurveProviderHolder *holder = check_Type_CurveProvider(L, 1);
	holder->~CurveProviderHolder();
	return 0;
}

static int make_curve_provider(lua_State *L)
{
	CurveProviderHolder *iptr = (CurveProviderHolder *)lua_newuserdata(L, sizeof(CurveProviderHolder));
	new (iptr) CurveProviderHolder();
	luaL_getmetatable(L, Type_CurveProvider);
	lua_setmetatable(L, -2);
	// All curves associated with the provider
	// are referenced by a table that we hold as uservalue.
	// This is to avoid the curves being garbage collected
	lua_newtable(L);
	lua_setuservalue(L, -2);
	return 1;
}

static int add_curve_provider_mapping(lua_State *L)
{
	CurveProviderHolder *holder = check_Type_CurveProvider(L, 1);
	lua_Integer curve_id = luaL_checkinteger(L, 2);
	IRCurveHolder *curve = check_Type_IRCurve(L, 3);
	int top = lua_gettop(L);
	PricingCurve c1 = PricingCurve((uint32_t)curve_id);

	luaL_argcheck(L, c1.is_valid(), 2, "invalid pricing curve id");

	holder->provider.add_mapping(c1, curve);

	// We want the provider to hold a reference to the curve
	// We do this by putting all curves in a table and
	// associate this with uservalue of the provider
	// Furthermore we index the curves by their curve id
	// so that by simply returning the table, the Lua
	// client can lookup the curves!
	int type = lua_getuservalue(L, -3); // index 1
	lua_assert(type == LUA_TTABLE);
	// table at -1, curve at -2
	lua_pushinteger(L, curve_id);
	// curve_id at -1, table at -2, curve at -3
	lua_pushvalue(L, -3);
	// curve at -1, curve_id at -2, table at -3
	lua_rawset(L, -3);
	// table at -1
	// pop the table
	lua_pop(L, 1);

	lua_assert(top == lua_gettop(L));
	return 0;
}

static int curve_provider_get_curves(lua_State *L)
{
	CurveProviderHolder *holder = check_Type_CurveProvider(L, 1);
	// FIXME we should return a copy as else user may
	// remove some curves causing a potential memory
	// management issue
	lua_getuservalue(L, 1);
	return 1;
}

/* __gc function for Type_IRCurve */
static int collect_Type_IRCurve(lua_State *L)
{
	IRCurveHolder *holder = check_Type_IRCurve(L, 1);
	holder->~IRCurveHolder();
	return 0;
}

static int create_ircurve(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1,
		      "table expected"); // will do longjmp

	// We have to be very careful below in our error handling
	// as Lua error handling causes a longjmp which will lead to
	// leaks with C++ objects

	char str[80] = {0};
	int temp = 0;
	int order = 0;
	InterpolatorType it = InterpolatorType::LINEAR; // default
	Currency ccy = Currency::CURRENCY_UNSPECIFIED;
	IndexFamily indexFamily = IndexFamily::INDEX_FAMILY_UNSPECIFIED;
	Tenor tenor = Tenor::TENOR_UNSPECIFIED;
	MarketDataQualifier qual = MarketDataQualifier::MDQ_NORMAL;
	IRRateType rate_type = IRRateType::ZERO_RATE;
	DayCountFraction dfc = DayCountFraction::ACT_365_FIXED;
	PricingCurveType curve_type = PricingCurveType::PRICING_CURVE_TYPE_FORWARD;
	auto converter = get_default_converter();
	int as_of_date = 0;

	if (table_get_integer(L, 1, "order", &temp) && temp >= 1 && temp <= 2) {
		order = temp;
	}
	if (table_get_integer(L, 1, "reference_date", &temp) && is_valid_date(temp)) {
		as_of_date = temp;
	}
	if (table_get_string(L, 1, "interpolator", str, sizeof str)) {
		it = converter->interpolator_type_from_string(str);
	}
	if (table_get_string(L, 1, "currency", str, sizeof str)) {
		ccy = converter->currency_from_string(str);
	}
	if (table_get_string(L, 1, "index_family", str, sizeof str)) {
		indexFamily = converter->index_family_from_string(str);
	}
	if (table_get_string(L, 1, "day_count_fraction", str, sizeof str)) {
		dfc = converter->day_count_fraction_from_string(str);
	}
	if (table_get_string(L, 1, "curve_type", str, sizeof str)) {
		curve_type = converter->pricing_curve_type_from_string(str);
	}
	if (table_get_string(L, 1, "value_type", str, sizeof str)) {
		rate_type = converter->rate_type_from_string(str);
	}
	if (table_get_string(L, 1, "tenor", str, sizeof str) && str[0]) {
		Period p;
		if (converter->period_from_string(str, &p)) {
			tenor = converter->tenor_from_period_unit_and_len(p.units(), p.length());
		}
	}

	constexpr size_t LEN = 50;
	Date xvalues[LEN];
	double yvalues[LEN];

	auto xn = table_get_intarray(L, 1, "maturities", xvalues, LEN);
	auto yn = table_get_doublearray(L, 1, "values", yvalues, LEN);

	const char *errmsg = nullptr;
	if (xn < 4) {
		errmsg = "at least 4 values are required";
	} else if (xn != yn) {
		errmsg = "maturities and values must have same number of elements";
	}

	if (!errmsg) {
		int *x = xvalues;
		for (int i = 0; i < xn; i++) {
			if (!is_valid_date(x[i])) {
				errmsg = "invalid maturity date found";
				break;
			}
			if ((i == 0 && x[0] <= as_of_date) || (i > 0 && x[i] <= x[i - 1])) {
				errmsg = "maturities must be in ascending order";
				break;
			}
			if (rate_type == IRRateType::DISCOUNT_FACTOR) {
				if (yvalues[i] <= 0 || yvalues[i] > 1.00001) {
					errmsg = "invalid discount factor found";
				}
			} else if (rate_type == IRRateType::ZERO_RATE) {
				if (yvalues[i] < -0.02 || yvalues[i] > 0.2) {
					errmsg = "suspect zero rate found";
				}
			}
		}
	}

	if (!is_valid_date(as_of_date) || as_of_date > xvalues[0] || as_of_date < (xvalues[0] - 366)) {
		errmsg = "invalid as of date";
	}

	if (!errmsg) {
		auto curve_id = make_curve_id(curve_type, ccy, indexFamily, tenor, as_of_date, 0, qual, 0);
		std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve =
		    make_curve(&GlobalAllocator, curve_id, as_of_date, xvalues, yvalues, xn, it, rate_type, order, dfc);
		if (curve) {
			IRCurveHolder *iptr = (IRCurveHolder *)lua_newuserdata(L, sizeof(IRCurveHolder));
			new (iptr) IRCurveHolder(std::move(curve));
			luaL_getmetatable(L, Type_IRCurve);
			lua_setmetatable(L, -2);
		} else {
			errmsg = "failed to construct the curve requested";
		}
	}
	if (errmsg)
		luaL_argerror(L, 1, errmsg);
	return errmsg ? 0 : 1;
}

static int curve_tostring(lua_State *L)
{
	IRCurveHolder *iptr = check_Type_IRCurve(L, 1);
	char name[100];
	snprintf(name, sizeof name, "%s %s", Type_IRCurve, iptr->curve->name().c_str());
	lua_pushstring(L, name);
	return 1;
}

static int curve_getvalues(lua_State *L)
{
	IRCurveHolder *iptr = check_Type_IRCurve(L, 1);
	size_t n = iptr->curve->last_pillar();
	size_t len;
	ravi_create_integer_array(L, n, 0);
	lua_Integer *m = ravi_get_integer_array_rawdata(L, -1, &len);
	ravi_create_number_array(L, n, 0.0);
	lua_Number *rates = ravi_get_number_array_rawdata(L, -1, &len);
	ravi_create_number_array(L, n, 0.0);
	lua_Number *factors = ravi_get_number_array_rawdata(L, -1, &len);
	for (int i = 1; i <= n; i++) {
		m[i] = iptr->curve->maturity_date(i);
		double t = iptr->curve->maturity_time(i);
		rates[i] = iptr->curve->zero_rate(t);
		factors[i] = iptr->curve->discount(t);
	}
	return 3;
}

static int curve_sensitivities(lua_State *L)
{
	IRCurveHolder *iptr = check_Type_IRCurve(L, 1);
	lua_Number n;
	const char *errmsg = nullptr;
	if (lua_isinteger(L, 2)) {
		int x = (int)lua_tointeger(L, 2);
		if (is_valid_date(x) && x >= iptr->curve->as_of_date() && x <= iptr->curve->last_maturity()) {
			n = iptr->curve->time_from_reference((Date)x);
		} else {
			errmsg = "invalid date";
		}

	} else {
		n = luaL_checknumber(L, 2);
		if (n <= 0.0 || n > iptr->curve->last_maturity_time()) {
			errmsg = "invalid date";
		}
	}
	int rc = 0;
	if (!errmsg) {
		auto allocSet = get_threadspecific_allocators();
		auto A = allocSet->tempspace_allocator;
		FixedRegionAllocatorGuard g(A);
		auto ad = iptr->curve->get_sensitivities(n, A);
		if (ad) {
			auto result = alloc_adouble(L, ad->vars_, ad->order_, -1, 0.0);
			redukti_adouble_assign(result, ad.get());
			rc = 1;
		}
	}
	if (errmsg)
		luaL_argerror(L, 2, errmsg);
	return rc;
}

static int get_zero_rate(lua_State *L)
{
	double rate = 0;
	IRCurveHolder *iptr = check_Type_IRCurve(L, 1);
	if (lua_isinteger(L, 2)) {
		int x = (int)lua_tointeger(L, 2);
		if (is_valid_date(x) && x >= iptr->curve->as_of_date() && x <= iptr->curve->last_maturity()) {
			rate = iptr->curve->zero_rate((Date)x);
		} else {
			luaL_argerror(L, 2, "invalid date");
		}
	} else {
		lua_Number x = luaL_checknumber(L, 2);
		if (x > 0.0 && x <= iptr->curve->last_maturity_time()) {
			rate = iptr->curve->zero_rate(x);
		} else {
			luaL_argerror(L, 2, "invalid time from reference");
		}
	}
	lua_pushnumber(L, rate);
	return 1;
}

static int get_discount_factor(lua_State *L)
{
	double rate = 0;
	IRCurveHolder *iptr = check_Type_IRCurve(L, 1);
	if (lua_isinteger(L, 2)) {
		int x = (int)lua_tointeger(L, 2);
		if (is_valid_date(x) && x >= iptr->curve->as_of_date() && x <= iptr->curve->last_maturity()) {
			rate = iptr->curve->discount((Date)x);
		} else {
			luaL_argerror(L, 2, "invalid date");
		}
	} else {
		lua_Number x = luaL_checknumber(L, 2);
		if (x > 0.0 && x <= iptr->curve->last_maturity_time()) {
			rate = iptr->curve->discount(x);
		} else {
			luaL_argerror(L, 2, "invalid time from reference");
		}
	}
	lua_pushnumber(L, rate);
	return 1;
}

static int make_pricing_curve(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1,
		      "table expected"); // will do longjmp

	PricingCurveType curve_type = PricingCurveType::PRICING_CURVE_TYPE_UNSPECIFIED;
	Currency ccy = Currency::CURRENCY_UNSPECIFIED;
	IndexFamily index_family = IndexFamily::INDEX_FAMILY_UNSPECIFIED;
	Tenor tenor = Tenor::TENOR_UNSPECIFIED;

	char str[80] = {0};
	auto converter = get_default_converter();
	if (table_get_string(L, 1, "currency", str, sizeof str)) {
		ccy = converter->currency_from_string(str);
	}
	if (table_get_string(L, 1, "index_family", str, sizeof str)) {
		index_family = converter->index_family_from_string(str);
	}
	if (table_get_string(L, 1, "curve_type", str, sizeof str)) {
		curve_type = converter->pricing_curve_type_from_string(str);
	}
	if (table_get_string(L, 1, "tenor", str, sizeof str)) {
		Period p;
		if (converter->period_from_string(str, &p)) {
			tenor = converter->tenor_from_period_unit_and_len(p.units(), p.length());
		}
	}

	luaL_argcheck(L, curve_type != PricingCurveType::PRICING_CURVE_TYPE_UNSPECIFIED, 1, "invalid curve_type");
	luaL_argcheck(L, ccy != Currency::CURRENCY_UNSPECIFIED, 1, "invalid currency");
	luaL_argcheck(L, tenor <= Tenor::TENOR_12M, 1, "invalid tenor, mut be <= 12M");

	lua_pushinteger(L, PricingCurve(curve_type, ccy, index_family, tenor).id());
	return 1;
}

/* __gc for Type_CurveMapper */
static int collect_Type_PricingCashflows(lua_State *L)
{
	PricingCashflowsHolder *holder = check_Type_PricingCashflows(L, 1);
	holder->~PricingCashflowsHolder();
	return 0;
}

static int make_pricing_cashflows(lua_State *L)
{
	ValuationContextHolder *vc = check_Type_ValuationContext(L, 1);
	CurveMapperHolder *cm = check_Type_CurveMapper(L, 2);
	CFCollectionHolder *cf = check_Type_CFCollection(L, 3);

	PricingCashflowsHolder *iptr = (PricingCashflowsHolder *)lua_newuserdata(L, sizeof(PricingCashflowsHolder));
	new (iptr) PricingCashflowsHolder();
	luaL_getmetatable(L, Type_PricingCashflows);
	lua_setmetatable(L, -2);

	// Following does not retain references to the mapper or valuation context so we don't have to
	// worry about garbage collection of those.
	iptr->cashflows = construct_cashflows(&iptr->stackregion_allocator, cf->cashflows, vc->context, &cm->mapper);
	return 1;
}

static int pricing_cashflows_isvalid(lua_State *L)
{
	PricingCashflowsHolder *iptr = check_Type_PricingCashflows(L, 1);
	lua_pushboolean(L, iptr->cashflows != nullptr);
	return 1;
}

/* __gc for Type_PricingResult */
static int collect_Type_PricingResult(lua_State *L)
{
	PricingResult *holder = check_Type_PricingResult(L, 1);
	holder->~PricingResult();
	return 0;
}

// Calculate PV and if ValuationContext.derivative_order > 0
// then also delta and gamma
// extern double compute_present_value(FixedRegionAllocator *A,
//	const ValuationContext &ctx,
//	const Cashflows *flows,
//	const CurveProvider *mapping_provider,
//	Sensitivities &sensitivities,
//	StatusCode &status);
static int calculate_present_value(lua_State *L)
{
	ValuationContextHolder *context = check_Type_ValuationContext(L, 1);
	PricingCashflowsHolder *cashflows = check_Type_PricingCashflows(L, 2);
	CurveProviderHolder *provider = check_Type_CurveProvider(L, 3);
	PricingResult *result = (PricingResult *)lua_newuserdata(L, sizeof(PricingResult));
	new (result) PricingResult();
	luaL_getmetatable(L, Type_PricingResult);
	lua_setmetatable(L, -2);
	compute_present_value(&result->pricing_allocator, context->context, cashflows->cashflows, &provider->provider,
			      result->sensitivities, result->status);
	if (result->status != StatusCode::kOk) {
		error("Pricing failed because %s\n", error_message(result->status));
	} else {
		if (is_debug_enabled()) {
			result->sensitivities.dump(stdout);
		}
	}
	int top = lua_gettop(L);
	// Since the senitivities retain references to curves we need to ensure that
	// the curve provider is not garbage collected before the pricing result
	lua_newtable(L);
	lua_pushvalue(L, 3);		       // curve provider
	lua_setfield(L, -2, "curve_provider"); // store in the table
	lua_setuservalue(L, -2);	       // store table as uservalue  in pricing result
	assert(lua_gettop(L) == top);
	return 1;
}

// Returns the curves for which there are sensitivities
static int pricing_result_curve_ids(lua_State *L)
{
	PricingResult *result = check_Type_PricingResult(L, 1);
	luaL_argcheck(L, result->status == StatusCode::kOk, 1, "invalid state");

	std::vector<CurveId> curves;
	result->sensitivities.get_curve_ids(curves);

	lua_newtable(L);
	for (int i = 0; i < curves.size(); i++) {
		auto name = curve_id_to_string(curves[i]);
		lua_pushstring(L, name.c_str());
		lua_seti(L, -2, curves[i]);
	}
	return 1;
}

static int pricing_result_get_1d_sensitivities(lua_State *L)
{
	PricingResult *result = check_Type_PricingResult(L, 1);
	luaL_argcheck(L, result->status == StatusCode::kOk, 1, "invalid state");
	CurveId id = luaL_checkinteger(L, 2);

	auto sens = result->sensitivities.find_first_order_sensitivities(id);
	if (sens == nullptr || sens->count() == 0)
		return 0;

	ravi_create_number_array(L, sens->count() - 1, 0.0);
	size_t len;
	double *ptr = ravi_get_number_array_rawdata(L, -1, &len);
	for (int i = 1; i < sens->count(); i++) {
		ptr[i] = sens->at(i);
	}
	return 1;
}

static int pricing_result_isvalid(lua_State *L)
{
	PricingResult *result = check_Type_PricingResult(L, 1);
	lua_pushboolean(L, result->status == StatusCode::kOk);
	return 1;
}

static const char *build_curve_definitions(lua_State *L, int idx, BootstrapCurvesRequest *request, Arena *arena)
{
	char str[256] = {0};
	int temp = 0;
	const char *errmsg = nullptr;
	auto converter = get_default_converter();
	for (int i = 0; !errmsg && i < lua_rawlen(L, idx); i++) {
		lua_rawgeti(L, idx, i + 1);
		if (lua_istable(L, -1)) {
			int id = 0;
			InterpolatorType it = InterpolatorType::LINEAR; // default
			Currency ccy = Currency::CURRENCY_UNSPECIFIED;
			IndexFamily indexFamily = IndexFamily::INDEX_FAMILY_UNSPECIFIED;
			Tenor tenor = Tenor::TENOR_UNSPECIFIED;
			MarketDataQualifier qual = MarketDataQualifier::MDQ_NORMAL;
			IRRateType interpolated_on = IRRateType::ZERO_RATE;
			MaturityGenerationRule rule =
			    MaturityGenerationRule::MATURITY_GENERATION_RULE_DERIVE_FROM_INSTRUMENTS;
			CurveGroup group = CurveGroup::CURVE_GROUP_UNSPECIFIED;
			buffer_ptr<char> buf(1024);
			std::vector<const char *> tenors;
			if (table_get_string(L, -1, "group", str, sizeof str)) {
				group = converter->curve_group_from_string(str);
			}
			if (table_get_integer(L, -1, "curve_id", &temp)) {
				id = temp;
			}
			if (table_get_string(L, -1, "interpolator", str, sizeof str)) {
				it = converter->interpolator_type_from_string(str);
			}
			if (table_get_string(L, -1, "currency", str, sizeof str)) {
				ccy = converter->currency_from_string(str);
			}
			if (table_get_string(L, -1, "index_family", str, sizeof str)) {
				indexFamily = converter->index_family_from_string(str);
			}
			if (table_get_string(L, -1, "interpolated_on", str, sizeof str)) {
				interpolated_on = converter->rate_type_from_string(str);
			}
			if (table_get_string(L, -1, "curve_tenor", str, sizeof str) && str[0]) {
				Period p;
				if (converter->period_from_string(str, &p)) {
					tenor = converter->tenor_from_period_unit_and_len(p.units(), p.length());
				}
			}
			if (table_get_string(L, -1, "maturity_generation_rule", str, sizeof str)) {
				rule = converter->maturity_generation_rule_from_string(str);
			}
			if (rule == MaturityGenerationRule::MATURITY_GENERATION_RULE_FIXED_TENORS) {
				if (table_get_string(L, -1, "tenors", str, sizeof str)) {
					parse_delimited(str, str + strlen(str), tenors, buf, ":");
				}
			}
			if (ccy == Currency::CURRENCY_UNSPECIFIED ||
			    indexFamily == IndexFamily::INDEX_FAMILY_UNSPECIFIED ||
			    group == CurveGroup::CURVE_GROUP_UNSPECIFIED || id == 0) {
				errmsg = "validation of curve definition failed";
			} else {
				auto defn = request->add_curve_definitions();
				defn->set_currency(ccy);
				defn->set_curve_group(group);
				defn->set_id(id);
				defn->set_index_family(indexFamily);
				defn->set_tenor(tenor);
				defn->set_interpolated_on(interpolated_on);
				defn->set_interpolator_type(it);
				defn->set_maturity_generation_rule(rule);
				for (int i = 0; i < tenors.size(); i++) {
					Period p;
					if (converter->period_from_string(tenors[i], &p)) {
						auto t =
						    converter->tenor_from_period_unit_and_len(p.units(), p.length());
						if (t == TENOR_UNSPECIFIED) {
							errmsg = "failed to parse "
								 "tenor";
							break;
						} else {
							defn->add_tenors(t);
						}
					} else {
						errmsg = "failed to parse tenor";
						break;
					}
				}
			}
		} else {
			errmsg = "unexpected value, table expected";
		}
		lua_pop(L, 1);
	}
	return errmsg;
}

static const char *build_parcurve_set(lua_State *L, int idx, ParCurveSet *parset, Arena *arena)
{
	char str[256] = {0};
	int temp = 0;
	const char *errmsg = nullptr;
	auto converter = get_default_converter();
	ParCurve *curve = nullptr;
	int current_curve_id = -1;
	for (int i = 0; !errmsg && i < lua_rawlen(L, idx); i++) {
		lua_rawgeti(L, idx, i + 1);
		if (lua_istable(L, -1)) {
			int id = 0;
			int forward_curve_id = 0;
			int discount_curve_id = 0;
			double par_rate = 0;
			char instrument_type[80] = {0};
			char instrument_id[80] = {0};
			Tenor tenor = TENOR_UNSPECIFIED;
			if (table_get_integer(L, -1, "curve_id", &temp)) {
				id = temp;
			}
			if (table_get_string(L, -1, "instrument_type", str, sizeof str)) {
				string_copy(instrument_type, str, sizeof instrument_type);
			}
			if (table_get_string(L, -1, "instrument_id", str, sizeof str)) {
				string_copy(instrument_id, str, sizeof instrument_id);
			}
			if (!table_get_number(L, -1, "par_rate", &par_rate, -1.0)) {
				errmsg = "failed to get par rate";
			}
			if (table_get_integer(L, -1, "forward_curve_id", &temp)) {
				forward_curve_id = temp;
			}
			if (table_get_integer(L, -1, "discount_curve_id", &temp)) {
				discount_curve_id = temp;
			}
			if (table_get_string(L, -1, "floating_tenor", str, sizeof str) && str[0]) {
				Period p;
				if (converter->period_from_string(str, &p)) {
					tenor = converter->tenor_from_period_unit_and_len(p.units(), p.length());
				}
			}
			if (id != current_curve_id) {
				curve = parset->add_par_curves();
				curve->set_curve_definition_id(id);
				auto par_rates = Arena::CreateMessage<ParRates>(arena);
				curve->set_allocated_par_rates(par_rates);
				current_curve_id = id;
			}
			assert(curve);
			auto inst = curve->add_instruments();
			inst->set_instrument_key(instrument_id);
			inst->set_instrument_type(instrument_type);
			inst->set_floating_tenor(tenor);
			inst->set_forward_curve_definition_id(forward_curve_id);
			inst->set_discount_curve_definition_id(discount_curve_id);
			curve->mutable_par_rates()->add_values(par_rate);
		} else {
			errmsg = "unexpected value, table expected";
		}
		lua_pop(L, 1);
	}
	return errmsg;
}

static int build_curves(lua_State *L)
{
	int as_of_date = luaL_checkinteger(L, 1);
	luaL_argcheck(L, lua_istable(L, 2), 2, "table expected");
	luaL_argcheck(L, lua_istable(L, 3), 3, "table expected");
	luaL_argcheck(L, is_valid_date(as_of_date), 1, "invalid date");
	const char *errmsg = nullptr;
	char temp[256];
	char pricing_script[256];
	if (lua_gettop(L) == 4 && lua_isstring(L, 4)) {
		snprintf(pricing_script, sizeof pricing_script, "%s", lua_tostring(L, 4));
	} else {
		string_copy(pricing_script, PRICING_SCRIPT, sizeof pricing_script);
	}

	ZeroCurveSetHolder *ptr = (ZeroCurveSetHolder *)lua_newuserdata(L, sizeof(ZeroCurveSetHolder));
	new (ptr) ZeroCurveSetHolder();
	luaL_getmetatable(L, Type_ZeroCurveSet);
	lua_setmetatable(L, -2);

	auto bootstrap_request = Arena::CreateMessage<BootstrapCurvesRequest>(&ptr->arena);
	ptr->request = bootstrap_request; // Save so that we cab=n later on access the definitions etc

	errmsg = build_curve_definitions(L, 2, bootstrap_request, &ptr->arena);
	if (!errmsg) {
		auto parset = Arena::CreateMessage<ParCurveSet>(&ptr->arena);
		bootstrap_request->set_allocated_par_curve_set(parset);
		parset->set_as_of_date(as_of_date);
		errmsg = build_parcurve_set(L, 3, parset, &ptr->arena);
	}
	if (!errmsg) {
		bootstrap_request->set_business_date(as_of_date);
		// FIXME these should be obtained from parameters (options?)
		bootstrap_request->set_generate_par_sensitivities(false);
		bootstrap_request->set_max_solver_iterations(30);
		bootstrap_request->set_solver_type(SolverType::SOLVER_TYPE_LEVENBERG_MARQUARDT);

		std::unique_ptr<CurveBuilderService> bootstrapper =
		    get_curve_builder_service(std::string(pricing_script));
		auto reply = Arena::CreateMessage<BootstrapCurvesReply>(&ptr->arena);
		bootstrapper->handle_bootstrap_request(bootstrap_request, reply);
		assert(reply->has_header());
		if (reply->header().response_code() != StandardResponseCode::SRC_OK) {
			string_copy(temp, reply->header().response_message().c_str(), sizeof temp);
			errmsg = temp;
		} else {
			ptr->zero_curve_set = Arena::CreateMessage<ZeroCurveSet>(&ptr->arena);
			// FIXME we should call one of the arena swap / unsafe methods here
			ptr->zero_curve_set->mutable_zero_curves()->CopyFrom(reply->curves());
			ptr->zero_curve_set->mutable_zero_par_sensitivities()->CopyFrom(reply->par_sensitivities());
		}
	}

	int rc = 1;
	if (errmsg) {
		lua_pushstring(L, errmsg);
		rc++;
	}

	return rc;
}

/* __gc for Type_ZeroCurveSet */
static int collect_Type_ZeroCurveSet(lua_State *L)
{
	ZeroCurveSetHolder *holder = check_Type_ZeroCurveSet(L, 1);
	holder->~ZeroCurveSetHolder();
	return 0;
}

static int zerocurveset_getcurves(lua_State *L)
{
	ZeroCurveSetHolder *holder = check_Type_ZeroCurveSet(L, 1);
	luaL_argcheck(L, holder->zero_curve_set != nullptr, 1, "ZeroCurveSet is not valid");
	lua_newtable(L);
	for (int i = 0; i < holder->zero_curve_set->zero_curves_size(); i++) {
		auto id = holder->zero_curve_set->zero_curves(i).curve_definition_id();
		const IRCurveDefinition *defn = nullptr;
		for (int j = 0; j < holder->request->curve_definitions_size(); j++) {
			if (holder->request->curve_definitions(j).id() == id) {
				defn = &holder->request->curve_definitions(j);
				break;
			}
		}
		if (defn == nullptr)
			// FIXME should be an error
			continue;
		std::unique_ptr<YieldCurve, Deleter<YieldCurve>> curve =
		    make_curve(holder->request->business_date(), defn, holder->zero_curve_set->zero_curves(i), 0);
		if (curve) {
			IRCurveHolder *iptr = (IRCurveHolder *)lua_newuserdata(L, sizeof(IRCurveHolder));
			new (iptr) IRCurveHolder(std::move(curve));
			luaL_getmetatable(L, Type_IRCurve);
			lua_setmetatable(L, -2);
			lua_seti(L, -2, id);
		} else {
			// FIXME error handling
		}
	}
	return 1;
}

static int zerocurveset_isvalid(lua_State *L)
{
	ZeroCurveSetHolder *holder = check_Type_ZeroCurveSet(L, 1);
	lua_pushboolean(L, holder->zero_curve_set != nullptr);
	return 1;
}

static int make_fixingservice(lua_State *L)
{
	luaL_argcheck(L, lua_istable(L, 1), 1, "table expected");

	FixingServiceHolder *iptr = (FixingServiceHolder *)lua_newuserdata(L, sizeof(FixingServiceHolder));
	new (iptr) FixingServiceHolder();
	luaL_getmetatable(L, Type_FixingService);
	lua_setmetatable(L, -2);

	int pos = lua_gettop(L);

	lua_pushnil(L);
	while (lua_next(L, -3)) {
		// stack now contains: -1 => value; -2 => key; -3 => table
		// key must be index
		// value must be a table
		auto index = test_Type_Index(L, -2);
		if (index != nullptr && lua_istable(L, -1)) {
			std::vector<Value> values;
			// no loop through the values
			lua_pushnil(L);
			while (lua_next(L, -2)) {
				// stack now should have
				// -1 => fixing value
				// -2 => fixing key (date)
				if (lua_isinteger(L, -2) && lua_isnumber(L, -1)) {
					int date = (int)lua_tointeger(L, -2);
					lua_Number value = lua_tonumber(L, -1);
					if (is_valid_date(date) && value >= -0.01 && value <= 0.2) {
						values.push_back(Value(date, value));
					}
				}
				// pop the value
				lua_pop(L, 1);
			}
			std::unique_ptr<TimeSeries> ts = std::unique_ptr<TimeSeries>(new TimeSeries());
			ts->add(values.size(), values.data());
			iptr->fixing_service.set_fixings(index->index->id(), std::move(ts));
		}
		// pop the value
		lua_pop(L, 1);
	}

	assert(lua_gettop(L) == pos);
	// We have fixing service at -1
	// table parameter at -2
	// Put a reference to the table in fixing service
	assert(lua_istable(L, -2));
	lua_pushvalue(L, -2); // now table is at -1, fixing service at -2
	lua_setuservalue(L, -2);
	assert(lua_gettop(L) == pos);

	return 1;
}

static int get_fixing(lua_State *L)
{
	FixingServiceHolder *holder = check_Type_FixingService(L, 1);
	IndexHolder *iptr = check_Type_Index(L, 2);
	int date = (int)luaL_checkinteger(L, 3);

	auto ts = holder->fixing_service.get_fixings(iptr->index->id());

	if (ts) {
		double value;
		if (ts->find(date, value)) {
			lua_pushnumber(L, value);
		} else {
			lua_pushnil(L);
		}
	} else {
		lua_pushnil(L);
	}
	return 1;
}

static int collect_Type_FixingService(lua_State *L)
{
	FixingServiceHolder *holder = check_Type_FixingService(L, 1);
	holder->~FixingServiceHolder();
	return 0;
}

#define Kint(v, s)                                                                                                     \
	lua_pushinteger(L, v);                                                                                         \
	lua_setfield(L, -2, s)

#define Kstring(v, s)                                                                                                  \
	lua_pushstring(L, s);                                                                                          \
	lua_setfield(L, -2, s)

// redukti table must be at the top of the stack
static void setup_constants(lua_State *L)
{
	int top = lua_gettop(L);
	lua_newtable(L);

	Kint(PRICING_CURVE_TYPE_FORWARD, "FORWARD_CURVE");
	Kint(PRICING_CURVE_TYPE_DISCOUNT, "DISCOUNT_CURVE");

	Kint(CURRENCY_UNSPECIFIED, "CURRENCY_UNSPECIFIED");
	Kint(USD, "USD");
	Kint(GBP, "GBP");
	Kint(EUR, "EUR");
	Kint(JPY, "JPY");
	Kint(AUD, "AUD");
	Kint(CAD, "CAD");
	Kint(CHF, "CHF");
	Kint(CZK, "CZK");
	Kint(DKK, "DKK");
	Kint(HKD, "HKD");
	Kint(NOK, "NOK");
	Kint(NZD, "NZD");
	Kint(SEK, "SEK");
	Kint(SGD, "SGD");
	Kint(MXN, "MXN");
	Kint(ZAR, "ZAR");
	Kint(HUF, "HUF");
	Kint(PLN, "PLN");

	Kint(INDEX_FAMILY_UNSPECIFIED, "INDEX_FAMILY_UNSPECIFIED");
	Kint(LIBOR, "LIBOR");
	Kint(FEDFUND, "FEDFUND");
	Kint(EURIBOR, "EURIBOR");
	Kint(SONIA, "SONIA");
	Kint(AONIA, "AONIA");
	Kint(BBSW, "BBSW");
	Kint(CDOR, "CDOR");
	Kint(CORRA, "CORRA");
	Kint(TOIS, "TOIS");
	Kint(PRIBOR, "PRIBOR");
	Kint(CIBOR, "CIBOR");
	Kint(HIBOR, "HIBOR");
	Kint(BUBOR, "BUBOR");
	Kint(TONA, "TONA");
	Kint(NIBOR, "NIBOR");
	Kint(BKBM, "BKBM");
	Kint(WIBOR, "WIBOR");
	Kint(STIBOR, "STIBOR");
	Kint(SOR, "SOR");
	Kint(JIBAR, "JIBAR");
	Kint(EONIA, "EONIA");
	Kint(BBR, "BBR");
	Kint(BA, "BA");
	Kint(USDMXNBASIS, "USDMXNBASIS");
	Kint(TIIE, "TIIE");

	Kint(TENOR_UNSPECIFIED, "TENOR_UNSPECIFIED");
	Kint(TENOR_1D, "TENOR_1D");
	Kint(TENOR_2D, "TENOR_2D");
	Kint(TENOR_3D, "TENOR_3D");
	Kint(TENOR_4D, "TENOR_4D");
	Kint(TENOR_5D, "TENOR_5D");
	Kint(TENOR_6D, "TENOR_6D");
	Kint(TENOR_1W, "TENOR_1W");
	Kint(TENOR_2W, "TENOR_2W");
	Kint(TENOR_3W, "TENOR_3W");
	Kint(TENOR_1M, "TENOR_1M");
	Kint(TENOR_2M, "TENOR_2M");
	Kint(TENOR_3M, "TENOR_3M");
	Kint(TENOR_4M, "TENOR_4M");
	Kint(TENOR_5M, "TENOR_5M");
	Kint(TENOR_6M, "TENOR_6M");
	Kint(TENOR_7M, "TENOR_7M");
	Kint(TENOR_8M, "TENOR_8M");
	Kint(TENOR_9M, "TENOR_9M");
	Kint(TENOR_10M, "TENOR_10M");
	Kint(TENOR_11M, "TENOR_11M");
	Kint(TENOR_12M, "TENOR_12M");

	lua_setfield(L, -2, "Constants");
}

static const luaL_Reg reduktilib[] = {{"date", make_date},
				      {"date_parts", date_parts},
				      {"addperiod", add_period},
				      {"schedule", build_schedule},
				      {"cashflows", build_cfcollection},
				      {"index", get_index},
				      {"dayfraction", get_day_count_fraction},
				      {"calendar", get_calendar},
				      {"future_imm_dates", future_imm_dates},
				      {"fra_tenors", fra_tenors},
				      {"loadcsv", read_csv},
				      {"loadmatrix", read_matrix},
				      {"stringreplace", string_replace},
				      {"mean", stats_mean},
				      {"stddev", stats_stddev},
				      {"variance", stats_variance},
				      {"lagged_difference", stats_lagged_difference},
				      {"volatility", stats_volatility},
				      {"adouble1", create_adouble_order1},
				      {"adouble2", create_adouble_order2},
				      {"sqrt", adouble_sqrt},
				      {"pow", adouble_power},
				      {"log", adouble_log},
				      {"exp", adouble_exp},
				      {"abs", adouble_abs},
				      {"cos", adouble_cos},
				      {"sin", adouble_sin},
				      {"tan", adouble_tan},
				      {"interpolator", create_interpolator},
				      {"curve", create_ircurve},
				      {"pricing_curve_id", make_pricing_curve},
				      {"curve_mapper", make_curve_mapper},
				      {"valuation_context", make_valuation_context},
				      {"prepare_cashflows_for_pricing", make_pricing_cashflows},
				      {"curve_provider", make_curve_provider},
				      {"present_value", calculate_present_value},
				      {"build_curves", build_curves},
				      {"fixing_service", make_fixingservice},
				      {NULL, NULL}};

static const luaL_Reg index_methods[] = {{"maturity", index_maturity},
					 {"date_components", index_date_components},
					 {"adjust_date", index_adjust_date},
					 {NULL, NULL}};

static const luaL_Reg dayfraction_methods[] = {{"fraction", day_fraction_evaluate}, {NULL, NULL}};

static const luaL_Reg calendar_methods[] = {{"advance", calendar_advance}, {NULL, NULL}};

static const luaL_Reg interpolator_methods[] = {{"interpolate", interpolator_get}, {NULL, NULL}};

static const luaL_Reg pricing_cashflows_methods[] = {{"ok", pricing_cashflows_isvalid}, {NULL, NULL}};

static const luaL_Reg pricing_result_methods[] = {{"ok", pricing_result_isvalid},
						  {"curve_ids", pricing_result_curve_ids},
						  {"delta", pricing_result_get_1d_sensitivities},
						  {NULL, NULL}};

static const luaL_Reg fixing_service_methods[] = {{"fixing", get_fixing}, {NULL, NULL}};

static const luaL_Reg zerocurveset_methods[] = {
    {"ok", zerocurveset_isvalid}, {"curves", zerocurveset_getcurves}, {NULL, NULL}};

static const luaL_Reg curve_mapper_methods[] = {
    {"add_mapping", add_curve_mapping}, {"get_mapping", get_curve_mapping}, {NULL, NULL}};

static const luaL_Reg curve_provider_methods[] = {
    {"add_mapping", add_curve_provider_mapping}, {"curves", curve_provider_get_curves}, {NULL, NULL}};

static const luaL_Reg ircurve_methods[] = {{"zero_rate", get_zero_rate},
					   {"discount_factor", get_discount_factor},
					   {"sensitivities", curve_sensitivities},
					   {"values", curve_getvalues},
					   {NULL, NULL}};

static const luaL_Reg adouble_methods[] = {{"gradient", adouble_get_gradient},
					   {"hessian", adouble_get_hessian},
					   {"abs", adouble_abs},
					   {"pow", adouble_power},
					   {"sqrt", adouble_sqrt},
					   {"exp", adouble_exp},
					   {"log", adouble_log},
					   {"cos", adouble_cos},
					   {"sin", adouble_sin},
					   {"tan", adouble_tan},
					   {NULL, NULL}};

LUAMOD_API int raviopen_redukti(lua_State *L)
{
	luaL_newmetatable(L, Type_CFCollection);
	lua_pushcfunction(L, collect_Type_CFCollection);
	lua_setfield(L, -2, "__gc");
	lua_pushcfunction(L, cfcollection_tostring);
	lua_setfield(L, -2, "__tostring");
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_ValuationContext);
	lua_pushcfunction(L, collect_Type_ValuationContext);
	lua_setfield(L, -2, "__gc");
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_FixingService);
	lua_pushcfunction(L, collect_Type_FixingService);
	lua_setfield(L, -2, "__gc");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, fixing_service_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_Interpolator);
	lua_pushcfunction(L, collect_Type_Interpolator);
	lua_setfield(L, -2, "__gc");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, interpolator_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_PricingCashflows);
	lua_pushcfunction(L, collect_Type_PricingCashflows);
	lua_setfield(L, -2, "__gc");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, pricing_cashflows_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_CurveMapper);
	lua_pushcfunction(L, collect_Type_CurveMapper);
	lua_setfield(L, -2, "__gc");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, curve_mapper_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_CurveProvider);
	lua_pushcfunction(L, collect_Type_CurveProvider);
	lua_setfield(L, -2, "__gc");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, curve_provider_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_PricingResult);
	lua_pushcfunction(L, collect_Type_PricingResult);
	lua_setfield(L, -2, "__gc");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, pricing_result_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_IRCurve);
	lua_pushcfunction(L, collect_Type_IRCurve);
	lua_setfield(L, -2, "__gc");
	lua_pushcfunction(L, curve_tostring);
	lua_setfield(L, -2, "__tostring");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, ircurve_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_ZeroCurveSet);
	lua_pushcfunction(L, collect_Type_ZeroCurveSet);
	lua_setfield(L, -2, "__gc");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, zerocurveset_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_Index);
	lua_pushcfunction(L, index_tostring);
	lua_setfield(L, -2, "__tostring");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, index_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_DayFraction);
	lua_pushcfunction(L, day_fraction_tostring);
	lua_setfield(L, -2, "__tostring");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, dayfraction_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_Calendar);
	lua_pushcfunction(L, calendar_tostring);
	lua_setfield(L, -2, "__tostring");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, calendar_methods, 0);
	lua_pop(L, 1);

	luaL_newmetatable(L, Type_ADouble);
	lua_pushcfunction(L, adouble_tostring);
	lua_setfield(L, -2, "__tostring");
	lua_pushcfunction(L, adouble_add);
	lua_setfield(L, -2, "__add");
	lua_pushcfunction(L, adouble_sub);
	lua_setfield(L, -2, "__sub");
	lua_pushcfunction(L, adouble_mul);
	lua_setfield(L, -2, "__mul");
	lua_pushcfunction(L, adouble_div);
	lua_setfield(L, -2, "__div");
	lua_pushcfunction(L, adouble_unm);
	lua_setfield(L, -2, "__unm");
	lua_pushcfunction(L, adouble_eq);
	lua_setfield(L, -2, "__eq");
	lua_pushcfunction(L, adouble_lt);
	lua_setfield(L, -2, "__lt");
	lua_pushcfunction(L, adouble_le);
	lua_setfield(L, -2, "__le");
	lua_pushcfunction(L, adouble_get);
	lua_setfield(L, -2, "__call");
	lua_pushvalue(L, -1);		/* push metatable */
	lua_setfield(L, -2, "__index"); /* metatable.__index = metatable */
	luaL_setfuncs(L, adouble_methods, 0);
	lua_pop(L, 1);

	luaL_newlib(L, reduktilib);
	setup_constants(L);
	return 1;
}
}

namespace redukti
{

std::unique_ptr<CFCollection> raviapi_get_CFCollection(lua_State *L, int idx)
{
	CFCollectionHolder *holder = test_Type_CFCollection(L, idx);
	CFCollection *ptr = nullptr;
	if (holder != nullptr) {
		ptr = holder->cashflows;
		holder->cashflows = nullptr;
		return std::unique_ptr<CFCollection>(ptr);
	}
	return std::unique_ptr<CFCollection>();
}

} // namespace redukti
