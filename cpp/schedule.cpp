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

#include <schedule.h>
#include <schedule.pb.h>

#include <calendars.h>
#include <date.h>
#include <logger.h>

#include <algorithm>
#include <stdexcept>

namespace redukti
{

Date adjust_date(Date d, RollConvention rc, int roll_date) noexcept
{
	auto ymd = date_components(d);
	switch (rc) {
	case EOM:
		return end_of_month(ymd.y, ymd.m);
	case IMM:
		// IMM - Third Wednesday of the month
		return make_date(nth_weekday(3, Wednesday, ymd.m, ymd.y));
	case IMMAUD:
		// IMMAUD - 1 Sydney day before 2nd Friday
		// (http://www.sfe.com.au/content/sfe/trading/con_specs.pdf)
		return make_date(nth_weekday(2, Friday, ymd.m, ymd.y)) - 1;
	case IMMCAD:
		// IMMCAD - 2 days before 3rd Wednesday
		return make_date(nth_weekday(3, Wednesday, ymd.m, ymd.y)) - 2;
	case IMMNZD:
		// IMMNZD - 1st Wednesday following 9th of the month
		// (http://www.sfe.com.au/content/sfe/trading/con_specs.pdf)
		return next_weekday(make_date(10, ymd.m, ymd.y), Wednesday);
	default:
		if (roll_date == 31)
			return end_of_month(ymd.y, ymd.m);
		int temp = last_day_of_month(ymd.y, ymd.m);
		if (roll_date > temp)
			roll_date = temp;
		if (roll_date == ymd.d)
			return d;
		else
			return make_date(roll_date, ymd.m, ymd.y);
	}
}

class CalculationPeriod
{
	private:
	Date unadjustedStart_;
	Date unadjustedEnd_;
	int flags_;
	Date adjustedStart_;
	Date adjustedEnd_;
	Date paymentDate_;
	Date adjustedPayment_;

	public:
	CalculationPeriod(Date start, Date end, bool stub)
	    : unadjustedStart_(start), unadjustedEnd_(end), flags_(0), adjustedStart_(start), adjustedEnd_(end),
	      paymentDate_(end), adjustedPayment_(end)
	{
		if (stub)
			flags_ |= 1;
	}
	Date unadjusted_start() const noexcept { return unadjustedStart_; }
	Date unadjusted_end() const noexcept { return unadjustedEnd_; }
	bool is_stub() const noexcept { return is_set(1); }
	bool is_payment() const noexcept { return is_set(2); }
	CalculationPeriod &payment(Date d) noexcept
	{
		paymentDate_ = d;
		return set(2);
	}
	Date payment_date() const noexcept { return paymentDate_; }
	CalculationPeriod &adjusted_start(Date d) noexcept
	{
		adjustedStart_ = d;
		return set(4);
	}
	Date adjusted_start() const noexcept { return adjustedStart_; }
	CalculationPeriod &adjusted_end(Date d) noexcept
	{
		adjustedEnd_ = d;
		return set(8);
	}
	Date adjusted_end() const noexcept { return adjustedEnd_; }
	CalculationPeriod &adjusted_payment(Date d) noexcept
	{
		adjustedPayment_ = d;
		return set(16);
	}
	Date adjusted_payment() const noexcept { return adjustedPayment_; }

	private:
	CalculationPeriod &set(int value) noexcept
	{
		flags_ |= value;
		return *this;
	}
	bool is_set(int value) const noexcept { return (flags_ & value) != 0; }
};

class ScheduleHelper
{
	public:
	StatusCode build_zero_coupon(ScheduleParameters &params, std::vector<CalculationPeriod> &periods) noexcept
	{
		periods.push_back(CalculationPeriod(params.effective_date(), params.termination_date(), false));
		return StatusCode::kOk;
	}

	/**
	 * A zero coupon schedule has one payment period.
	 */
	bool is_zero_coupon(ScheduleParameters &params) noexcept
	{
		if (params.payment_frequency() == TENOR_1T) {
			return true;
		}
		if (params.payment_frequency() != TENOR_UNSPECIFIED) {
			Period p = Period::tenor_to_period(params.payment_frequency());
			Date end = add(params.effective_date(), p);
			if (end >= params.termination_date()) {
				params.set_payment_frequency(TENOR_1T);
				return true;
			}
		}
		return false;
	}

	bool has_front_stub(ScheduleParameters &params) noexcept
	{
		if (is_valid_date(params.first_regular_period_start_date())) {
			return true;
		}
		if (params.stub_location() == SHORT_FRONT_STUB || params.stub_location() == LONG_FRONT_STUB) {
			return true;
		}
		return false;
	}

	bool has_back_stub(ScheduleParameters &params) noexcept
	{
		if (is_valid_date(params.last_regular_period_end_date())) {
			return true;
		}
		if (params.stub_location() == SHORT_BACK_STUB || params.stub_location() == LONG_BACK_STUB) {
			return true;
		}
		return false;
	}

	bool has_no_stubs(ScheduleParameters &params) noexcept
	{
		if (!is_valid_date(params.first_regular_period_start_date()) &&
		    !is_valid_date(params.last_regular_period_end_date())) {
			return true;
		}
		return false;
	}

	StatusCode build(ScheduleParameters &params, std::vector<CalculationPeriod> &periods) noexcept
	{
		periods.clear();
		auto code = validate(params);
		if (code != StatusCode::kOk) {
			return code;
		}
		StatusCode result = StatusCode::kOk;
		if (is_zero_coupon(params)) {
			result = build_zero_coupon(params, periods);
		} else if (has_front_stub(params)) {
			result = generate_backward(params, periods);
		} else if (has_back_stub(params)) {
			result = generate_forward(params, periods);
		} else {
			result = generate_backward(params, periods);
		}
		// Fix up stub period - i.e. make long stub if requested
		return result;
	}

	/**
	 * Generate schedule from front to back - useful when
	 * front stub is present
	 */
	StatusCode generate_forward(ScheduleParameters &params, std::vector<CalculationPeriod> &periods) noexcept
	{
		Date regStart = params.effective_date();
		if (is_valid_date(params.first_regular_period_start_date()) &&
		    params.first_regular_period_start_date() > regStart) {
			regStart = params.first_regular_period_start_date();
		}
		if (regStart > params.effective_date()) {
			// initial stub
			periods.push_back(CalculationPeriod(params.effective_date(), regStart, true));
		}
		if (params.roll_convention() != RollConvention::ROLL_CONVENTION_NONE &&
		    params.roll_convention() != RollConvention::ROLL_CONVENTION_DEFAULT) {
			if (regStart != adjust_date(regStart, params.roll_convention(), 0)) {
				return StatusCode::kSCH_FirstRegularStartDateInconsistent;
			}
		}
		int roll_date = date_components(regStart).d;
		Date regEnd = params.termination_date();
		if (is_valid_date(params.last_regular_period_end_date()) &&
		    params.last_regular_period_end_date() < regEnd) {
			regEnd = params.last_regular_period_end_date();
		}
		while (regStart < regEnd) {
			Date periodEnd = add(regStart, Period::tenor_to_period(params.calculation_frequency()));
			periodEnd = adjust_date(periodEnd, params.roll_convention(), roll_date);
			if (periodEnd > regEnd) {
				break;
			}
			periods.push_back(CalculationPeriod(regStart, periodEnd, false));
			regStart = periodEnd;
		}
		if (regStart < params.termination_date()) {
			periods.push_back(CalculationPeriod(regStart, params.termination_date(), true));
		}
		return StatusCode::kOk;
	}

	/**
	 * Generate schedule from back to front - useful when
	 * back stub is present
	 */
	StatusCode generate_backward(ScheduleParameters &params, std::vector<CalculationPeriod> &periods) noexcept
	{
		Date regEnd = params.termination_date();
		if (is_valid_date(params.last_regular_period_end_date()) &&
		    params.last_regular_period_end_date() < regEnd) {
			regEnd = params.last_regular_period_end_date();
		}
		if (regEnd < params.termination_date()) {
			periods.push_back(CalculationPeriod(regEnd, params.termination_date(), true));
		}
		if (params.roll_convention() != RollConvention::ROLL_CONVENTION_NONE &&
		    params.roll_convention() != RollConvention::ROLL_CONVENTION_DEFAULT) {
			if (regEnd != adjust_date(regEnd, params.roll_convention(), 0)) {
				return StatusCode::kSCH_LastRegularEndDateInconsistent;
			}
		}
		int roll_date = date_components(regEnd).d;
		Date regStart = params.effective_date();
		if (is_valid_date(params.first_regular_period_start_date()) &&
		    params.first_regular_period_start_date() > regStart) {
			regStart = params.first_regular_period_start_date();
		}
		while (regEnd > regStart) {
			Date periodStart = sub(regEnd, Period::tenor_to_period(params.calculation_frequency()));
			periodStart = adjust_date(periodStart, params.roll_convention(), roll_date);
			if (periodStart < regStart) {
				break;
			}
			periods.push_back(CalculationPeriod(periodStart, regEnd, false));
			regEnd = periodStart;
		}
		if (regEnd > params.effective_date()) {
			periods.push_back(CalculationPeriod(params.effective_date(), regEnd, true));
		}
		std::reverse(periods.begin(), periods.end());
		return StatusCode::kOk;
	}

	StatusCode validate(ScheduleParameters &params) noexcept
	{
		if (!is_valid_date(params.effective_date())) {
			return StatusCode::kSCH_EffectiveDateRequired;
		}
		if (!is_valid_date(params.termination_date()) && params.term() == TENOR_UNSPECIFIED) {
			return StatusCode::kSCH_TerminateDateOrTermRequired;
		}
		if (params.term() == TENOR_1T && !is_valid_date(params.termination_date())) {
			return StatusCode::kSCH_TerminationDateRequired;
		}
		if (!is_valid_date(params.termination_date()) && params.term() != TENOR_UNSPECIFIED &&
		    params.term() != TENOR_1T) {
			params.set_termination_date(
			    add(params.effective_date(), Period::tenor_to_period(params.term())));
		}
		if (params.calculation_frequency() == TENOR_UNSPECIFIED) {
			if (!is_zero_coupon(params) && params.payment_frequency() == TENOR_UNSPECIFIED) {
				return kSCH_CalculationFrequencyRequired;
			}
			params.set_calculation_frequency(params.payment_frequency());
		}
		return StatusCode::kOk;
	}

	static const Calendar *build_calendar(const google::protobuf::RepeatedField<google::protobuf::int32> &values)
	{
		switch (values.size()) {
		case 1:
			return get_calendar_factory()->get_calendar((BusinessCenter)values.Get(0));
		case 2:
			return get_calendar_factory()->get_calendar(
			    {(BusinessCenter)values.Get(0), (BusinessCenter)values.Get(1)});
		case 3:
			return get_calendar_factory()->get_calendar({(BusinessCenter)values.Get(0),
								     (BusinessCenter)values.Get(1),
								     (BusinessCenter)values.Get(2)});
		case 4:
			return get_calendar_factory()->get_calendar(
			    {(BusinessCenter)values.Get(0), (BusinessCenter)values.Get(1),
			     (BusinessCenter)values.Get(2), (BusinessCenter)values.Get(3)});
		default:
			return nullptr;
		}
	}

	void adjust_periods(ScheduleParameters &params, std::vector<CalculationPeriod> &periods) noexcept
	{
		const Calendar *period_calendar;
		const Calendar *payment_calendar;

		payment_calendar = build_calendar(params.payment_calendars());
		period_calendar = build_calendar(params.period_calendars());
		if (period_calendar == nullptr)
			period_calendar = payment_calendar;

		if (period_calendar != nullptr && params.period_convention() != BusinessDayConvention::UNADJUSTED &&
		    params.period_convention() != BusinessDayConvention::BUSINESS_DAY_CONVENTION_UNSPECIFIED) {
			for (CalculationPeriod &p : periods) {
				p.adjusted_start(
				    period_calendar->adjust(p.unadjusted_start(), params.period_convention()));
				p.adjusted_end(period_calendar->adjust(p.unadjusted_end(), params.period_convention()));
			}
		}

		if (payment_calendar != nullptr && params.payment_convention() != BusinessDayConvention::UNADJUSTED &&
		    params.payment_convention() != BusinessDayConvention::BUSINESS_DAY_CONVENTION_UNSPECIFIED) {
			for (CalculationPeriod &p : periods) {
				if (p.is_payment()) {
					auto paydate =
					    payment_calendar->adjust(p.adjusted_end(), params.payment_convention());
					if (params.payment_lag() != 0)
						paydate = payment_calendar->advance(
						    paydate, Period(params.payment_lag(), PeriodUnit::DAYS));
					p.adjusted_payment(paydate);
				}
			}
		}
	}

	StatusCode generate_payment_schedule(ScheduleParameters &params, std::vector<CalculationPeriod> &periods)
	{
		if (params.payment_frequency() == TENOR_1T) {
			periods.back().payment(periods.back().unadjusted_end());
			return StatusCode::kOk;
		}

		Period payFreq = Period::tenor_to_period(params.payment_frequency()).normalised();
		Period calcFreq = params.calculation_frequency() != TENOR_UNSPECIFIED
				      ? Period::tenor_to_period(params.calculation_frequency()).normalised()
				      : payFreq;

		if (payFreq.units() != calcFreq.units()) {
			return StatusCode::kSCH_IncompatiblePaymentAndCalculationFrequencies;
		}
		if (payFreq.length() < calcFreq.length()) {
			return StatusCode::kSCH_PaymentFrequencyLessThanCalculationFrequency;
		}
		if (payFreq.length() % calcFreq.length()) {
			return StatusCode::kSCH_PaymentFrequencyNotMultipleOfCalculationFrequency;
		}
		int calPeriods = payFreq.length() / calcFreq.length();
		assert(calPeriods > 0);
		if (is_valid_date(params.first_payment_date())) {
			unsigned int i = 0;
			for (; i < periods.size(); i++) {
				auto &p = periods[i];
				if (p.unadjusted_end() >= params.first_payment_date()) {
					break;
				}
			}
			for (; i < periods.size(); i += calPeriods) {
				auto &p = periods[i];
				p.payment(p.unadjusted_end());
			}
			// Last period is always a payment
			auto &p = periods.back();
			p.payment(p.unadjusted_end());
		} else {
			int i = periods.size() - 1;
			// Last period is always a payment
			auto &p = periods[i];
			p.payment(p.unadjusted_end());
			if (p.is_stub()) {
				if (is_valid_date(params.last_regular_payment_date())) {
					for (i--; i >= 0; i--) {
						auto &p = periods[i];
						if (p.unadjusted_end() <= params.last_regular_payment_date()) {
							break;
						}
					}
				} else
					i--; // skip stub
			}
			for (; i >= 0; i -= calPeriods) {
				auto &p = periods[i];
				p.payment(p.unadjusted_end());
			}
		}
		return StatusCode::kOk;
	}
};

StatusCode build_schedule(ScheduleParameters &params, std::vector<CalculationPeriod> &periods) noexcept
{
	ScheduleHelper sh;
	auto code = sh.build(params, periods);
	if (code != StatusCode::kOk)
		return code;
	if (params.payment_frequency() != TENOR_UNSPECIFIED) {
		code = sh.generate_payment_schedule(params, periods);
		if (code != StatusCode::kOk)
			return code;
	}
	sh.adjust_periods(params, periods);
	return StatusCode::kOk;
}

StatusCode build_schedule(ScheduleParameters &params, Schedule &schedule) noexcept
{
	std::vector<CalculationPeriod> periods;
	auto code = build_schedule(params, periods);
	if (code != StatusCode::kOk)
		return code;
	for (int i = 0; i < periods.size(); i++) {
		if (i == 0 && periods.size() > 1 && periods[i].is_stub()) {
			schedule.set_has_front_stub(true);
		} else if (i == periods.size() - 1 && periods.size() > 1 && periods[i].is_stub()) {
			schedule.set_has_back_stub(true);
		}
		schedule.add_adjusted_start_dates(periods[i].adjusted_start());
		schedule.add_adjusted_end_dates(periods[i].adjusted_end());
		schedule.add_adjusted_payment_dates(periods[i].is_payment() ? periods[i].adjusted_payment() : 0);
	}
	return StatusCode::kOk;
}

//////////////////////////////////// Tests

static int test1()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(1, Month::August, 2013));
	params.set_termination_date(make_date(1, Month::August, 2016));
	params.set_payment_frequency(TENOR_1T);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	if (periods.size() != 1)
		return 1;
	CalculationPeriod &p = periods[0];
	if (p.unadjusted_start() != make_date(1, Month::August, 2013) ||
	    p.unadjusted_end() != make_date(1, Month::August, 2016) || p.is_stub()) {
		return 1;
	}
	return 0;
}

static int test2()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(1, Month::August, 2013));
	params.set_term(TENOR_3Y);
	params.set_payment_frequency(TENOR_1T);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	if (params.termination_date() != make_date(1, Month::August, 2016)) {
		return 1;
	}
	if (periods.size() != 1) {
		return 1;
	}
	CalculationPeriod &p = periods[0];
	if (p.unadjusted_start() != make_date(1, Month::August, 2013) ||
	    p.unadjusted_end() != make_date(1, Month::August, 2016) || p.is_stub()) {
		return 1;
	}
	return 0;
}

static int test3()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(1, Month::August, 2013));
	params.set_termination_date(make_date(1, Month::November, 2013));
	params.set_payment_frequency(TENOR_12M);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	if (periods.size() != 1)
		return 1;
	CalculationPeriod &p = periods[0];
	if (p.unadjusted_start() != make_date(1, Month::August, 2013) ||
	    p.unadjusted_end() != make_date(1, Month::November, 2013) || p.is_stub()) {
		return 1;
	}
	return 0;
}

static int test4()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(1, Month::August, 2013));
	params.set_termination_date(make_date(1, Month::August, 2016));
	params.set_calculation_frequency(TENOR_12M);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	Date expectedStartDates[] = {make_date(1, Month::August, 2013), make_date(1, Month::August, 2014),
				     make_date(1, Month::August, 2015)};
	Date expectedEndDates[] = {make_date(1, Month::August, 2014), make_date(1, Month::August, 2015),
				   make_date(1, Month::August, 2016)};

	auto len = std::end(expectedStartDates) - std::begin(expectedStartDates);
	if (len != (int)periods.size()) {
		return 1;
	}

	for (int i = 0; i < len; i++) {
		if (periods[i].unadjusted_start() != expectedStartDates[i] ||
		    periods[i].unadjusted_end() != expectedEndDates[i]) {
			return 1;
		}
	}
	return 0;
}

static int test5()
{
	// Test auto final stub
	ScheduleParameters params;
	params.set_effective_date(make_date(1, Month::August, 2013));
	params.set_termination_date(make_date(1, Month::November, 2016));
	params.set_calculation_frequency(TENOR_12M);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	Date expectedStartDates[] = {make_date(1, Month::August, 2013), make_date(1, Month::August, 2014),
				     make_date(1, Month::August, 2015), make_date(1, Month::August, 2016)};
	Date expectedEndDates[] = {make_date(1, Month::August, 2014), make_date(1, Month::August, 2015),
				   make_date(1, Month::August, 2016), make_date(1, Month::November, 2016)};

	auto len = std::end(expectedStartDates) - std::begin(expectedStartDates);
	if (len != (int)periods.size()) {
		return 1;
	}

	for (int i = 0; i < len; i++) {
		if (periods[i].unadjusted_start() != expectedStartDates[i] ||
		    periods[i].unadjusted_end() != expectedEndDates[i]) {
			return 1;
		}
	}
	return 0;
}

static int test6()
{
	// Test initial stub
	ScheduleParameters params;
	params.set_effective_date(make_date(1, Month::August, 2013));
	params.set_termination_date(make_date(1, Month::November, 2016));
	params.set_stub_location(SHORT_FRONT_STUB);
	params.set_calculation_frequency(TENOR_12M);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	Date expectedStartDates[] = {make_date(1, Month::August, 2013), make_date(1, Month::November, 2013),
				     make_date(1, Month::November, 2014), make_date(1, Month::November, 2015)};
	Date expectedEndDates[] = {make_date(1, Month::November, 2013), make_date(1, Month::November, 2014),
				   make_date(1, Month::November, 2015), make_date(1, Month::November, 2016)};

	auto len = std::end(expectedStartDates) - std::begin(expectedStartDates);
	if (len != (int)periods.size()) {
		return 1;
	}

	for (int i = 0; i < len; i++) {
		if (periods[i].unadjusted_start() != expectedStartDates[i] ||
		    periods[i].unadjusted_end() != expectedEndDates[i]) {
			return 1;
		}
	}
	return 0;
}

static int test7()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(27, 8, 2013));
	params.set_termination_date(make_date(27, 12, 2013));
	params.set_calculation_frequency(TENOR_1M);
	params.set_payment_frequency(TENOR_2M);
	params.add_period_calendars(GBLO);
	params.add_payment_calendars(GBLO);
	params.set_period_convention(BusinessDayConvention::MODIFIED_FOLLOWING);
	params.set_payment_convention(BusinessDayConvention::MODIFIED_FOLLOWING);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	Date expectedStartDates[] = {make_date(27, 8, 2013), make_date(27, 9, 2013), make_date(27, 10, 2013),
				     make_date(27, 11, 2013)};
	Date expectedEndDates[] = {make_date(27, 9, 2013), make_date(27, 10, 2013), make_date(27, 11, 2013),
				   make_date(27, 12, 2013)};
	Date expectedPayment[] = {make_date(28, 10, 2013), make_date(27, 12, 2013)};

	auto len = std::end(expectedStartDates) - std::begin(expectedStartDates);
	if (len != (int)periods.size()) {
		return 1;
	}

	for (int i = 0; i < len; i++) {
		if (periods[i].unadjusted_start() != expectedStartDates[i] ||
		    periods[i].unadjusted_end() != expectedEndDates[i]) {
			return 1;
		}
		if (periods[i].is_payment()) {
			if (std::find(std::begin(expectedPayment), std::end(expectedPayment),
				      periods[i].adjusted_payment()) == std::end(expectedPayment)) {
				return 1;
			}
		}
	}
	return 0;
}

static int test8()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(27, 8, 2013));
	params.set_termination_date(make_date(20, 1, 2014));
	params.set_calculation_frequency(TENOR_1M);
	params.set_payment_frequency(TENOR_2M);
	params.set_stub_location(SHORT_BACK_STUB);
	params.add_period_calendars(GBLO);
	params.add_payment_calendars(GBLO);
	params.set_period_convention(BusinessDayConvention::MODIFIED_FOLLOWING);
	params.set_payment_convention(BusinessDayConvention::MODIFIED_FOLLOWING);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	Date expectedStartDates[] = {make_date(27, 8, 2013), make_date(27, 9, 2013), make_date(27, 10, 2013),
				     make_date(27, 11, 2013), make_date(27, 12, 2013)};
	Date expectedEndDates[] = {make_date(27, 9, 2013), make_date(27, 10, 2013), make_date(27, 11, 2013),
				   make_date(27, 12, 2013), make_date(20, 1, 2014)};
	Date expectedPayment[] = {make_date(28, 10, 2013), make_date(27, 12, 2013), make_date(20, 1, 2014)};

	auto len = std::end(expectedStartDates) - std::begin(expectedStartDates);
	if (len != (int)periods.size()) {
		return 1;
	}

	for (int i = 0; i < len; i++) {
		if (periods[i].unadjusted_start() != expectedStartDates[i] ||
		    periods[i].unadjusted_end() != expectedEndDates[i]) {
			return 1;
		}
		if (periods[i].is_payment()) {
			if (std::find(std::begin(expectedPayment), std::end(expectedPayment),
				      periods[i].adjusted_payment()) == std::end(expectedPayment)) {
				return 1;
			}
		}
	}
	return 0;
}

static int test9()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(27, 8, 2013));
	params.set_termination_date(make_date(20, 1, 2014));
	params.set_calculation_frequency(TENOR_1M);
	params.set_payment_frequency(TENOR_2M);
	params.add_period_calendars(GBLO);
	params.add_payment_calendars(GBLO);
	params.set_period_convention(BusinessDayConvention::MODIFIED_FOLLOWING);
	params.set_payment_convention(BusinessDayConvention::MODIFIED_FOLLOWING);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	Date expectedStartDates[] = {make_date(27, 8, 2013), make_date(27, 9, 2013), make_date(27, 10, 2013),
				     make_date(27, 11, 2013), make_date(27, 12, 2013)};
	Date expectedEndDates[] = {make_date(27, 9, 2013), make_date(27, 10, 2013), make_date(27, 11, 2013),
				   make_date(27, 12, 2013), make_date(20, 1, 2014)};
	Date expectedPayment[] = {make_date(28, 10, 2013), make_date(27, 12, 2013), make_date(20, 1, 2014)};

	auto len = std::end(expectedStartDates) - std::begin(expectedStartDates);
	if (len != (int)periods.size()) {
		return 1;
	}

	for (int i = 0; i < len; i++) {
		if (periods[i].unadjusted_start() != expectedStartDates[i] ||
		    periods[i].unadjusted_end() != expectedEndDates[i]) {
			return 1;
		}
		if (periods[i].is_payment()) {
			if (std::find(std::begin(expectedPayment), std::end(expectedPayment),
				      periods[i].adjusted_payment()) == std::end(expectedPayment)) {
				return 1;
			}
		}
	}
	return 0;
}

static int test10()
{
	ScheduleParameters params;
	params.set_effective_date(make_date(27, 8, 2013));
	params.set_termination_date(make_date(27, 12, 2013));
	params.set_calculation_frequency(TENOR_2M);
	params.set_payment_frequency(TENOR_2M);
	params.add_period_calendars(GBLO);
	params.add_payment_calendars(GBLO);
	params.set_period_convention(BusinessDayConvention::MODIFIED_FOLLOWING);
	params.set_payment_convention(BusinessDayConvention::MODIFIED_FOLLOWING);

	std::vector<CalculationPeriod> periods;

	build_schedule(params, periods);

	Date expectedStartDates[] = {make_date(27, 8, 2013), make_date(27, 10, 2013)};
	Date expectedEndDates[] = {make_date(27, 10, 2013), make_date(27, 12, 2013)};
	Date expectedPayment[] = {make_date(28, 10, 2013), make_date(27, 12, 2013)};

	auto len = std::end(expectedStartDates) - std::begin(expectedStartDates);
	if (len != (int)periods.size()) {
		return 1;
	}

	for (int i = 0; i < len; i++) {
		if (periods[i].unadjusted_start() != expectedStartDates[i] ||
		    periods[i].unadjusted_end() != expectedEndDates[i]) {
			return 1;
		}
		if (periods[i].is_payment()) {
			if (std::find(std::begin(expectedPayment), std::end(expectedPayment),
				      periods[i].adjusted_payment()) == std::end(expectedPayment)) {
				return 1;
			}
		}
	}
	return 0;
}

int test_schedule_generation()
{
	int failure_count = 0;

	failure_count += test1();
	failure_count += test2();
	failure_count += test3();
	failure_count += test4();
	failure_count += test5();
	failure_count += test6();
	failure_count += test7();
	failure_count += test8();
	failure_count += test9();
	failure_count += test10();
	if (failure_count == 0)
		printf("Schedule Generation Tests OK\n");
	else
		printf("Schedule Generation Tests FAILED\n");

	return failure_count;
}

} // namespace redukti