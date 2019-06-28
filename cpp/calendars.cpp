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
#include <enums.pb.h>

#include <allocators.h>
#include <calendars.h>
#include <hashtable.h>

#include <cstring>

#include <algorithm>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>

namespace redukti
{

// Western calendars only
// FIXME
bool is_weekend(Weekday w) noexcept { return w == Saturday || w == Sunday; }

class CalendarImpl : public Calendar
{
	public:
	CalendarImpl(BusinessCenter id, const std::set<Date> &holidays) noexcept : id_(id), holidays_(holidays) {}
	virtual int id() const noexcept override { return id_; }
	virtual bool is_holiday(Date d) const noexcept override
	{
		// FIXME this will only work for places where Sat/Sun are holidays
		return is_weekend(weekday(d)) || holidays_.find(d) != holidays_.end();
	}

	private:
	BusinessCenter id_;
	std::set<Date> holidays_;

	private:
	CalendarImpl(const CalendarImpl &) = delete;
	CalendarImpl &operator=(const CalendarImpl &) = delete;
};

Date Calendar::adjust(Date d, BusinessDayConvention c) const noexcept
{
	if (c == BusinessDayConvention::UNADJUSTED)
		return d;

	Date d1 = d;
	if (c == BusinessDayConvention::FOLLOWING || c == BusinessDayConvention::MODIFIED_FOLLOWING) {
		while (is_holiday(d1))
			d1 = d1 + 1;
		if (d == d1)
			return d;
		if (c == BusinessDayConvention::MODIFIED_FOLLOWING) {
			auto m1 = date_components(d1).m;
			auto m = date_components(d).m;
			if (m1 != m) {
				return adjust(d, BusinessDayConvention::PRECEDING);
			}
		}
	} else if (c == BusinessDayConvention::PRECEDING || c == BusinessDayConvention::MODIFIED_PRECEDING) {
		while (is_holiday(d1))
			d1 = d1 - 1;
		if (d1 == d)
			return d;
		auto m1 = date_components(d1).m;
		auto m = date_components(d).m;
		if (c == BusinessDayConvention::MODIFIED_PRECEDING && m1 != m) {
			return adjust(d, BusinessDayConvention::FOLLOWING);
		}
	} else {
		assert(false);
	}
	return d1;
}

Date Calendar::advance(Date d, int n, PeriodUnit unit, BusinessDayConvention c, bool endOfMonth) const noexcept
{
	if (n == 0) {
		return adjust(d, c);
	} else if (unit == PeriodUnit::DAYS) {
		Date d1 = d;
		if (n > 0) {
			while (n > 0) {
				d1 = d1 + 1;
				while (is_holiday(d1))
					d1 = d1 + 1;
				n--;
			}
		} else {
			while (n < 0) {
				d1 = d1 - 1;
				while (is_holiday(d1))
					d1 = d1 - 1;
				n++;
			}
		}
		return d1;
	} else if (unit == PeriodUnit::WEEKS) {
		Date d1 = add(d, Period(n, unit));
		return adjust(d1, c);
	} else {
		Date d1 = add(d, Period(n, unit));
		// we are sure the unit is Months or Years
		if (endOfMonth && is_end_of_month(d))
			return Calendar::end_of_month(d1);
		return adjust(d1, c);
	}
}

bool Calendar::is_end_of_month(Date d) const noexcept
{
	auto m1 = date_components(d).m;
	auto m2 = date_components(adjust(d + 1)).m;
	return m1 != m2;
}

Date Calendar::advance(Date d, const Period &p, BusinessDayConvention c, bool endOfMonth) const noexcept
{
	return advance(d, p.length(), p.units(), c, endOfMonth);
}

int Calendar::business_days_between(Date from, Date to, bool includeFirst, bool includeLast) const noexcept
{
	int wd = 0;
	if (from != to) {
		if (from < to) {
			// the last one is treated separately to avoid
			// incrementing Date::maxDate()
			for (Date d = from; d < to; d += 1) {
				if (is_businessday(d))
					++wd;
			}
			if (is_businessday(to))
				++wd;
		} else if (from > to) {
			for (Date d = to; d < from; ++d) {
				if (is_businessday(d))
					++wd;
			}
			if (is_businessday(from))
				++wd;
		}
		if (is_businessday(from) && !includeFirst)
			wd--;
		if (is_businessday(to) && !includeLast)
			wd--;
		if (from > to)
			wd = -wd;
	}
	return wd;
}

// Western calendars

int easter_monday(int y) noexcept
{
	// clang-format off
	static const int EasterMonday[] = {
	    98,  90,  103, 95,  114, 106, 91,  111, 102,      // 1901-1909
	    87,  107, 99,  83,  103, 95,  115, 99,  91,  111, // 1910-1919
	    96,  87,  107, 92,  112, 103, 95,  108, 100, 91,  // 1920-1929
	    111, 96,  88,  107, 92,  112, 104, 88,  108, 100, // 1930-1939
	    85,  104, 96,  116, 101, 92,  112, 97,  89,  108, // 1940-1949
	    100, 85,  105, 96,  109, 101, 93,  112, 97,  89,  // 1950-1959
	    109, 93,  113, 105, 90,  109, 101, 86,  106, 97,  // 1960-1969
	    89,  102, 94,  113, 105, 90,  110, 101, 86,  106, // 1970-1979
	    98,  110, 102, 94,  114, 98,  90,  110, 95,  86,  // 1980-1989
	    106, 91,  111, 102, 94,  107, 99,  90,  103, 95,  // 1990-1999
	    115, 106, 91,  111, 103, 87,  107, 99,  84,  103, // 2000-2009
	    95,  115, 100, 91,  111, 96,  88,  107, 92,  112, // 2010-2019
	    104, 95,  108, 100, 92,  111, 96,  88,  108, 92,  // 2020-2029
	    112, 104, 89,  108, 100, 85,  105, 96,  116, 101, // 2030-2039
	    93,  112, 97,  89,  109, 100, 85,  105, 97,  109, // 2040-2049
	    101, 93,  113, 97,  89,  109, 94,  113, 105, 90,  // 2050-2059
	    110, 101, 86,  106, 98,  89,  102, 94,  114, 105, // 2060-2069
	    90,  110, 102, 86,  106, 98,  111, 102, 94,  114, // 2070-2079
	    99,  90,  110, 95,  87,  106, 91,  111, 103, 94,  // 2080-2089
	    107, 99,  91,  103, 95,  115, 107, 91,  111, 103, // 2090-2099
	    88,  108, 100, 85,  105, 96,  109, 101, 93,  112, // 2100-2109
	    97,  89,  109, 93,  113, 105, 90,  109, 101, 86,  // 2110-2119
	    106, 97,  89,  102, 94,  113, 105, 90,  110, 101, // 2120-2129
	    86,  106, 98,  110, 102, 94,  114, 98,  90,  110, // 2130-2139
	    95,  86,  106, 91,  111, 102, 94,  107, 99,  90,  // 2140-2149
	    103, 95,  115, 106, 91,  111, 103, 87,  107, 99,  // 2150-2159
	    84,  103, 95,  115, 100, 91,  111, 96,  88,  107, // 2160-2169
	    92,  112, 104, 95,  108, 100, 92,  111, 96,  88,  // 2170-2179
	    108, 92,  112, 104, 89,  108, 100, 85,  105, 96,  // 2180-2189
	    116, 101, 93,  112, 97,  89,  109, 100, 85,  105  // 2190-2199
	};
	// clang-format on
	return EasterMonday[y - 1901];
}

class GBLOCalendarImpl : public Calendar
{
	public:
	GBLOCalendarImpl() noexcept {}
	virtual int id() const noexcept override { return BusinessCenter::GBLO; }
	virtual bool is_holiday(Date date) const noexcept override
	{
		auto w = weekday(date);
		YearMonthDay ymd = date_components(date);
		int dd = day_of_year(ymd);
		int d = ymd.d, m = ymd.m, y = ymd.y;
		int em = easter_monday(y);
		// clang-format off
		if (is_weekend(w)
		    // New Year's Day (possibly moved to Monday)
		    || ((d == 1 || ((d == 2 || d == 3) && w == Monday)) && m == January)
		    // Good Friday
		    || (dd == em - 3)
		    // Easter Monday
		    || (dd == em)
		    // first Monday of May (Early May Bank Holiday)
		    || (d <= 7 && w == Monday && m == May)
		    // last Monday of May (Spring Bank Holiday)
		    || (d >= 25 && w == Monday && m == May && y != 2002)
		    // last Monday of August (Summer Bank Holiday)
		    || (d >= 25 && w == Monday && m == August)
		    // Christmas (possibly moved to Monday or Tuesday)
		    || ((d == 25 || (d == 27 && (w == Monday || w == Tuesday))) && m == December)
		    // Boxing Day (possibly moved to Monday or Tuesday)
		    || ((d == 26 || (d == 28 && (w == Monday || w == Tuesday))) && m == December)
		    // June 3rd, 2002 only (Golden Jubilee Bank Holiday)
		    // June 4rd, 2002 only (special Spring Bank Holiday)
		    || ((d == 3 || d == 4) && m == June && y == 2002)
		    // December 31st, 1999 only
		    || (d == 31 && m == December && y == 1999))
			return true;
		// clang-format on
		return false;
	}

	private:
	GBLOCalendarImpl(const GBLOCalendarImpl &);
	GBLOCalendarImpl &operator=(const GBLOCalendarImpl &);
};

class USNYCalendarImpl : public Calendar
{
	public:
	USNYCalendarImpl() noexcept {}
	virtual int id() const noexcept override { return BusinessCenter::USNY; }
	virtual bool is_holiday(Date date) const noexcept override
	{
		auto w = weekday(date);
		YearMonthDay ymd = date_components(date);
		int dd = day_of_year(ymd);
		int d = ymd.d, m = ymd.m, y = ymd.y;
		int em = easter_monday(y);
		// clang-format off
		if (is_weekend(w)
		    // New Year's Day (possibly moved to Monday if on Sunday)
		    || ((d == 1 || (d == 2 && w == Monday)) && m == January)
		    // Martin Luther King's birthday (third Monday in January)
		    || ((d >= 15 && d <= 21) && w == Monday && m == January)
		    // Washington's birthday (third Monday in February)
		    || ((d >= 15 && d <= 21) && w == Monday && m == February)
		    // Good Friday
		    || (dd == em - 3)
		    // Memorial Day (last Monday in May)
		    || (d >= 25 && w == Monday && m == May)
		    // Independence Day (Monday if Sunday or Friday if Saturday)
		    || ((d == 4 || (d == 5 && w == Monday) || (d == 3 && w == Friday)) && m == July)
		    // Labor Day (first Monday in September)
		    || (d <= 7 && w == Monday && m == September)
		    // Columbus Day (second Monday in October)
		    || ((d >= 8 && d <= 14) && w == Monday && m == October)
		    // Veteran's Day (Monday if Sunday or Friday if Saturday)
		    || ((d == 11 || (d == 12 && w == Monday) || (d == 10 && w == Friday)) && m == November)
		    // Thanksgiving Day (fourth Thursday in November)
		    || ((d >= 22 && d <= 28) && w == Thursday && m == November)
		    // Christmas (Monday if Sunday or Friday if Saturday)
		    || ((d == 25 || (d == 26 && w == Monday) || (d == 24 && w == Friday)) && m == December))
			return true;
		// clang-format on
		return false;
	}

	private:
	USNYCalendarImpl(const USNYCalendarImpl &);
	USNYCalendarImpl &operator=(const USNYCalendarImpl &);
};

class EUTACalendarImpl : public Calendar
{
	public:
	EUTACalendarImpl() noexcept {}
	virtual int id() const noexcept override { return BusinessCenter::EUTA; }
	virtual bool is_holiday(Date date) const noexcept override
	{
		auto w = weekday(date);
		YearMonthDay ymd = date_components(date);
		int dd = day_of_year(ymd);
		int d = ymd.d, m = ymd.m, y = ymd.y;
		int em = easter_monday(y);
		// clang-format off
		if (is_weekend(w)
		    // New Year's Day
		    || (d == 1 && m == January)
		    // Good Friday
		    || (dd == em - 3 && y >= 2000)
		    // Easter Monday
		    || (dd == em && y >= 2000)
		    // Labour Day
		    || (d == 1 && m == May && y >= 2000)
		    // Christmas
		    || (d == 25 && m == December)
		    // Day of Goodwill
		    || (d == 26 && m == December && y >= 2000)
		    // December 31st, 1998, 1999, and 2001 only
		    || (d == 31 && m == December && (y == 1998 || y == 1999 || y == 2001)))
			return true;
		// clang-format on
		return false;
	}

	private:
	EUTACalendarImpl(const EUTACalendarImpl &);
	EUTACalendarImpl &operator=(const EUTACalendarImpl &);
};

/*
    Holidays for the Bovespa stock exchange
    <ul>
    <li>Saturdays</li>
    <li>Sundays</li>
    <li>New Year's Day, January 1st</li>
    <li>Sao Paulo City Day, January 25th</li>
    <li>Tiradentes's Day, April 21th</li>
    <li>Labour Day, May 1st</li>
    <li>Revolution Day, July 9th</li>
    <li>Independence Day, September 7th</li>
    <li>Nossa Sra. Aparecida Day, October 12th</li>
    <li>All Souls Day, November 2nd</li>
    <li>Republic Day, November 15th</li>
    <li>Black Consciousness Day, November 20th (since 2007)</li>
    <li>Christmas, December 25th</li>
    <li>Passion of Christ</li>
    <li>Carnival</li>
    <li>Corpus Christi</li>
    <li>the last business day of the year</li>
    </ul>
    */
class BRSPCalendarImpl : public Calendar
{
	public:
	BRSPCalendarImpl() noexcept {}
	virtual int id() const noexcept override { return BusinessCenter::BRSP; }
	virtual bool is_holiday(Date date) const noexcept override
	{
		auto w = weekday(date);
		YearMonthDay ymd = date_components(date);
		int dd = day_of_year(ymd);
		int d = ymd.d, m = ymd.m, y = ymd.y;
		int em = easter_monday(y);
		// clang-format off
		if (is_weekend(w)
		    // New Year's Day
		    || (d == 1 && m == January)
		    // Sao Paulo City Day
		    //|| (d == 25 && m == January)
		    // Tiradentes Day
		    || (d == 21 && m == April)
		    // Labor Day
		    || (d == 1 && m == May)
		    // Revolution Day
		    //|| (d == 9 && m == July)
		    // Independence Day
		    || (d == 7 && m == September)
		    // Nossa Sra. Aparecida Day
		    || (d == 12 && m == October)
		    // All Souls Day
		    || (d == 2 && m == November)
		    // Republic Day
		    || (d == 15 && m == November)
		    // Black Consciousness Day
		    //|| (d == 20 && m == November && y >= 2007)
		    // Christmas
		    || (d == 25 && m == December)
		    // Passion of Christ
		    || (dd == em - 3)
		    // Carnival
		    || (dd == em - 49 || dd == em - 48)
		    // Corpus Christi
		    || (dd == em + 59)
		    // last business day of the year
		    //|| (m == December && (d == 31 || (d >= 29 && w == Friday)))
		)
			return true;
		// clang-format on
		return false;
	}

	private:
	BRSPCalendarImpl(const BRSPCalendarImpl &);
	BRSPCalendarImpl &operator=(const BRSPCalendarImpl &);
};

class AUSYCalendarImpl : public Calendar
{
	public:
	AUSYCalendarImpl() noexcept {}
	virtual int id() const noexcept override { return BusinessCenter::AUSY; }
	virtual bool is_holiday(Date date) const noexcept override
	{
		auto w = weekday(date);
		YearMonthDay ymd = date_components(date);
		int dd = day_of_year(ymd);
		int d = ymd.d, m = ymd.m, y = ymd.y;
		int em = easter_monday(y);
		// clang-format off
		if (is_weekend(w)
		    // New Year's Day (possibly moved to Monday)
		    || (d == 1 && m == January)
		    // Australia Day, January 26th (possibly moved to Monday)
		    || ((d == 26 || ((d == 27 || d == 28) && w == Monday)) && m == January)
		    // Good Friday
		    || (dd == em - 3)
		    // Easter Monday
		    || (dd == em)
		    // ANZAC Day, April 25th (possibly moved to Monday)
		    || ((d == 25 || (d == 26 && w == Monday)) && m == April)
		    // Queen's Birthday, second Monday in June
		    || ((d > 7 && d <= 14) && w == Monday && m == June)
		    // Bank Holiday, first Monday in August
		    || (d <= 7 && w == Monday && m == August)
		    // Labour Day, first Monday in October
		    || (d <= 7 && w == Monday && m == October)
		    // Christmas, December 25th (possibly Monday or Tuesday)
		    || ((d == 25 || (d == 27 && (w == Monday || w == Tuesday))) && m == December)
		    // Boxing Day, December 26th (possibly Monday or Tuesday)
		    || ((d == 26 || (d == 28 && (w == Monday || w == Tuesday))) && m == December))
			return true;
		// clang-format on
		return false;
	}

	private:
	AUSYCalendarImpl(const BRSPCalendarImpl &);
	AUSYCalendarImpl &operator=(const AUSYCalendarImpl &);
};

class JPTOCalendarImpl : public Calendar
{
	public:
	JPTOCalendarImpl() noexcept {}
	virtual int id() const noexcept override { return BusinessCenter::JPTO; }
	virtual bool is_holiday(Date date) const noexcept override
	{
		auto w = weekday(date);
		YearMonthDay ymd = date_components(date);
		// int dd = day_of_year(ymd);
		int d = ymd.d, m = ymd.m, y = ymd.y;
		// equinox calculation
		const double exact_vernal_equinox_time = 20.69115;
		const double exact_autumnal_equinox_time = 23.09;
		const double diff_per_year = 0.242194;
		const double moving_amount = (y - 2000) * diff_per_year;
		int number_of_leap_years = (y - 2000) / 4 + (y - 2000) / 100 - (y - 2000) / 400;
		int ve = // vernal equinox day
		    int(exact_vernal_equinox_time + moving_amount - number_of_leap_years);
		int ae = // autumnal equinox day
		    int(exact_autumnal_equinox_time + moving_amount - number_of_leap_years);
		// checks
		// clang-format off
		if (is_weekend(w)
		    // New Year's Day
		    || (d == 1 && m == January)
		    // Bank Holiday
		    || (d == 2 && m == January)
		    // Bank Holiday
		    || (d == 3 && m == January)
		    // Coming of Age Day (2nd Monday in January),
		    // was January 15th until 2000
		    || (w == Monday && (d >= 8 && d <= 14) && m == January && y >= 2000) ||
		    ((d == 15 || (d == 16 && w == Monday)) && m == January && y < 2000)
		    // National Foundation Day
		    || ((d == 11 || (d == 12 && w == Monday)) && m == February)
		    // Vernal Equinox
		    || ((d == ve || (d == ve + 1 && w == Monday)) && m == March)
		    // Greenery Day
		    || ((d == 29 || (d == 30 && w == Monday)) && m == April)
		    // Constitution Memorial Day
		    || (d == 3 && m == May)
		    // Holiday for a Nation
		    || (d == 4 && m == May)
		    // Children's Day
		    || (d == 5 && m == May)
		    // any of the three above observed later if on Saturday or Sunday
		    || (d == 6 && m == May && (w == Monday || w == Tuesday || w == Wednesday))
		    // Marine Day (3rd Monday in July),
		    // was July 20th until 2003, not a holiday before 1996
		    || (w == Monday && (d >= 15 && d <= 21) && m == July && y >= 2003) ||
		    ((d == 20 || (d == 21 && w == Monday)) && m == July && y >= 1996 && y < 2003)
		    // Respect for the Aged Day (3rd Monday in September),
		    // was September 15th until 2003
		    || (w == Monday && (d >= 15 && d <= 21) && m == September && y >= 2003) ||
		    ((d == 15 || (d == 16 && w == Monday)) && m == September && y < 2003)
		    // If a single day falls between Respect for the Aged Day
		    // and the Autumnal Equinox, it is holiday
		    || (w == Tuesday && d + 1 == ae && d >= 16 && d <= 22 && m == September && y >= 2003)
		    // Autumnal Equinox
		    || ((d == ae || (d == ae + 1 && w == Monday)) && m == September)
		    // Health and Sports Day (2nd Monday in October),
		    // was October 10th until 2000
		    || (w == Monday && (d >= 8 && d <= 14) && m == October && y >= 2000) ||
		    ((d == 10 || (d == 11 && w == Monday)) && m == October && y < 2000)
		    // National Culture Day
		    || ((d == 3 || (d == 4 && w == Monday)) && m == November)
		    // Labor Thanksgiving Day
		    || ((d == 23 || (d == 24 && w == Monday)) && m == November)
		    // Emperor's Birthday
		    || ((d == 23 || (d == 24 && w == Monday)) && m == December && y >= 1989)
		    // Bank Holiday
		    || (d == 31 && m == December)
		    // one-shot holidays
		    // Marriage of Prince Akihito
		    || (d == 10 && m == April && y == 1959)
		    // Rites of Imperial Funeral
		    || (d == 24 && m == February && y == 1989)
		    // Enthronement Ceremony
		    || (d == 12 && m == November && y == 1990)
		    // Marriage of Prince Naruhito
		    || (d == 9 && m == June && y == 1993))
			return true;
		// clang-format on
		return false;
	}

	private:
	JPTOCalendarImpl(const JPTOCalendarImpl &);
	JPTOCalendarImpl &operator=(const JPTOCalendarImpl &);
};

static GBLOCalendarImpl gblo;
static USNYCalendarImpl usny;
static EUTACalendarImpl euta;
static BRSPCalendarImpl brsp;
static AUSYCalendarImpl ausy;
static JPTOCalendarImpl jpto;

/**
 * Make a joint calendar from the supplied
 * calendars
 */
Calendar *make_joint_calendar(int id, JointCalendarRule rule, const Calendar *c1, const Calendar *c2,
			      const Calendar *c3 = nullptr, const Calendar *c4 = nullptr) noexcept;

template <> struct NilHelper<int> {
	static void set_nil(int *v) { *v = 0; }
	static bool is_nil(const int *v) { return *v == 0; }
};

template <> struct NilHelper<Calendar *> {
	static void set_nil(Calendar **v) { *v = nullptr; }
	static bool is_nil(const Calendar *const *v) { return nullptr == *v; }
};

class CalendarFactoryImpl : public CalendarService
{
	public:
	CalendarFactoryImpl() : lock_(), inuse_(false)
	{
		std::fill(std::begin(default_calendars_), std::end(default_calendars_), nullptr);
		std::fill(std::begin(default_calendars_allocstatus_), std::end(default_calendars_allocstatus_), false);
		default_calendars_[AUSY] = &ausy;
		default_calendars_[EUTA] = &euta;
		default_calendars_[GBLO] = &gblo;
		default_calendars_[USNY] = &usny;
		default_calendars_[JPTO] = &jpto;
		default_calendars_[BRSP] = &brsp;
	}
	virtual ~CalendarFactoryImpl() noexcept;
	virtual const Calendar *get_calendar(BusinessCenter id) noexcept
	{
		std::lock_guard<std::mutex> guard(lock_);
		inuse_ = true;
		if (id == BusinessCenter::BUSINESS_CENTER_UNSPECIFIED || !BusinessCenter_IsValid(id))
			return nullptr;
		return default_calendars_[id];
	}
	// Set a calendar to given instance.
	// The service will take ownership of the instance
	// May fail if calendar instance already set and has been
	// accessed by a client - i.e. new calendars can only be set prior to
	// any use.
	virtual bool set_calendar(BusinessCenter id, std::unique_ptr<Calendar> calendar) noexcept
	{
		std::lock_guard<std::mutex> guard(lock_);
		if (id == BusinessCenter::BUSINESS_CENTER_UNSPECIFIED || !BusinessCenter_IsValid(id))
			return false;
		if (default_calendars_[id] && inuse_)
			return false;
		default_calendars_[id] = calendar.release();
		default_calendars_allocstatus_[id] = true;
		return true;
	}

	// Create a calendar from a set of holidays and assign it to the business center
	// If the assignment is successful the service will take ownership of the instance
	// May fail if calendar instance already set and has been
	// accessed by a client - i.e. new calendars can only be set prior to
	// any use.
	virtual bool set_calendar(BusinessCenter id, const Date *holidays, size_t n) noexcept final
	{
		// This is kind of wasteful but we don't care as it is rare event
		std::set<Date> holiday_set;
		for (int i = 0; i < n; i++)
			holiday_set.insert(holidays[i]);
		auto calendar = std::make_unique<CalendarImpl>(id, holiday_set);
		return set_calendar(id, std::move(calendar));
	}

	virtual Calendar *get_calendar(JointCalendarParameters calendars, JointCalendarRule rule) noexcept;

	private:
	CalendarFactoryImpl(const CalendarFactoryImpl &) = delete;
	CalendarFactoryImpl &operator=(const CalendarFactoryImpl &) = delete;

	private:
	std::mutex lock_;
	std::array<Calendar *, BusinessCenter_ARRAYSIZE> default_calendars_;
	std::array<bool, BusinessCenter_ARRAYSIZE> default_calendars_allocstatus_;
	Table<int, Calendar *> joint_calendars_cache;
	bool inuse_;
};

CalendarFactoryImpl::~CalendarFactoryImpl() noexcept
{
	for (auto iter = std::begin(joint_calendars_cache); iter != std::end(joint_calendars_cache); iter++) {
		Calendar *c = iter->value();
		if (c != nullptr) {
			// printf("Deleting joint calendar %d\n", c->id());
			delete c;
		}
	}
	for (int i = 0; i < default_calendars_.size(); i++) {
		if (default_calendars_allocstatus_[i])
			delete default_calendars_[i];
	}
}

Calendar *CalendarFactoryImpl::get_calendar(JointCalendarParameters calendars, JointCalendarRule rule) noexcept
{
	// sort the calendars
	std::sort(std::begin(calendars.centers), std::end(calendars.centers));
	// create a concatenated id
	int id = ((calendars.centers[3] << 1) | rule) | (calendars.centers[2] << 8) | (calendars.centers[1] << 16) |
		 (calendars.centers[0] << 24);

	{
		std::lock_guard<std::mutex> guard(lock_);
		inuse_ = true;
		auto val = joint_calendars_cache.get(id);
		if (val)
			return val->val();
	}
	std::array<const Calendar *, 4> raw;
	std::fill(std::begin(raw), std::end(raw), nullptr);
	int x = 0;
	for (int i = 0; i < calendars.centers.size(); i++) {
		if (calendars.centers[i] == BusinessCenter::BUSINESS_CENTER_UNSPECIFIED)
			continue;
		const Calendar *c = get_calendar(calendars.centers[i]);
		if (c == nullptr) {
			// FIXME generate error
			assert(false);
		}
		raw[x++] = c;
	}
	auto joint_calendar = make_joint_calendar(id, rule, raw[0], raw[1], raw[2], raw[3]);
	{
		std::lock_guard<std::mutex> guard(lock_);
		// Maybe another thread got there first?
		auto val = joint_calendars_cache.get(id);
		if (val)
			return val->val();
		val = joint_calendars_cache.newkey(id);
		val->setval(joint_calendar);
		// printf("Added joint calendar %d\n", val->val()->id());
		return joint_calendar;
	}
}

static std::unique_ptr<CalendarService> defaultCalendarFactory(new CalendarFactoryImpl());

CalendarService *get_calendar_factory() noexcept { return defaultCalendarFactory.get(); }

/**
 * Composite calendar
 */
class CalendarSet : public Calendar
{
	public:
	CalendarSet(int id, JointCalendarRule rule, const Calendar *c1, const Calendar *c2,
		    const Calendar *c3 = nullptr, const Calendar *c4 = nullptr) noexcept
	    : id_(id), rule_(rule)
	{
		calendars_[0] = c1;
		calendars_[1] = c2;
		calendars_[2] = c3;
		calendars_[3] = c4;
	}
	virtual ~CalendarSet() noexcept {}
	/**
	 * Check if the given date is a holiday. All
	 * holiday calendars checked.
	 */
	virtual bool is_holiday(Date d) const noexcept override;
	virtual int id() const noexcept override { return id_; }
	virtual void get_ids(std::array<BusinessCenter, 4> &centers) const noexcept override
	{
		std::fill(std::begin(centers), std::end(centers), BusinessCenter::BUSINESS_CENTER_UNSPECIFIED);
		auto n = std::end(calendars_) - std::begin(calendars_);
		for (int i = 0; i < n && calendars_[i] != nullptr; i++) {
			centers[i] = (BusinessCenter)calendars_[i]->id();
		}
	}

	private:
	int id_;
	const Calendar *calendars_[4];
	JointCalendarRule rule_;

	private:
	CalendarSet(const CalendarSet &) = delete;
	CalendarSet &operator=(const CalendarSet &) = delete;
};

bool CalendarSet::is_holiday(Date date) const noexcept
{
	auto n = std::end(calendars_) - std::begin(calendars_);
	switch (rule_) {
	case JOIN_BUSINESS_DAYS:
		/**
		 * A date is a business day
		 * for the joint calendar
		 * if it is a business day
		 * for any of the given
		 * calendars
		 */
		for (int i = 0; i < n && calendars_[i] != nullptr; i++) {
			if (calendars_[i]->is_businessday(date))
				return false;
		}
		return true;
	case JOIN_HOLIDAYS:
	default:
		/**
		 * A date is a holiday
		 * for the joint calendar
		 * if it is a holiday
		 * for any of the given
		 * calendars
		 */
		for (int i = 0; i < n && calendars_[i] != nullptr; i++) {
			if (calendars_[i]->is_holiday(date))
				return true;
		}
		return false;
	}
}

Calendar *make_joint_calendar(int id, JointCalendarRule rule, const Calendar *c1, const Calendar *c2,
			      const Calendar *c3, const Calendar *c4) noexcept
{
	return new CalendarSet(id, rule, c1, c2, c3, c4);
}

const Calendar *build_calendar(CalendarService *calendar_service,
			       const google::protobuf::RepeatedField<google::protobuf::int32> &values,
			       JointCalendarRule rule)
{
	switch (values.size()) {
	case 1:
		return calendar_service->get_calendar((BusinessCenter)values.Get(0));
	case 2:
		return calendar_service->get_calendar({(BusinessCenter)values.Get(0), (BusinessCenter)values.Get(1)},
						      rule);
	case 3:
		return calendar_service->get_calendar(
		    {(BusinessCenter)values.Get(0), (BusinessCenter)values.Get(1), (BusinessCenter)values.Get(2)},
		    rule);
	case 4:
		return calendar_service->get_calendar({(BusinessCenter)values.Get(0), (BusinessCenter)values.Get(1),
						       (BusinessCenter)values.Get(2), (BusinessCenter)values.Get(3)},
						      rule);
	default:
		return nullptr;
	}
}

const Calendar *build_calendar(CalendarService *calendar_service, const std::vector<BusinessCenter> &values,
			       JointCalendarRule rule)
{
	switch (values.size()) {
	case 1:
		return calendar_service->get_calendar(values[0]);
	case 2:
		return calendar_service->get_calendar({values[0], values[1]}, rule);
	case 3:
		return calendar_service->get_calendar({values[0], values[1], values[2]}, rule);
	case 4:
		return calendar_service->get_calendar({values[0], values[1], values[2], values[3]}, rule);
	default:
		return nullptr;
	}
}

int test_calendars()
{
	int failure_count = 0;

	auto factory = get_calendar_factory();
	auto USNY = factory->get_calendar(BusinessCenter::USNY);
	auto GBLO = factory->get_calendar(BusinessCenter::GBLO);
	auto USNY_GBLO = factory->get_calendar({BusinessCenter::USNY, BusinessCenter::GBLO}, JOIN_HOLIDAYS);
	auto USNY_GBLO2 = factory->get_calendar({BusinessCenter::GBLO, BusinessCenter::USNY}, JOIN_HOLIDAYS);
	if (USNY_GBLO != USNY_GBLO2)
		failure_count++;
	auto USNY_GBLO3 = factory->get_calendar({BusinessCenter::GBLO, BusinessCenter::USNY}, JOIN_BUSINESS_DAYS);
	if (USNY_GBLO == USNY_GBLO3)
		failure_count++;
	failure_count += USNY->is_holiday(make_date(1, Month::January, 2013)) ? 0 : 1;
	failure_count += USNY->is_businessday(make_date(2, Month::January, 2013)) ? 0 : 1;
	failure_count += USNY_GBLO->is_holiday(make_date(21, Month::January, 2013)) ? 0 : 1;
	failure_count += USNY->is_holiday(make_date(21, Month::January, 2013)) ? 0 : 1;
	failure_count += GBLO->is_holiday(make_date(21, Month::January, 2013)) ? 1 : 0;
	failure_count += USNY->is_holiday(make_date(18, Month::February, 2013)) ? 0 : 1;
	failure_count += GBLO->is_holiday(make_date(18, Month::February, 2013)) ? 1 : 0;
	failure_count += USNY_GBLO->is_holiday(make_date(18, Month::February, 2013)) ? 0 : 1;
	failure_count += USNY->is_holiday(make_date(29, Month::March, 2013)) ? 0 : 1;
	failure_count += USNY->is_holiday(make_date(1, Month::April, 2013)) ? 1 : 0;
	failure_count += USNY_GBLO->is_holiday(make_date(4, Month::July, 2013)) ? 0 : 1;
	failure_count += USNY_GBLO->is_holiday(make_date(2, Month::September, 2013)) ? 0 : 1;
	failure_count += USNY_GBLO->is_holiday(make_date(14, Month::October, 2013)) ? 0 : 1;
	failure_count += USNY_GBLO->is_holiday(make_date(11, Month::November, 2013)) ? 0 : 1;
	failure_count += USNY_GBLO->is_holiday(make_date(28, Month::November, 2013)) ? 0 : 1;
	if (failure_count == 0)
		printf("Calendar Tests OK\n");
	else
		printf("Calendar Tests FAILED\n");
	return failure_count;
}
} // namespace redukti
