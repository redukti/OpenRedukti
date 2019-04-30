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
 * Portions derived from:
 * http://howardhinnant.github.io/date_algorithms.html
 * The MIT License (MIT)
 *
 * Copyright (c) 2015, 2016, 2017 Howard Hinnant
 * Copyright (c) 2016 Adrian Colomitchi
 * Copyright (c) 2017 Florian Dang
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef _REDUKTI_DATE_H_
#define _REDUKTI_DATE_H_

#include <enums.pb.h>

#include <inttypes.h>
#include <limits>

#undef min
#undef max

namespace redukti
{

// This class provides a Period (length + TimeUnit) class
// and implements a limited algebra.
// Must be standard layout for C compatibility
class Period
{
	public:
	// Construct a period from length and unit
	Period(int n, PeriodUnit unit) noexcept : length_(n), unit_(unit) {}
	// Default constructor : 0D period
	Period() noexcept : length_(0), unit_(PeriodUnit::DAYS) {}
	Period(const Period &o) = default;
	Period &operator=(const Period &o) = default;
	~Period() = default;
	int length() const noexcept { return length_; }
	PeriodUnit units() const noexcept { return unit_; }
	bool operator==(const Period &p) const noexcept { return length_ == p.length_ && unit_ == p.unit_; }
	// Normalisation converts weeks to days and
	// years to months
	Period normalised() const noexcept
	{
		if (unit_ == PeriodUnit::WEEKS) {
			return Period(length_ * 7, PeriodUnit::DAYS);
		}
		if (unit_ == PeriodUnit::YEARS) {
			return Period(length_ * 12, PeriodUnit::MONTHS);
		}
		return *this;
	}
	// Converts a tenor to period representation
	// Must be updated if definition of Tenor changes.
	static Period tenor_to_period(Tenor tenor);

	private:
	int32_t length_;
	PeriodUnit unit_;
};

enum Weekday {
	Sunday = 0,
	Monday = 1,
	Tuesday = 2,
	Wednesday = 3,
	Thursday = 4,
	Friday = 5,
	Saturday = 6,
	Sun = 0,
	Mon = 1,
	Tue = 2,
	Wed = 3,
	Thu = 4,
	Fri = 5,
	Sat = 6
};

// Month names
enum Month {
	January = 1,
	February = 2,
	March = 3,
	April = 4,
	May = 5,
	June = 6,
	July = 7,
	August = 8,
	September = 9,
	October = 10,
	November = 11,
	December = 12,
	Jan = 1,
	Feb = 2,
	Mar = 3,
	Apr = 4,
	Jun = 6,
	Jul = 7,
	Aug = 8,
	Sep = 9,
	Oct = 10,
	Nov = 11,
	Dec = 12
};

// Date type. Uses an int to
// represent a serial number.
// this implementation is immutable - hence
// thread-safe.
typedef int32_t Date;

struct YearMonthDay {
	short y;
	unsigned char m;
	unsigned char d;
};

// Returns number of days since civil 1899-12-31.  Negative values indicate
//    days prior to 1899-12-31.
// Preconditions:  y-m-d represents a date in the civil (Gregorian) calendar
//                 m is in [1, 12]
//                 d is in [1, last_day_of_month(y, m)]
//                 y is "approximately" in
//                   [numeric_limits<Int>::min()/366,
//                   numeric_limits<Int>::max()/366]
//                 Exact range of validity is:
//                 [civil_from_days(numeric_limits<Int>::min()),
//                  civil_from_days(numeric_limits<Int>::max()-719468+25569)]
// Notes: The original algorithm has been modified to make
// the serial date match Excel dates. This is done by making the start
// date 31/Dec/1899 rather than 1/Jan/1970.
constexpr Date make_date(unsigned d, unsigned m, int y) noexcept
{
	y -= m <= 2;
	const int era = (y >= 0 ? y : y - 399) / 400;
	const unsigned yoe = static_cast<unsigned>(y - era * 400);	   // [0, 399]
	const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1; // [0, 365]
	const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;	  // [0, 146096]
	return era * 146097 + static_cast<int>(doe) - 719468 + 25569; // +25569 adjusts the serial number to match Excel
}

constexpr Date make_date(YearMonthDay ymd) { return make_date(ymd.d, ymd.m, ymd.y); }

// Returns year/month/day triple in civil calendar
// Preconditions:  z is number of days since 1899-12-31 and is in the range:
//                   [numeric_limits<Int>::min(),
//                   numeric_limits<Int>::max()-719468+25569].
// Notes: The original algorithm has been modified to make
// the serial date match Excel dates. This is done by making the start
// date 31/Dec/1899 rather than 1/Jan/1970.
constexpr YearMonthDay date_components(Date z)
{
	z += (719468 - 25569);
	const int era = (z >= 0 ? z : z - 146096) / 146097;
	const unsigned doe = static_cast<unsigned>(z - era * 146097);		    // [0, 146096]
	const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; // [0, 399]
	int y = static_cast<int>(yoe) + era * 400;
	const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
	const unsigned mp = (5 * doy + 2) / 153;		      // [0, 11]
	int d = doy - (153 * mp + 2) / 5 + 1;			      // [1, 31]
	int m = mp + (mp < 10 ? 3 : -9);			      // [1, 12]
	y += (m <= 2);
	return YearMonthDay{static_cast<short>(y), static_cast<unsigned char>(m), static_cast<unsigned char>(d)};
}

// Day of the year, where Jan 1 is 1, Jan 2 is 2, Feb 1 is 32 and so on.
constexpr int day_of_year(YearMonthDay ymd) { return make_date(ymd) - make_date(31, 12, ymd.y - 1); }

// Returns day of week in civil calendar [0, 6] -> [Sun, Sat]
// Preconditions:  z is number of days since 1899-12-31 and is in the range:
//                   [numeric_limits<Int>::min(), numeric_limits<Int>::max()-4].
// Notes: The original algorithm has been modified to make
// the serial date match Excel dates. This is done by making the start
// date 31/Dec/1899 rather than 1/Jan/1970.
constexpr unsigned char weekday(Date z) noexcept
{
	z -= 25569; // adjust to 1970-01-01
	return static_cast<unsigned char>(static_cast<unsigned>(z >= -4 ? (z + 4) % 7 : (z + 5) % 7 + 6));
}

// Preconditions: m is in [1, 12]
// Returns: The number of days in the month m of common year
// The result is always in the range [28, 31].
constexpr unsigned last_day_of_month_common_year(unsigned m) noexcept
{
	constexpr unsigned char a[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	return a[m - 1];
}

// Returns: true if y is a leap year in the civil calendar, else false
constexpr bool is_leap(int y) noexcept { return y % 4 == 0 && (y % 100 != 0 || y % 400 == 0); }

// Preconditions: m is in [1, 12]
// Returns: The number of days in the month m of year y
// The result is always in the range [28, 31].
constexpr unsigned last_day_of_month(int y, unsigned m) noexcept
{
	return m != 2 || !is_leap(y) ? last_day_of_month_common_year(m) : 29u;
}

// Add/subtract periods from dates
extern Date add(Date date, const Period &) noexcept;
extern Date sub(Date date, const Period &) noexcept;

// Construct an end of month date for the
// given year and month
constexpr Date end_of_month(int y, unsigned m) noexcept { return make_date(last_day_of_month(y, m), m, y); }

// Test whether given date is the calendar end of the month
constexpr bool is_end_of_month(YearMonthDay ymd) noexcept { return ymd.d == last_day_of_month(ymd.y, ymd.m); }

// Preconditions: x <= 6 && y <= 6
// Returns: The number of days from the weekday y to the weekday x.
// The result is always in the range [0, 6].
constexpr unsigned weekday_difference(unsigned x, unsigned y) noexcept
{
	x -= y;
	return x <= 6 ? x : x + 7;
}

// Preconditions: wd <= 6
// Returns: The weekday following wd
// The result is always in the range [0, 6].
constexpr unsigned next_weekday(unsigned wd) noexcept { return wd < 6 ? wd + 1 : 0; }

// Preconditions: wd <= 6
// Returns: The weekday prior to wd
// The result is always in the range [0, 6].
inline constexpr unsigned prev_weekday(unsigned wd) noexcept { return wd > 0 ? wd - 1 : 6; }

// next given weekday following or equal to the given date
// E.g., the Friday following Tuesday, January 15th, 2002
//   was January 18th, 2002.
// see also http://www.cpearson.com/excel/DateTimeWS.htm
constexpr Date next_weekday(Date d, Weekday desired_weekday) noexcept
{
	auto yesterdays_weekday = weekday(d - 1);
	return d + weekday_difference(desired_weekday, next_weekday(yesterdays_weekday));
}

// n-th given weekday in the given month and year
// E.g., the 4th Thursday of March, 1998 was March 26th,
// 1998.
extern YearMonthDay nth_weekday(unsigned n, unsigned wd, unsigned month, int year);

constexpr bool is_weekend(unsigned wd) { return wd == Sunday || wd == Saturday; }

// Min allowed date (Jan 1st, 1901)
// This is imposed by OpenRedukti
// This is helpful because then 0 can be used to represent an invalid date
constexpr Date minimum_date() noexcept
{
	return 367; // Jan 1st, 1901
}

// We limit the max date so that we can ensure date values
// fit in 24 bits
// Dec 31st, 2199
constexpr Date maximum_date() noexcept
{
	return 109574; // Dec 31st, 2199
}

// Parse a string representation of date
// It will detect seperator character '/' or '-'.
// The formats acceptable are yyyy/mm/dd, dd/mm/yyyy, yyyy-mm-dd, or dd-mm-yyyy
// No exceptions thrown
// Returns true on success
extern bool parse_date(const char *s, Date *d) noexcept;

// We need to ensure that 0 is not a valid date as this
// helps us with protobuf representation of dates as integers
// where unspecified value is 0.
// Another requirement is to limit the max date so that
// date values can fit into 24 bits.
static inline bool is_valid_date(Date date) noexcept { return date >= minimum_date() && date <= maximum_date(); }

extern int test_date();

} // namespace redukti

#endif
