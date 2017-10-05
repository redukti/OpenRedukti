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
 *
 * Portions derived from Quantlib.
 * License: http://quantlib.org/license.shtml
 */
#include <date.h>
#include <logger.h>

#include <cstring>
#include <type_traits>

#ifdef _MSC_VER
#define strcasecmp _stricmp
#endif

namespace redukti
{

static_assert(std::is_standard_layout<Period>::value, "Period is not standard layout");
static_assert(std::is_standard_layout<Date>::value, "Date is not standard layout");

YearMonthDay nth_weekday(unsigned n, unsigned wd, unsigned month, int year)
{
	// Do validation here
	assert(n <= 4); // We could possibly allow n == 5, but it would sometimes fail
	assert(wd <= 6);
	assert(1 <= month && month <= 12);
	const unsigned wd_1st = weekday(make_date(1, month, year));
	const unsigned day = weekday_difference(wd, wd_1st) + 1 + (n - 1) * 7;
	return YearMonthDay{(short)year, (unsigned char)month, (unsigned char)day};
}

static Date advance(Date date, int n, PeriodUnit units) noexcept
{
	switch (units) {
	case PeriodUnit::DAYS:
		return date + n;
	case PeriodUnit::WEEKS:
		return date + 7 * n;
	case PeriodUnit::MONTHS: {
		YearMonthDay ymd = date_components(date);
		int d = ymd.d;
		int m = ymd.m + n;
		int y = ymd.y;
		while (m > 12) {
			m -= 12;
			y += 1;
		}
		while (m < 1) {
			m += 12;
			y -= 1;
		}
		int length = last_day_of_month(y, m);
		if (d > length)
			d = length;
		return make_date(d, m, y);
	}
	case PeriodUnit::YEARS: {
		YearMonthDay ymd = date_components(date);
		int d = ymd.d;
		int m = ymd.m;
		int y = ymd.y + n;

		if (d == 29 && m == (int)February && !is_leap(y))
			d = 28;
		return make_date(d, m, y);
	}
	default:
		assert(false);
		return date;
	}
}

Date add(Date date, const Period &p) noexcept { return advance(date, p.length(), p.units()); }

Date sub(Date date, const Period &p) noexcept { return advance(date, -p.length(), p.units()); }

// Parse a string representation of date
// It will detect seperator character '/' or '-'.
// The formats acceptable are yyyy/mm/dd, dd/mm/yyyy, yyyy-mm-dd, or dd-mm-yyyy
bool parse_date(const char *s, Date *d) noexcept
{
	char temp[11];
	strncpy(temp, s, sizeof temp);
	temp[sizeof temp - 1] = 0;
	const char *limit = temp + strlen(temp);

	char *start = temp;
	char *end = start;

	int year = -1;
	int month = -1;
	int day = -1;

	const char *parts[3] = {0};

	int i = 0;
	while (end < limit) {
		if (*end == '/' || *end == '-') {
			*end = 0;
			if (i < 3) {
				parts[i++] = start;
			}
			start = end + 1;
		}
		end++;
	}
	if (i < 3 && start < end) {
		parts[i++] = start;
	}
	if (i < 3) {
		debug("failed to parse date [%s]: invalid date\n", s);
		return false;
	}
	if (strlen(parts[0]) == 4) {
		year = atoi(parts[0]);
		month = atoi(parts[1]);
		day = atoi(parts[2]);
	} else {
		day = atoi(parts[0]);
		month = atoi(parts[1]);
		year = atoi(parts[2]);
	}
	*d = make_date(day, month, year);
	return true;
}

Period Period::tenor_to_period(Tenor tenor)
{
	if (tenor <= TENOR_6D)
		return Period(tenor, PeriodUnit::DAYS);
	else if (tenor <= TENOR_3W)
		return Period(tenor - TENOR_1W + 1, PeriodUnit::WEEKS);
	else if (tenor <= TENOR_23M)
		return Period(tenor - TENOR_1M + 1, PeriodUnit::MONTHS);
	else if (tenor <= TENOR_65Y)
		return Period(tenor - TENOR_2Y + 2, PeriodUnit::YEARS);
	else
		return Period(tenor, PeriodUnit::TERMS);
}

///////////////////////////////////// TESTS

static int test_consistency()
{
	Date minDate = make_date(1, 1, 1900), maxDate = make_date(30, 12, 2199);
	YearMonthDay minDateMinusOne = date_components(minDate - 1);
	int dyold = day_of_year(minDateMinusOne), dold = minDateMinusOne.d, mold = minDateMinusOne.m,
	    yold = minDateMinusOne.y, wdold = weekday(minDate - 1);

	int failure_count = 0;
	for (Date t = minDate; t <= maxDate; t++) {
		YearMonthDay t_ymd = date_components(t);
		int dy = day_of_year(t_ymd), d = t_ymd.d, m = t_ymd.m, y = t_ymd.y, wd = weekday(t);

		// check if skipping any date
		if (!((dy == dyold + 1) || (dy == 1 && dyold == 365 && !is_leap(yold)) ||
		      (dy == 1 && dyold == 366 && is_leap(yold)))) {
			failure_count++;
		}
		dyold = dy;

		if (!((d == dold + 1 && m == mold && y == yold) || (d == 1 && m == mold + 1 && y == yold) ||
		      (d == 1 && m == 1 && y == yold + 1))) {
			failure_count++;
		}
		dold = d;
		mold = m;
		yold = y;

		// check month definition
		if (m < 1 || m > 12) {
			failure_count++;
		}
		// check day definition
		if (d < 1) {
			failure_count++;
		}
		if (!((m == 1 && d <= 31) || (m == 2 && d <= 28) || (m == 2 && d == 29 && is_leap(y)) ||
		      (m == 3 && d <= 31) || (m == 4 && d <= 30) || (m == 5 && d <= 31) || (m == 6 && d <= 30) ||
		      (m == 7 && d <= 31) || (m == 8 && d <= 31) || (m == 9 && d <= 30) || (m == 10 && d <= 31) ||
		      (m == 11 && d <= 30) || (m == 12 && d <= 31))) {
			failure_count++;
		}
		// check weekday definition
		if (!((int(wd) == int(wdold + 1)) || (int(wd) == 0 && int(wdold) == 6))) {
			failure_count++;
		}
		wdold = wd;

		// create the same date with a different constructor
		Date s = make_date(d, m, y);
		// check serial number consistency
		if (s != t) {
			failure_count++;
		}
	}
	return failure_count;
}

int test_parsing()
{
	int failure_count = 0;
	Date dt0 = make_date(5, Month::Feb, 2013);
	const char *str1 = "2013/2/5";
	const char *str2 = "05/02/2013";

	Date dt1;
	parse_date(str1, &dt1);
	Date dt2;
	parse_date(str2, &dt2);

	auto wkday1 = weekday(dt1);
	if (wkday1 != Tuesday)
		failure_count++;

	auto ymd = date_components(dt2);

	if (dt0 != dt1 || dt2 != dt1) {
		failure_count++;
	}
	if (ymd.d != 5 || ymd.m != 2 || ymd.y != 2013) {
		failure_count++;
	}
	return failure_count;
}

static int test_basics()
{
	// auto ymd = date_components(1);
	// Excel returns 1 here due to bug in Excel
	// As correct answer is 2.
	// For us 1 is 31/Dec/1899
	auto d0 = make_date(1, 1, 1900);
	if (d0 <= 0)
		return 1;
	// Check Excel compatibility
	d0 = make_date(1, 3, 1900);
	if (d0 != 61)
		return 1;
	if (make_date(1, 1, 1901) != 367)
		return 1;
	if (make_date(31, 12, 2199) != 109574)
		return 1;
	auto d1 = make_date(5, 7, 2017);
	if (d1 != 42921)
		return 1;
	auto ymd1 = date_components(d1);
	if (ymd1.d != 5 || ymd1.m != 7 || ymd1.y != 2017)
		return 1;
	if (weekday(d1) != Wednesday)
		return 1;
	if (next_weekday(d1 + 1, Wednesday) != make_date(12, 7, 2017))
		return 1;
	if (next_weekday(d1, Wednesday) != d1)
		return 1;
	if (nth_weekday(1, Wednesday, 7, 2017).d != 5)
		return 1;
	if (nth_weekday(2, Wednesday, 7, 2017).d != 12)
		return 1;
	if (nth_weekday(3, Wednesday, 7, 2017).d != 19)
		return 1;
	if (nth_weekday(4, Wednesday, 7, 2017).d != 26)
		return 1;
	if (day_of_year(ymd1) != 186)
		return 1;
	return 0;
}

int test_date()
{
	int failure_count = 0;

	failure_count += test_basics();
	failure_count += test_parsing();
	failure_count += test_consistency();
	if (failure_count == 0)
		printf("Date Tests OK\n");
	else
		printf("Date Tests FAILED\n");
	return failure_count;
}

} // namespace redukti
