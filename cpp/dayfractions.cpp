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
#include <dayfractions.h>

#include <algorithm>
#include <map>
#include <mutex>
#include <stdexcept>

#include <cmath>
#include <cstring>

namespace redukti
{

class Actual360 : public DayFraction
{
	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::ACT_360; }
	virtual double year_fraction(Date d1, Date d2) const override final { return (d2 - d1) / 360.0; }
	virtual double year_fraction(Date d1, Date d2, bool eom) const override final { return year_fraction(d1, d2); }
	virtual double year_fraction(Date d1, Date d2, Date periodStart, Date periodEnd) const override final
	{
		return year_fraction(d1, d2);
	}
};

class Actual365Fixed : public DayFraction
{
	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::ACT_365_FIXED; }
	virtual double year_fraction(Date d1, Date d2) const override final { return (d2 - d1) / 365.0; }
	virtual double year_fraction(Date d1, Date d2, bool eom) const override final { return year_fraction(d1, d2); }
	virtual double year_fraction(Date d1, Date d2, Date periodStart, Date periodEnd) const override final
	{
		return year_fraction(d1, d2);
	}
};

/**
 * 30/360
 */
class Thirty360_ : public DayFraction
{
	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::THIRTY_360; }
	virtual double year_fraction(Date d1, Date d2) const override final
	{
		//        int dd1 = d1.day_of_month(), dd2 = d2.day_of_month();
		//        int mm1 = d1.month(), mm2 = d2.month();
		//        int yy1 = d1.year(), yy2 = d2.year();
		auto ymd1 = date_components(d1);
		int dd1 = ymd1.d, mm1 = ymd1.m, yy1 = ymd1.y;
		auto ymd2 = date_components(d2);
		int dd2 = ymd2.d, mm2 = ymd2.m, yy2 = ymd2.y;
		if (dd2 == 31 && dd1 < 30) {
			dd2 = 1;
			mm2++;
		}
		return (360 * (yy2 - yy1) + 30 * (mm2 - mm1 - 1) + std::max(int(0), 30 - dd1) +
			std::min(int(30), dd2)) /
		       360.0;
	}
	virtual double year_fraction(Date d1, Date d2, bool eom) const override final { return year_fraction(d1, d2); }
	virtual double year_fraction(Date d1, Date d2, Date periodStart, Date periodEnd) const override final
	{
		return year_fraction(d1, d2);
	}
};

/**
 * 30E/360
 */
class Thirty360_E : public DayFraction
{
	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::THIRTYE_360; }
	virtual double year_fraction(Date d1, Date d2) const override final
	{
		// int dd1 = d1.day_of_month(), dd2 = d2.day_of_month();
		// int mm1 = d1.month(), mm2 = d2.month();
		// int yy1 = d1.year(), yy2 = d2.year();
		auto ymd1 = date_components(d1);
		int dd1 = ymd1.d, mm1 = ymd1.m, yy1 = ymd1.y;
		auto ymd2 = date_components(d2);
		int dd2 = ymd2.d, mm2 = ymd2.m, yy2 = ymd2.y;

		return (360 * (yy2 - yy1) + 30 * (mm2 - mm1 - 1) + std::max(int(0), 30 - dd1) +
			std::min(int(30), dd2)) /
		       360.0;
	}
	virtual double year_fraction(Date d1, Date d2, bool eom) const override final { return year_fraction(d1, d2); }
	virtual double year_fraction(Date d1, Date d2, Date periodStart, Date periodEnd) const override final
	{
		return year_fraction(d1, d2);
	}
};

/**
 * 30E/360.ISDA
 */
class Thirty360_E_ISDA : public DayFraction
{
	int adjust_day_of_month(YearMonthDay ymd, bool finalPeriod) const
	{
		int dayOfMonth = ymd.d;
		if (ymd.m == Feb && last_day_of_month(ymd.y, ymd.m) == dayOfMonth && !finalPeriod)
			return 30;
		if (dayOfMonth > 30)
			return 30;
		return dayOfMonth;
	}

	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::THIRTYE_360_ISDA; }
	virtual double year_fraction(Date d1, Date d2) const override final
	{
		// TODO log warning
		return year_fraction(d1, d2, false);
	}
	virtual double year_fraction(Date d1, Date d2, bool finalPeriod) const override final
	{
		auto ymd1 = date_components(d1);
		auto ymd2 = date_components(d2);
		int dd1 = adjust_day_of_month(ymd1, false), dd2 = adjust_day_of_month(ymd2, finalPeriod);
		int mm1 = ymd1.m, mm2 = ymd2.m;
		int yy1 = ymd1.y, yy2 = ymd2.y;
		return (360 * (yy2 - yy1) + 30 * (mm2 - mm1 - 1) + std::max(int(0), 30 - dd1) + dd2) / 360.0;
	}
	virtual double year_fraction(Date d1, Date d2, Date periodStart, Date periodEnd) const override final
	{
		// TODO log warning
		return year_fraction(d1, d2, false);
	}
};

/**
 * ACT/ACT
 */
class ActualActual : public DayFraction
{
	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::ACT_ACT_ISDA; }
	virtual double year_fraction(Date d1, Date d2) const override final
	{
		if (d1 == d2)
			return 0.0;
		if (d1 > d2)
			return -year_fraction(d2, d1);
		auto ymd1 = date_components(d1);
		auto ymd2 = date_components(d2);
		int y1 = ymd1.y, y2 = ymd2.y;
		double dib1 = (is_leap(y1) ? 366.0 : 365.0), dib2 = (is_leap(y2) ? 366.0 : 365.0);
		double sum = y2 - y1 - 1;
		sum += (dib1 - day_of_year(ymd1) + 1) / dib1;
		sum += (day_of_year(ymd2) - 1) / dib2;
		return sum;
	}
	virtual double year_fraction(Date d1, Date d2, bool eom) const override final { return year_fraction(d1, d2); }
	virtual double year_fraction(Date d1, Date d2, Date periodStart, Date periodEnd) const override final
	{
		return year_fraction(d1, d2);
	}
};

/**
 * ACT/ACT
 */
class ActualActualISMA : public DayFraction
{
	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::ACT_ACT_ISMA; }
	virtual double year_fraction(Date d1, Date d2) const override final
	{
		// TODO warn
		return year_fraction(d1, d2, d1, d2);
	}
	virtual double year_fraction(Date d1, Date d2, bool eom) const override final
	{
		// TODO warn
		return year_fraction(d1, d2, d1, d2);
	}
	virtual double year_fraction(Date d1, Date d2, Date d3, Date d4) const override final
	{
		if (d1 == d2)
			return 0.0;

		if (d1 > d2)
			return -year_fraction(d2, d1, d3, d4);

		Date refPeriodStart = d3;
		Date refPeriodEnd = d4;

		//        QL_REQUIRE(refPeriodEnd > refPeriodStart &&
		//        refPeriodEnd > d1,
		//                   "invalid reference period: "
		//                   << "date 1: " << d1
		//                   << ", date 2: " << d2
		//                   << ", reference period start: " <<
		//                   refPeriodStart
		//                   << ", reference period end: " <<
		//                   refPeriodEnd);

		// estimate roughly the length in months of a period
		int months = int(0.5 + 12 * double(refPeriodEnd - refPeriodStart) / 365);

		// for short periods...
		if (months == 0) {
			// ...take the reference period as 1 year from d1
			refPeriodStart = d1;
			refPeriodEnd = add(d1, Period(1, PeriodUnit::YEARS));
			months = 12;
		}

		double period = double(months) / 12.0;

		if (d2 <= refPeriodEnd) {
			// here refPeriodEnd is a future (notional?) payment
			// date
			if (d1 >= refPeriodStart) {
				// here refPeriodStart is the last (maybe
				// notional) payment date. refPeriodStart <= d1
				// <= d2 <= refPeriodEnd [maybe the equality
				// should be enforced, since refPeriodStart < d1
				// <= d2 < refPeriodEnd could give wrong
				// results] ???
				return period * double((d2 - d1)) / double(refPeriodEnd - refPeriodStart);
			} else {
				// here refPeriodStart is the next (maybe
				// notional) payment date and refPeriodEnd is
				// the second next (maybe notional) payment
				// date. d1 < refPeriodStart < refPeriodEnd AND
				// d2 <= refPeriodEnd this case is long first
				// coupon

				// the last notional payment date
				Date previousRef = sub(refPeriodStart, Period(months, PeriodUnit::MONTHS));
				if (d2 > refPeriodStart)
					return year_fraction(d1, refPeriodStart, previousRef, refPeriodStart) +
					       year_fraction(refPeriodStart, d2, refPeriodStart, refPeriodEnd);
				else
					return year_fraction(d1, d2, previousRef, refPeriodStart);
			}
		} else {
			// here refPeriodEnd is the last (notional?) payment
			// date d1 < refPeriodEnd < d2 AND refPeriodStart <
			// refPeriodEnd
			//            QL_REQUIRE(refPeriodStart<=d1,
			//                       "invalid dates: "
			//                       "d1 < refPeriodStart <
			//                       refPeriodEnd < d2");
			// now it is: refPeriodStart <= d1 < refPeriodEnd < d2

			// the part from d1 to refPeriodEnd
			double sum = year_fraction(d1, refPeriodEnd, refPeriodStart, refPeriodEnd);

			// the part from refPeriodEnd to d2
			// count how many regular periods are in [refPeriodEnd,
			// d2], then add the remaining time
			int i = 0;
			Date newRefStart = refPeriodStart, newRefEnd = refPeriodEnd;
			for (;;) {
				newRefStart = add(refPeriodEnd, Period(months * i, PeriodUnit::MONTHS));
				newRefEnd = add(refPeriodEnd, Period(months * (i + 1), PeriodUnit::MONTHS));
				if (d2 < newRefEnd) {
					break;
				} else {
					sum += period;
					i++;
				}
			}
			sum += year_fraction(newRefStart, d2, newRefStart, newRefEnd);
			return sum;
		}
	}
};

static const int minYear_ = 1900;
static const int maxYear_ = 2199; // We have a bug that the cache loop attempts to access 1 Jan 2200
static const int numYears = (maxYear_ - minYear_ + 1);
static const int numMonths = numYears * 12 + 1;

class Business252 : public DayFraction
{
	const Calendar *calendar_;

	// cache data for every month
	int *monthlyFigures_;
	// cache data for every year
	int *yearlyFigures_;

	int offset(int y, int m) const
	{
		y -= minYear_;
		int n = y * 12 + m;
		return n;
	}

	bool same_year(Date d1, Date d2) const { return date_components(d1).y == date_components(d2).y; }

	bool same_month(Date d1, Date d2) const
	{
		return same_year(d1, d2) && date_components(d1).m == date_components(d2).m;
	}

	int business_days(int month, int year) const { return monthlyFigures_[offset(year, month)]; }

	int business_days(int year) const { return yearlyFigures_[year - minYear_]; }

	int count_days(Date d1, Date d2) const
	{
		Period oneMonth(1, PeriodUnit::MONTHS);
		Period oneYear(1, PeriodUnit::YEARS);

		if (same_month(d1, d2) || d1 >= d2) {
			// we treat the case of d1 > d2 here, since we'd need a
			// second cache to get it right (our cached figures are
			// for first included, last excluded and might have to
			// be changed going the other way.)
			return calendar_->business_days_between(d1, d2);
		} else if (same_year(d1, d2)) {
			int total = 0;
			Date d = make_date(1, date_components(d1).m, date_components(d1).y);
			// first, we get to the beginning of next month.
			d = add(d, oneMonth);
			total += calendar_->business_days_between(d1, d);
			// then, we add any whole months (whose figures might be
			// cached already) in the middle of our period.
			while (!same_month(d, d2)) {
				total += business_days(date_components(d).m, date_components(d).y);
				d = add(d, oneMonth);
			}
			// finally, we get to the end of the period.
			total += calendar_->business_days_between(d, d2);
			return total;
		} else {
			int total = 0;
			Date d = make_date(1, date_components(d1).m, date_components(d1).y);
			// first, we get to the beginning of next year.
			// The first bit gets us to the end of this month...
			d = add(d, oneMonth);
			total += calendar_->business_days_between(d1, d);
			// ...then we add any remaining months, possibly cached
			for (int m = date_components(d1).m + 1; m <= 12; ++m) {
				total += business_days(m, date_components(d).y);
			}
			// then, we add any whole year in the middle of our
			// period.
			d = make_date(1, January, date_components(d1).y + 1);
			while (!same_year(d, d2)) {
				total += business_days(date_components(d).y);
				d = add(d, oneYear);
			}
			// finally, we get to the end of the period.
			// First, we add whole months...
			for (int m = 1; m < date_components(d2).m; ++m) {
				total += business_days(m, date_components(d2).y);
			}
			// ...then the last bit.
			d = make_date(1, date_components(d2).m, date_components(d2).y);
			total += calendar_->business_days_between(d, d2);
			return total;
		}
	}

	public:
	virtual DayCountFraction id() const override final { return DayCountFraction::BUS_252; }
	Business252(const Calendar *c) : calendar_(c)
	{
		monthlyFigures_ = new int[numMonths];
		yearlyFigures_ = new int[numYears];
		for (int y = minYear_; y <= maxYear_; y++) {
			// calculate and store.
			int total = 0;
			for (int m = 1; m <= 12; m++) {
				// calculate and store.
				Date d1 = make_date(1, m, y);
				Date d2 = add(d1, Period(1, PeriodUnit::MONTHS));
				int mdays = c->business_days_between(d1, d2);
				int pos = offset(y, m);
				assert(pos < numMonths);
				monthlyFigures_[pos] = mdays;
				total += mdays;
			}
			int pos = y - minYear_;
			assert(pos < numYears);
			yearlyFigures_[pos] = total;
		}
	}
	virtual ~Business252()
	{
		delete[] yearlyFigures_;
		delete[] monthlyFigures_;
	}
	virtual double year_fraction(Date d1, Date d2) const override final { return count_days(d1, d2) / 252.0; }
	virtual double year_fraction(Date d1, Date d2, bool eom) const override final { return year_fraction(d1, d2); }
	virtual double year_fraction(Date d1, Date d2, Date periodStart, Date periodEnd) const override final
	{
		return year_fraction(d1, d2);
	}
};

static const Actual360 Act360;
static const Actual365Fixed Act365F;
static const Thirty360_ Thirty360;
static const Thirty360_E ThirtyE360;
static const Thirty360_E_ISDA ThirtyE360ISDA;
static const ActualActual ActActISDA;
static const ActualActualISMA ActActISMA;

static const DayFraction *fractions[]{
    nullptr,	 // DAY_COUNT_FRACTION_UNSPECIFIED = 0;
    &Act360,	 // ACT_360 = 1;
    &Act365F,	// ACT_365_FIXED = 2;
    &Thirty360,      // THIRTY_360 = 3;
    &ThirtyE360,     // THIRTYE_360 = 4;
    &ThirtyE360ISDA, // THIRTYE_360_ISDA = 5;
    nullptr,	 // ONE_ON_ONE = 6;
    &ActActISDA,     // ACT_ACT_ISDA = 7;
    &ActActISMA,     // ACT_ACT_ISMA = 8;
    nullptr,	 // BUS_252 = 9;
};

const DayFraction *get_day_fraction(DayCountFraction dfc)
{
	assert(dfc > 0 && dfc < DayCountFraction::BUS_252);
	return fractions[dfc];
}

// cache by calendar name
static std::mutex cache_bus_252_mutex;
static std::map<BusinessCenter, std::unique_ptr<Business252>> cached_bus_252;

const DayFraction *get_bus_252(CalendarService *calendarService, BusinessCenter center)
{
	// BUS/252 requires special handling
	const Calendar *calendar = calendarService->get_calendar(center);
	if (calendar == nullptr)
		return nullptr;
	std::lock_guard<std::mutex> lock(cache_bus_252_mutex);
	if (cached_bus_252.find(center) == std::end(cached_bus_252)) {
		cached_bus_252[center] = std::make_unique<Business252>(calendar);
	}
	return cached_bus_252[center].get();
}

static int SingleCase(const DayFraction &df, Date d1, Date d2, double expected)
{
	double value = df.year_fraction(d1, d2);
	if (std::fabs(value - expected) > 0.000001) {
		//    std::cerr << df.name() << " failed for " <<
		//      d1 << ", " << d2 << " expected " << expected
		//      << " got " << value;
		return 1;
	}
	return 0;
}

static int SingleCase(const DayFraction &df, Date d1, Date d2, Date d3, Date d4, double expected)
{
	double value = df.year_fraction(d1, d2, d3, d4);
	if (std::fabs(value - expected) > 0.00000001) {
		//    std::cerr << df.name() << " failed for " <<
		//      d1 << ", " << d2 << " expected " << expected
		//      << " got " << value;
		return 1;
	}
	return 0;
}

static int test_dayfraction()
{
	int failure_count = 0;
	failure_count +=
	    SingleCase(ActActISDA, make_date(1, Month::November, 2003), make_date(1, Month::May, 2004), 0.497724380567);
	failure_count +=
	    SingleCase(ActActISMA, make_date(1, Month::November, 2003), make_date(1, Month::May, 2004),
		       make_date(1, Month::November, 2003), make_date(1, Month::May, 2004), 0.500000000000);
	// short first calculation period (first period)
	failure_count += SingleCase(ActActISDA, make_date(1, Month::February, 1999), make_date(1, Month::July, 1999),
				    0.410958904110);
	failure_count += SingleCase(ActActISMA, make_date(1, Month::February, 1999), make_date(1, Month::July, 1999),
				    make_date(1, Month::July, 1998), make_date(1, Month::July, 1999), 0.410958904110);
	// short first calculation period (second period)
	failure_count +=
	    SingleCase(ActActISDA, make_date(1, Month::July, 1999), make_date(1, Month::July, 2000), 1.001377348600);
	failure_count += SingleCase(ActActISMA, make_date(1, Month::July, 1999), make_date(1, Month::July, 2000),
				    make_date(1, Month::July, 1999), make_date(1, Month::July, 2000), 1.000000000000);
	// long first calculation period (first period)
	failure_count += SingleCase(ActActISDA, make_date(15, Month::August, 2002), make_date(15, Month::July, 2003),
				    0.915068493151);
	failure_count +=
	    SingleCase(ActActISMA, make_date(15, Month::August, 2002), make_date(15, Month::July, 2003),
		       make_date(15, Month::January, 2003), make_date(15, Month::July, 2003), 0.915760869565);
	// long first calculation period (second period)
	/* Warning: the ISDA case is in disagreement with mktc1198.pdf */
	failure_count += SingleCase(ActActISDA, make_date(15, Month::July, 2003), make_date(15, Month::January, 2004),
				    0.504004790778);
	failure_count +=
	    SingleCase(ActActISMA, make_date(15, Month::July, 2003), make_date(15, Month::January, 2004),
		       make_date(15, Month::July, 2003), make_date(15, Month::January, 2004), 0.500000000000);
	// short final calculation period (penultimate period)
	failure_count += SingleCase(ActActISDA, make_date(30, Month::July, 1999), make_date(30, Month::January, 2000),
				    0.503892506924);
	failure_count +=
	    SingleCase(ActActISMA, make_date(30, Month::July, 1999), make_date(30, Month::January, 2000),
		       make_date(30, Month::July, 1999), make_date(30, Month::January, 2000), 0.500000000000);
	// short final calculation period (final period)
	failure_count += SingleCase(ActActISDA, make_date(30, Month::January, 2000), make_date(30, Month::June, 2000),
				    0.415300546448);
	failure_count +=
	    SingleCase(ActActISMA, make_date(30, Month::January, 2000), make_date(30, Month::June, 2000),
		       make_date(30, Month::January, 2000), make_date(30, Month::July, 2000), 0.417582417582);
	return failure_count;
}

static int test_Business252()
{
	std::vector<Date> testDates;
	testDates.push_back(make_date(1, February, 2002));
	testDates.push_back(make_date(4, February, 2002));
	testDates.push_back(make_date(16, May, 2003));
	testDates.push_back(make_date(17, December, 2003));
	testDates.push_back(make_date(17, December, 2004));
	testDates.push_back(make_date(19, December, 2005));
	testDates.push_back(make_date(2, January, 2006));
	testDates.push_back(make_date(13, March, 2006));
	testDates.push_back(make_date(15, May, 2006));
	testDates.push_back(make_date(17, March, 2006));
	testDates.push_back(make_date(15, May, 2006));
	testDates.push_back(make_date(26, July, 2006));
	testDates.push_back(make_date(28, June, 2007));
	testDates.push_back(make_date(16, September, 2009));
	testDates.push_back(make_date(26, July, 2016));

	double expected[] = {0.0039682539683, 1.2738095238095, 0.6031746031746, 0.9960317460317,  1.0000000000000,
			     0.0396825396825, 0.1904761904762, 0.1666666666667, -0.1507936507937, 0.1507936507937,
			     0.2023809523810, 0.912698412698,  2.214285714286,  6.84126984127};

	const DayFraction *dayCounter1 = get_bus_252(get_calendar_factory(), BusinessCenter::BRSP);

	double calculated;
	int failure_count = 0;
	for (int i = 1; i < testDates.size(); i++) {
		calculated = dayCounter1->year_fraction(testDates[i - 1], testDates[i]);
		if (std::fabs(calculated - expected[i - 1]) > 1.0e-12) {
			//      BOOST_ERROR("from " << testDates[i - 1]
			//        << " to " << testDates[i] << ":\n"
			//        << std::setprecision(14)
			//        << "    calculated: " << calculated << "\n"
			//        << "    expected:   " << expected[i - 1]);
			failure_count++;
		}
	}

	return failure_count;
}

static int test_search_dayfraction()
{
	static DayCountFraction fractions[] = {DayCountFraction::ACT_360,	  DayCountFraction::ACT_365_FIXED,
					       DayCountFraction::THIRTY_360,       DayCountFraction::THIRTYE_360,
					       DayCountFraction::THIRTYE_360_ISDA, DayCountFraction::ACT_ACT_ISDA,
					       DayCountFraction::ACT_ACT_ISMA};
	int failure_count = 0;
	for (int i = 0; i < 7; i++) {
		DayCountFraction cp = get_day_fraction(fractions[i])->id();
		int failed = cp != fractions[i];
		// if (failed) std::cerr << "Failed search for " << names[i] <<
		// "\n";
		failure_count += failed;
	}
	return failure_count;
}

int test_dayfractions()
{
	int failure_count = 0;
	failure_count += test_search_dayfraction();
	failure_count += test_dayfraction();
	failure_count += test_Business252();
	if (failure_count == 0)
		printf("DayFraction Tests OK\n");
	else
		printf("DayFraction Tests FAILED\n");
	return failure_count;
}

} // namespace redukti
