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
#ifndef _REDUKTI_DAYFRACTIONS_H
#define _REDUKTI_DAYFRACTIONS_H

#include <calendars.h>
#include <date.h>

namespace redukti
{

// Compute the difference between dates as per Day Count Convention.
// The difference is measured in factional units of a year, where one year 1.0.
// Must be immutable and thread-safe.
// Clients must be able to hold references to these for the lifetime of
// the application.
class DayFraction
{
	public:
	virtual ~DayFraction() {}

	// Calculate the difference d2-d2 as per convention
	// for the DayFraction; value is a decimal expressed as a year fraction.
	// So 1.0 means 1 year.
	virtual double year_fraction(Date d1, Date d2) const = 0;

	// Only used for ThirtyE360ISDA (30E/360.ISDA)
	// The finalPeriod flag indicates whether this fraction is for the
	// final period - i.e. d2 is maturity date. So typically,
	// when calculating the last calc period in a swap, this flag must be
	// set to true.
	virtual double year_fraction(Date d1, Date d2, bool finalPeriod) const = 0;

	// Used only for ACT/ACT.ISMA
	// refStart - If regular period or front stub then adjusted end date
	// 	minus calculation period frequency (roll convention NONE),
	// 	else adjusted start date
	// refEnd - If regular period or front stub then adjusted end date,
	// 	else adjusted start date plus calculation period
	// 	frequency (roll convention NONE)
	virtual double year_fraction(Date d1, Date d2, Date refStart, Date refEnd) const = 0;

	// Returns the ISDA name
	virtual DayCountFraction id() const = 0;
};

// Get an instance of a DayFraction
// Requirements:
// a) There must only be one instance associated with a particular DayCountFraction
// b) The DayCountFraction implementation must be immutable and hence thread-safe
extern const DayFraction *get_day_fraction(DayCountFraction dfc);

// The BUS252 day fraction requires a calendar.
// Requirements:
// a) There must only be one instance associated with a particular DayCountFraction
// b) The DayCountFraction implementation must be immutable and hence thread-safe
extern const DayFraction *get_bus_252(CalendarService *calendarService, BusinessCenter center);

extern int test_dayfractions();

} // namespace redukti

#endif
