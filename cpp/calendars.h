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
#ifndef _REDUKTI_CALENDARS_H
#define _REDUKTI_CALENDARS_H

#include <date.h>
#include <enums.pb.h>

#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace redukti
{

// The Calendar interface provides the means to determine whether
// a given date is a holiday for a business center or not. Also
// the interface provides methods for adjusting dates as per the
// holiday calendar.
// Immutable for thread safety.
class Calendar
{
	public:
	virtual ~Calendar() noexcept {}
	virtual int id() const noexcept = 0;
	// Returns all the ids - relevant for calendars made by combining
	// others
	virtual void get_ids(std::array<BusinessCenter, 4> &ids) const noexcept
	{
		std::fill(std::begin(ids), std::end(ids), BusinessCenter::BUSINESS_CENTER_UNSPECIFIED);
		ids[0] = (BusinessCenter)id();
	}
	virtual bool is_holiday(Date d) const noexcept = 0;
	bool is_businessday(Date d) const noexcept;
	bool is_end_of_month(Date d) const noexcept;

	// Adjust the given date to be the last business day of the month
	Date end_of_month(Date d) const noexcept;
	// Adjusts a non-business day to the appropriate near business day
	//  with respect to the given convention.
	Date adjust(Date date, BusinessDayConvention convention = BusinessDayConvention::FOLLOWING) const noexcept;
	//  Advances the given date of the given number of business days and
	//  returns the result. Note that if unit is Days then business day
	// convention and eom flags are not used as the date is move by the
	// specified business days. For other period units the date is moved as
	// per raw calendar and then adjusted if it falls on a holiday
	Date advance(Date date, int n, PeriodUnit unit,
		     BusinessDayConvention convention = BusinessDayConvention::FOLLOWING, bool end_of_month = false) const
	    noexcept;
	//  Advances the given date as specified by the given period and
	//  returns the result.
	//  The input date is not modified.
	Date advance(Date date, const Period &period,
		     BusinessDayConvention convention = BusinessDayConvention::FOLLOWING, bool end_of_month = false) const
	    noexcept;

	// Calculates the number of business days between two given
	// dates and returns the result.
	//
	int business_days_between(Date from, Date to, bool include_first = true, bool include_last = false) const
	    noexcept;
};

struct JointCalendarParameters {
	std::array<BusinessCenter, 4> centers;
	JointCalendarParameters(BusinessCenter center1 = BusinessCenter::BUSINESS_CENTER_UNSPECIFIED, // For Cython
				BusinessCenter center2 = BusinessCenter::BUSINESS_CENTER_UNSPECIFIED, // For Cython
				BusinessCenter center3 = BusinessCenter::BUSINESS_CENTER_UNSPECIFIED,
				BusinessCenter center4 = BusinessCenter::BUSINESS_CENTER_UNSPECIFIED)
	{
		assert(BusinessCenter_IsValid(center1));
		assert(BusinessCenter_IsValid(center2));
		assert(BusinessCenter_IsValid(center3));
		assert(BusinessCenter_IsValid(center4));
		centers[0] = center1;
		centers[1] = center2;
		centers[2] = center3;
		centers[3] = center4;
	}
};

// The Calendar Service manages calendar instances. It has to meet following
// requirements: a) It must always return the same Calendar instance for a given
// business center. Clients
//    can assume that the instance will not go away or change in any way as long
//    as the service is live.
// b) Ditto for joint calendar instances.
// c) Calendar instances must be immutable.
class CalendarService
{
	public:
	virtual ~CalendarService() {}
	// Return the calendar specified. Memory is managed by the
	// CalendarFactory so the caller must not delete.
	virtual const Calendar *get_calendar(BusinessCenter id) noexcept = 0;

	// Set a calendar to given instance.
	// The service will take ownership of the instance
	// May fail if calendar instance already set and has been
	// accessed by a client - i.e. new calendars can only be set prior to
	// any use.
	virtual bool set_calendar(BusinessCenter id, std::unique_ptr<Calendar> calendar) noexcept = 0;

	// Create a calendar from a set of holidays and assign it to the business center
	// If the assignment is successful the service will take ownership of the instance
	// May fail if calendar instance already set and has been
	// accessed by a client - i.e. new calendars can only be set prior to
	// any use.
	virtual bool set_calendar(BusinessCenter id, const Date *holidays, size_t n) noexcept = 0;

	// Create joint calendar
	// Note that the order in which the business centers are given
	// should not matter - i.e. the constituents must be sorted and then
	// combined so that for a given combination the returned instance is
	// always the same
	virtual Calendar *get_calendar(JointCalendarParameters calendars,
				       JointCalendarRule rule = JointCalendarRule::JOIN_HOLIDAYS) noexcept = 0;
};

// Get the calendar factory
extern CalendarService *get_calendar_factory() noexcept;

// Utility for constructing a joint calendar
extern const Calendar *build_calendar(CalendarService *calendar_service,
				      const google::protobuf::RepeatedField<google::protobuf::int32> &values,
				      JointCalendarRule rule = JointCalendarRule::JOIN_HOLIDAYS);

// Utility for constructing a joint calendar
const Calendar *build_calendar(CalendarService *calendar_service, const std::vector<BusinessCenter> &values,
			       JointCalendarRule rule = JointCalendarRule::JOIN_HOLIDAYS);

// Adjust the given date to be the last business day of the month
inline Date Calendar::end_of_month(Date d) const noexcept
{
	YearMonthDay ymd = date_components(d);
	return adjust(redukti::end_of_month(ymd.y, ymd.m), BusinessDayConvention::PRECEDING);
}

inline bool Calendar::is_businessday(Date d) const noexcept { return !is_holiday(d); }

extern int test_calendars();

} // namespace redukti

#endif
