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

#ifndef _REDUKTI_SCHEDULE_H
#define _REDUKTI_SCHEDULE_H

#include <date.h>
#include <enums.pb.h>
#include <status.h>

namespace redukti
{

class ScheduleParameters;
class Schedule;

// Build a schedule as per the schedule parameters
// If succesful returns true
//
// The minimum set of required parameters are:
// 1. effective_date
// 2. termination_date or term
// 3. payment_frequency
//
// Uses calendar service to obtain relevant calendars
extern StatusCode build_schedule(ScheduleParameters &params, Schedule &schedule) noexcept;

// Adjusts a date as per roll convention specified
extern Date adjust_date(Date d, RollConvention rc) noexcept;

extern int test_schedule_generation();
} // namespace redukti

#endif