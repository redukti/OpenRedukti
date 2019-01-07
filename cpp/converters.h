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

#ifndef _REDUKTI_CONVERTERS_H
#define _REDUKTI_CONVERTERS_H

#include <enums.pb.h>

namespace redukti
{

class Period;

class Converter
{
      public:
	virtual ~Converter() {}
	virtual BusinessCenter business_center_from_string(const char *value) const = 0;
	virtual const char *business_center_to_string(BusinessCenter value) const = 0;
	virtual BusinessDayConvention business_day_convention_from_string(const char *s) const = 0;
	virtual PeriodUnit period_unit_from_string(const char *s) const = 0;
	virtual bool period_from_string(const char *periodName, Period *p) const = 0;
	virtual DayCountFraction day_count_fraction_from_string(const char *value) const = 0;
	virtual const char *day_count_fraction_to_string(DayCountFraction value) const = 0;
	virtual Tenor tenor_from_period_unit_and_len(PeriodUnit unit, int value) const = 0;
	virtual bool tenor_to_period_unit_and_multiplier(Tenor value, PeriodUnit *unit, int *mult) const = 0;
	virtual std::string tenor_to_string(Tenor tenor) const = 0;
	virtual RollConvention roll_convention_from_string(const char *s) const = 0;
	virtual Currency currency_from_string(const char *s) const = 0;
	virtual const char *currency_to_string(Currency value) const = 0;
	virtual IsdaIndex isda_index_from_string(const char *s) const = 0;
	virtual const char *isda_index_to_string(IsdaIndex value) const = 0;
	virtual CompoundingMethod compounding_method_from_string(const char *value) const = 0;
	virtual IndexFamily index_family_from_string(const char *value) const = 0;
	virtual const char *index_family_to_string(IndexFamily value) const = 0;
	virtual const char *period_unit_to_string(PeriodUnit period_unit) const = 0;
	virtual int tenor_to_days(Tenor tenor) const = 0;
	virtual InterpolatorType interpolator_type_from_string(const char *s) const = 0;
	virtual PricingCurveType pricing_curve_type_from_string(const char *s) const = 0;
	virtual IRRateType rate_type_from_string(const char *s) const = 0;
	virtual CurveGroup curve_group_from_string(const char *value) const = 0;
	virtual MaturityGenerationRule maturity_generation_rule_from_string(const char *value) const = 0;
	virtual JointCalendarRule joint_calendar_rule_from_string(const char *value) const = 0;
};

extern const Converter *get_default_converter();
extern int test_conversions();

} // namespace redukti

#endif