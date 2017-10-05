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

#include <converters.h>
#include <date.h>
#include <logger.h>

#include <algorithm>
#include <cstring>
#include <utility>

#ifdef _MSC_VER
#define strcasecmp _stricmp
#endif

namespace redukti
{

static int map_tenor_to_days[] = {
    0,	// UNSPECIFIED
    1,	// 1D
    2,	// 2D
    3,	// 3D
    4,	// 4D
    5,	// 5D
    6,	// 6D
    7,	// 1W
    14,       // 2W
    21,       // 3W
    30,       // 1M
    60,       // 2M
    91,       // 3M
    121,      // 4M
    152,      // 5M
    182,      // 6M
    212,      // 7M
    243,      // 8M
    273,      // 9M
    304,      // 10M
    334,      // 11M
    365,      // 12M
    395,      // 13M
    425,      // 14M
    456,      // 15M
    486,      // 16M
    517,      // 17M
    547,      // 18M
    577,      // 19M
    608,      // 20M
    638,      // 21M
    669,      // 22M
    699,      // 23M
    365 * 2,  // 2Y
    365 * 3,  // 3Y
    365 * 4,  // 4Y
    365 * 5,  // 5Y
    365 * 6,  // 6Y
    365 * 7,  // 7Y
    365 * 8,  // 8Y
    365 * 9,  // 9Y
    365 * 10, // 10Y
    365 * 11, // 11Y
    365 * 12, // 12Y
    365 * 13, // 13Y
    365 * 14, // 14Y
    365 * 15, // 15Y
    365 * 16, // 16Y
    365 * 17, // 17Y
    365 * 18, // 18Y
    365 * 19, // 19Y
    365 * 20, // 20Y
    365 * 21, // 21Y
    365 * 22, // 22Y
    365 * 23, // 23Y
    365 * 24, // 24Y
    365 * 25, // 25Y
    365 * 26, // 26Y
    365 * 27, // 27Y
    365 * 28, // 28Y
    365 * 29, // 29Y
    365 * 30, // 30Y
    365 * 31, // 31Y
    365 * 32, // 32Y
    365 * 33, // 33Y
    365 * 34, // 34Y
    365 * 35, // 35Y
    365 * 36, // 36Y
    365 * 37, // 37Y
    365 * 38, // 38Y
    365 * 39, // 39Y
    365 * 40, // 40Y
    365 * 41, // 41Y
    365 * 42, // 42Y
    365 * 43, // 43Y
    365 * 44, // 44Y
    365 * 45, // 45Y
    365 * 46, // 46Y
    365 * 47, // 47Y
    365 * 48, // 48Y
    365 * 49, // 49Y
    365 * 50, // 50Y
    365 * 51, // 51Y
    365 * 52, // 52Y
    365 * 53, // 53Y
    365 * 54, // 54Y
    365 * 55, // 55Y
    365 * 56, // 56Y
    365 * 57, // 57Y
    365 * 58, // 58Y
    365 * 59, // 59Y
    365 * 60, // 60Y
    365 * 61, // 61Y
    365 * 62, // 62Y
    365 * 63, // 63Y
    365 * 64, // 64Y
    365 * 65, // 65Y
};

static std::pair<const char *, Currency> currencies[] = {
    {"USD", Currency::USD}, {"GBP", Currency::GBP}, {"EUR", Currency::EUR}, {"JPY", Currency::JPY},
    {"AUD", Currency::AUD}, {"CAD", Currency::CAD}, {"CHF", Currency::CHF}, {"CZK", Currency::CZK},
    {"DKK", Currency::DKK}, {"HKD", Currency::HKD}, {"NOK", Currency::NOK}, {"NZD", Currency::NZD},
    {"SEK", Currency::SEK}, {"SGD", Currency::SGD}, {"MXN", Currency::MXN}, {"ZAR", Currency::ZAR},
    {"HUF", Currency::HUF}, {"PLN", Currency::PLN}};

static std::pair<const char *, BusinessCenter> business_centers[] = {
    {"AUSY", BusinessCenter::AUSY}, {"BRSP", BusinessCenter::BRSP}, {"EUTA", BusinessCenter::EUTA},
    {"GBLO", BusinessCenter::GBLO}, {"JPTO", BusinessCenter::JPTO}, {"USNY", BusinessCenter::USNY},
};

static std::pair<const char *, IndexFamily> indices[] = {
    {"LIBOR", LIBOR}, {"FEDFUND", FEDFUND}, {"EURIBOR", EURIBOR}, {"SONIA", SONIA}, {"AONIA", AONIA},
    {"BBSW", BBSW},   {"CDOR", CDOR},       {"CORRA", CORRA},     {"TOIS", TOIS},   {"PRIBOR", PRIBOR},
    {"CIBOR", CIBOR}, {"HIBOR", HIBOR},     {"BUBOR", BUBOR},     {"TONA", TONA},   {"NIBOR", NIBOR},
    {"BKBM", BKBM},   {"WIBOR", WIBOR},     {"STIBOR", STIBOR},   {"SOR", SOR},     {"JIBAR", JIBAR},
    {"EONIA", EONIA}};

static std::pair<const char *, IsdaIndex> isdaindices[] = {
    {"AUD-BBR-BBSW", IsdaIndex::AUD_BBR_BBSW},
    {"AUD-LIBOR-BBA", IsdaIndex::AUD_LIBOR_BBA},
    {"AUD-AONIA-OIS-COMPOUND", IsdaIndex::AUD_AONIA_OIS_COMPOUND},
    {"CAD-BA-CDOR", IsdaIndex::CAD_BA_CDOR},
    {"CAD-LIBOR-BBA", IsdaIndex::CAD_LIBOR_BBA},
    {"CAD-CORRA-OIS-COMPOUND", IsdaIndex::CAD_CORRA_OIS_COMPOUND},
    {"CHF-LIBOR-BBA", IsdaIndex::CHF_LIBOR_BBA},
    {"CHF-TOIS-OIS-COMPOUND", IsdaIndex::CHF_TOIS_OIS_COMPOUND},
    {"CZK-PRIBOR-PRBO", IsdaIndex::CZK_PRIBOR_PRBO},
    {"DKK-CIBOR2-DKNA13", IsdaIndex::DKK_CIBOR2_DKNA13},
    {"DKK-CIBOR-DKNA13", IsdaIndex::DKK_CIBOR_DKNA13},
    {"EUR-EURIBOR-Reuters", IsdaIndex::EUR_EURIBOR_Reuters},
    {"EUR-EURIBOR-Telerate", IsdaIndex::EUR_EURIBOR_Telerate},
    {"EUR-LIBOR-BBA", IsdaIndex::EUR_LIBOR_BBA},
    {"EUR-EONIA-OIS-COMPOUND", IsdaIndex::EUR_EONIA_OIS_COMPOUND},
    {"GBP-LIBOR-BBA", IsdaIndex::GBP_LIBOR_BBA},
    {"GBP-WMBA-SONIA-COMPOUND", IsdaIndex::GBP_WMBA_SONIA_COMPOUND},
    {"HKD-HIBOR-HIBOR=", IsdaIndex::HKD_HIBOR_HIBOR},
    {"HKD-HIBOR-HKAB", IsdaIndex::HKD_HIBOR_HKAB},
    {"HKD-HIBOR-ISDC", IsdaIndex::HKD_HIBOR_ISDC},
    {"HUF-BUBOR-Reuters", IsdaIndex::HUF_BUBOR_Reuters},
    {"JPY-LIBOR-BBA", IsdaIndex::JPY_LIBOR_BBA},
    {"JPY-TONA-OIS-COMPOUND", IsdaIndex::JPY_TONA_OIS_COMPOUND},
    {"MXN-USDMXNBASIS", IsdaIndex::MXN_USDMXNBASIS},
    {"MXN-TIIE-Banxico", IsdaIndex::MXN_TIIE_Banxico},
    {"NOK-NIBOR-NIBR", IsdaIndex::NOK_NIBOR_NIBR},
    {"NZD-BBR-FRA", IsdaIndex::NZD_BBR_FRA},
    {"NZD-BBR-Telerate", IsdaIndex::NZD_BBR_Telerate},
    {"PLN-WIBOR-WIBO", IsdaIndex::PLN_WIBOR_WIBO},
    {"SEK-STIBOR-SIDE", IsdaIndex::SEK_STIBOR_SIDE},
    {"SGD-SOR-Reuters", IsdaIndex::SGD_SOR_Reuters},
    {"SGD-SOR-VWAP", IsdaIndex::SGD_SOR_VWAP},
    {"USD-LIBOR-BBA", IsdaIndex::USD_LIBOR_BBA},
    {"USD-Federal Funds-H.15-OIS-COMPOUND", IsdaIndex::USD_Federal_Funds_H_15_OIS_COMPOUND},
    {"USD-Federal Funds-H.15", IsdaIndex::USD_Federal_Funds_H_15},
    {"ZAR-JIBAR-SAFEX", IsdaIndex::ZAR_JIBAR_SAFEX},
    {"UK-RPI", IsdaIndex::UK_RPI},
    {"EUR-EXT-CPI", IsdaIndex::EUR_EXT_CPI},
    {"FRC-EXT-CPI", IsdaIndex::FRC_EXT_CPI},
    {"USA-CPI-U", IsdaIndex::USA_CPI_U}};

static std::pair<const char *, InterpolatorType> interpolators[] = {
    {"CubicSplineClamped", InterpolatorType::CUBIC_SPLINE_CLAMPED},
    {"CubicSplineNatural", InterpolatorType::CUBIC_SPLINE_NATURAL},
    {"CubicSplineNotAKnot", InterpolatorType::CUBIC_SPLINE_NOT_A_KNOT},
    {"FlatLeft", InterpolatorType::FLAT_LEFT},
    {"FlatRight", InterpolatorType::FLAT_RIGHT},
    {"Linear", InterpolatorType::LINEAR},
    {"LogCubicSplineNatural", InterpolatorType::LOG_CUBIC_SPLINE_NATURAL},
    {"LogLinear", InterpolatorType::LOG_LINEAR},
    {"MonotoneConvex", InterpolatorType::MONOTONE_CONVEX},
};

static std::pair<const char *, CurveGroup> curve_groups[] = {
    {"A", CurveGroup::CURVE_GROUP_A}, {"B", CurveGroup::CURVE_GROUP_B}, {"C", CurveGroup::CURVE_GROUP_C},
    {"D", CurveGroup::CURVE_GROUP_D}, {"E", CurveGroup::CURVE_GROUP_E}, {"F", CurveGroup::CURVE_GROUP_F},
};

static std::pair<const char *, MaturityGenerationRule> maturity_genrules[] = {
    {"DeriveFromInstruments", MaturityGenerationRule::MATURITY_GENERATION_RULE_DERIVE_FROM_INSTRUMENTS},
    {"FixedTenors", MaturityGenerationRule::MATURITY_GENERATION_RULE_FIXED_TENORS}};

static std::pair<const char *, DayCountFraction> day_count_fractions[] = {
    {"1/1", DayCountFraction::ONE_ON_ONE},
    {"30/360", DayCountFraction::THIRTY_360},
    {"30E/360", DayCountFraction::THIRTYE_360},
    {"30E/360.ISDA", DayCountFraction::THIRTYE_360_ISDA},
    {"ACT/360", DayCountFraction::ACT_360},
    {"ACT/365.FIXED", DayCountFraction::ACT_365_FIXED},
    {"ACT/ACT.ISDA", DayCountFraction::ACT_ACT_ISDA},
    {"ACT/ACT.ISMA", DayCountFraction::ACT_ACT_ISMA},
    {"BUS/252", DayCountFraction::BUS_252}};

class ValueConverter : public Converter
{
      public:
	BusinessCenter business_center_from_string(const char *value) const override final
	{
		auto iter = std::find_if(
		    std::begin(business_centers), std::end(business_centers),
		    [value](const std::pair<const char *, BusinessCenter> &v) { return strcmp(value, v.first) == 0; });
		if (iter == std::end(business_centers))
			return BusinessCenter::BUSINESS_CENTER_UNSPECIFIED;
		else
			return iter->second;
	}

	const char *business_center_to_string(BusinessCenter value) const override final
	{
		auto iter = std::find_if(
		    std::begin(business_centers), std::end(business_centers),
		    [value](const std::pair<const char *, BusinessCenter> &v) { return value == v.second; });
		if (iter == std::end(business_centers))
			return "";
		else
			return iter->first;
	}

	BusinessDayConvention business_day_convention_from_string(const char *s) const override final
	{
		if (s == nullptr || strlen(s) < 4)
			return BusinessDayConvention::UNADJUSTED;
		switch (*s) {
		case 'M': {
			switch (s[3]) {
			case 'F':
				return BusinessDayConvention::MODIFIED_FOLLOWING;
			case 'P':
				return BusinessDayConvention::MODIFIED_PRECEDING;
			}
		case 'F':
			return BusinessDayConvention::FOLLOWING;
		case 'P':
			return BusinessDayConvention::PRECEDING;
		case 'N':
			return BusinessDayConvention::UNADJUSTED;
		}
		}
		assert(false);
		return BusinessDayConvention::FOLLOWING;
	}

	PeriodUnit period_unit_from_string(const char *s) const override final
	{
		PeriodUnit unit = PeriodUnit::TERMS;
		switch (*s) {
		case 'd':
		case 'D':
			unit = PeriodUnit::DAYS;
			break;
		case 'w':
		case 'W':
			unit = PeriodUnit::WEEKS;
			break;
		case 'm':
		case 'M':
			unit = PeriodUnit::MONTHS;
			break;
		case 'y':
		case 'Y':
			unit = PeriodUnit::YEARS;
			break;
		}
		return unit;
	}

	bool period_from_string(const char *periodName, Period *p) const override final
	{
		int n = strlen(periodName);
		if (n == 0 || n > 10) {
			fprintf(stderr, "invalid period %s\n", periodName);
			return false;
		}
		PeriodUnit unit = period_unit_from_string(&periodName[n - 1]);
		int x = 0;
		if (unit == PeriodUnit::TERMS) {
			if (strcasecmp(periodName, "O/N") == 0 || strcasecmp(periodName, "ON") == 0) {
				unit = PeriodUnit::DAYS;
				x = 1;
			} else if (strcasecmp(periodName, "T/N") == 0 || strcasecmp(periodName, "TN") == 0) {
				unit = PeriodUnit::DAYS;
				x = 2;
			} else {
				// fprintf(stderr, "unknown period type %s\n",
				// periodName);
				return false;
			}
		} else {
			// Exploit the fact that atoi ignores any
			// non numeric characters at the end
			x = atoi(periodName);
		}
		*p = Period(x, unit);
		return true;
	}

	DayCountFraction day_count_fraction_from_string(const char *value) const override final
	{
		auto iter = std::find_if(std::begin(day_count_fractions), std::end(day_count_fractions),
					 [value](const std::pair<const char *, DayCountFraction> &v) {
						 return strcmp(value, v.first) == 0;
					 });
		if (iter == std::end(day_count_fractions))
			return DayCountFraction::DAY_COUNT_FRACTION_UNSPECIFIED;
		else
			return iter->second;
	}

	const char *day_count_fraction_to_string(DayCountFraction value) const final
	{
		auto iter = std::find_if(
		    std::begin(day_count_fractions), std::end(day_count_fractions),
		    [value](const std::pair<const char *, DayCountFraction> &v) { return value == v.second; });
		if (iter == std::end(day_count_fractions))
			return "";
		else
			return iter->first;
	}

	Currency currency_from_string(const char *value) const override final
	{
		auto iter = std::find_if(
		    std::begin(currencies), std::end(currencies),
		    [value](const std::pair<const char *, Currency> &v) { return strcmp(value, v.first) == 0; });
		if (iter == std::end(currencies))
			return Currency::CURRENCY_UNSPECIFIED;
		else
			return iter->second;
	}

	CurveGroup curve_group_from_string(const char *value) const override final
	{
		auto iter = std::find_if(
		    std::begin(curve_groups), std::end(curve_groups),
		    [value](const std::pair<const char *, CurveGroup> &v) { return strcmp(value, v.first) == 0; });
		if (iter == std::end(curve_groups))
			return CurveGroup::CURVE_GROUP_UNSPECIFIED;
		else
			return iter->second;
	}

	MaturityGenerationRule maturity_generation_rule_from_string(const char *value) const override final
	{
		auto iter = std::find_if(std::begin(maturity_genrules), std::end(maturity_genrules),
					 [value](const std::pair<const char *, MaturityGenerationRule> &v) {
						 return strcmp(value, v.first) == 0;
					 });
		if (iter == std::end(maturity_genrules))
			return MaturityGenerationRule::MATURITY_GENERATION_RULE_DERIVE_FROM_INSTRUMENTS;
		else
			return iter->second;
	}

	const char *currency_to_string(Currency value) const override final
	{
		auto iter =
		    std::find_if(std::begin(currencies), std::end(currencies),
				 [value](const std::pair<const char *, Currency> &v) { return value == v.second; });
		if (iter == std::end(currencies))
			return "";
		else
			return iter->first;
	}

	IsdaIndex isda_index_from_string(const char *value) const override final
	{
		auto iter = std::find_if(
		    std::begin(isdaindices), std::end(isdaindices),
		    [value](const std::pair<const char *, IsdaIndex> &v) { return strcmp(value, v.first) == 0; });
		if (iter == std::end(isdaindices))
			return IsdaIndex::ISDA_INDEX_UNSPECIFIED;
		else
			return iter->second;
	}

	const char *isda_index_to_string(IsdaIndex value) const override final
	{
		auto iter =
		    std::find_if(std::begin(isdaindices), std::end(isdaindices),
				 [value](const std::pair<const char *, IsdaIndex> &v) { return value == v.second; });
		if (iter == std::end(isdaindices))
			return "";
		else
			return iter->first;
	}

	IndexFamily index_family_from_string(const char *value) const override final
	{
		auto iter = std::find_if(
		    std::begin(indices), std::end(indices),
		    [value](const std::pair<const char *, IndexFamily> &v) { return strcmp(value, v.first) == 0; });
		if (iter == std::end(indices))
			return IndexFamily::INDEX_FAMILY_UNSPECIFIED;
		else
			return iter->second;
	}

	const char *index_family_to_string(IndexFamily value) const override final
	{
		auto iter =
		    std::find_if(std::begin(indices), std::end(indices),
				 [value](const std::pair<const char *, IndexFamily> &v) { return value == v.second; });
		if (iter == std::end(indices))
			return "";
		else
			return iter->first;
	}

	Tenor tenor_from_period_unit_and_len(PeriodUnit unit, int value) const override final
	{
		if (value <= 0)
			return Tenor::TENOR_UNSPECIFIED;
		// Here we try to keep months and days in the lower range
		// of numbers as these period are more common
		switch (unit) {
		case PeriodUnit::DAYS: {
			if (value >= 1 && value <= 6)
				return (Tenor)value;
			else if (value == 7)
				return TENOR_1W;
			else
				return TENOR_UNSPECIFIED;
		}
		case PeriodUnit::WEEKS: {
			if (value >= 1 && value <= 3)
				return (Tenor)(6 + value);
			else if (value == 4)
				return TENOR_1M;
			else
				return TENOR_UNSPECIFIED;
		}
		case PeriodUnit::MONTHS: {
			if (value >= 1 && value <= 23)
				return (Tenor)(9 + value);
			else if (value == 24)
				return TENOR_2Y;
			else
				return TENOR_UNSPECIFIED;
		}
		case PeriodUnit::YEARS: {
			if (value >= 2 && value <= 65)
				return (Tenor)(31 + value);
			else if (value == 1)
				return TENOR_12M;
			else
				return TENOR_UNSPECIFIED;
		}
		case PeriodUnit::TERMS:
			return TENOR_1T;
		default:
			return TENOR_UNSPECIFIED;
		}
	}

	RollConvention roll_convention_from_string(const char *s) const override final
	{
		if (s == nullptr)
			return RollConvention::ROLL_CONVENTION_DEFAULT;
		if (strcmp(s, "IMM") == 0)
			return IMM;
		else if (strcmp(s, "EOM") == 0)
			return EOM;
		else if (strcmp(s, "IMMAUD") == 0)
			return IMMAUD;
		else if (strcmp(s, "IMMCAD") == 0)
			return IMMCAD;
		else if (strcmp(s, "IMMNZD") == 0)
			return IMMNZD;
		else if (strcmp(s, "NONE") == 0)
			return RollConvention::ROLL_CONVENTION_NONE;
		else {
			return RollConvention::ROLL_CONVENTION_DEFAULT;
		}
	}

	CompoundingMethod compounding_method_from_string(const char *s) const override final
	{
		if (s == nullptr || *s == 0)
			return CompoundingMethod::CM_NONE;
		else if (strcmp(s, "FLAT") == 0)
			return CompoundingMethod::FLAT;
		else if (strcmp(s, "STRAIGHT") == 0)
			return CompoundingMethod::STRAIGHT;
		else if (strcmp(s, "SPREAD_EXCLUSIVE") == 0)
			return CompoundingMethod::SPREAD_EXCLUSIVE;
		else
			return CompoundingMethod::CM_NONE;
	}
	const char *period_unit_to_string(PeriodUnit value) const override final
	{
		switch (value) {
		case PeriodUnit::DAYS:
			return "D";
		case PeriodUnit::WEEKS:
			return "W";
		case PeriodUnit::MONTHS:
			return "M";
		case PeriodUnit::YEARS:
			return "Y";
		case PeriodUnit::TERMS:
			return "T";
		default:
			return "";
		}
	}
	bool tenor_to_period_unit_and_multiplier(Tenor value, PeriodUnit *unit, int *period) const override final
	{
		if ((int)value <= 0 || (int)value > 97) {
			*unit = PeriodUnit::PERIOD_UNIT_UNRECOGNIZED;
			*period = 0;
			return false;
		}
		switch ((int)value) {
		case 97: {
			*unit = PeriodUnit::TERMS;
			*period = 1;
			break;
		}
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6: {
			*unit = PeriodUnit::DAYS;
			*period = (int)value;
			break;
		}
		case 6 + 1:
		case 6 + 2:
		case 6 + 3: {
			*unit = PeriodUnit::WEEKS;
			*period = (int)value - 6;
			break;
		}
		case 9 + 1:
		case 9 + 2:
		case 9 + 3:
		case 9 + 4:
		case 9 + 5:
		case 9 + 6:
		case 9 + 7:
		case 9 + 8:
		case 9 + 9:
		case 9 + 10:
		case 9 + 11:
		case 9 + 12:
		case 9 + 13:
		case 9 + 14:
		case 9 + 15:
		case 9 + 16:
		case 9 + 17:
		case 9 + 18:
		case 9 + 19:
		case 9 + 20:
		case 9 + 21:
		case 9 + 22:
		case 9 + 23: {
			*unit = PeriodUnit::MONTHS;
			*period = (int)value - 9;
			break;
		}
		default: {
			*unit = PeriodUnit::YEARS;
			*period = (int)value - 31;
			break;
		}
		}
		return true;
	}
	std::string tenor_to_string(Tenor tenor) const override final
	{
		PeriodUnit unit;
		int period;
		if (tenor_to_period_unit_and_multiplier(tenor, &unit, &period)) {
			char buf[10];
			snprintf(buf, sizeof buf, "%d%s", period, period_unit_to_string(unit));
			return std::string(buf);
		} else
			return std::string();
	}
	int tenor_to_days(Tenor tenor) const override final
	{
		int t = (int)tenor;
		if (t < Tenor::TENOR_1D || t > Tenor::TENOR_65Y)
			return 0;
		return map_tenor_to_days[t - 1];
	}

	InterpolatorType interpolator_type_from_string(const char *value) const
	{
		// mappings must be in alphabetic order
		auto iter = std::find_if(std::begin(interpolators), std::end(interpolators),
					 [value](const std::pair<const char *, InterpolatorType> &v) {
						 return strcmp(value, v.first) == 0;
					 });
		if (iter == std::end(interpolators)) {
			warn("cannot find interpolation method [%s], defaulting "
			     "to LINEAR\n",
			     value);
			return InterpolatorType::LINEAR;
		}
		return iter->second;
	}

	PricingCurveType pricing_curve_type_from_string(const char *s) const
	{
		if (s) {
			if (*s == 'F' || *s == 'f')
				return PricingCurveType::PRICING_CURVE_TYPE_FORWARD;
			else
				return PricingCurveType::PRICING_CURVE_TYPE_DISCOUNT;
		}
		return PRICING_CURVE_TYPE_FORWARD;
	}

	IRRateType rate_type_from_string(const char *s) const
	{
		switch (*s) {
		case 'Z':
		case 'z':
			return IRRateType::ZERO_RATE;
		case 'f':
		case 'F':
			return IRRateType::FORWARD_RATE;
		case 'd':
		case 'D':
			return IRRateType::DISCOUNT_FACTOR;
		default:
			return IRRateType::ZERO_RATE;
		}
	}
};

static ValueConverter default_converter;

const Converter *get_default_converter() { return &default_converter; }

/////////////////////////// Tests

int test_conversions()
{
	int failure_count = 0;

	if (get_default_converter()->business_center_from_string("USNY") != BusinessCenter::USNY)
		failure_count++;
	if (get_default_converter()->business_day_convention_from_string("MODFOLLOWING") !=
	    BusinessDayConvention::MODIFIED_FOLLOWING)
		failure_count++;
	if (get_default_converter()->day_count_fraction_from_string("ACT/360") != DayCountFraction::ACT_360)
		failure_count++;
	if (get_default_converter()->period_unit_from_string("Y") != PeriodUnit::YEARS)
		failure_count++;
	Period p;
	if (!get_default_converter()->period_from_string("12M", &p) || p.units() != PeriodUnit::MONTHS ||
	    p.length() != 12)
		failure_count++;
	for (int i = TENOR_1D; i <= TENOR_1T /* Max value */; i++) {
		Tenor t1 = (Tenor)i;
		p = Period::tenor_to_period(t1);
		Tenor t2 = get_default_converter()->tenor_from_period_unit_and_len(p.units(), p.length());
		if (t1 != t2)
			failure_count++;
	}
	if (get_default_converter()->roll_convention_from_string("IMM") != IMM)
		failure_count++;
	if (get_default_converter()->isda_index_from_string("USD-Federal Funds-H.15-OIS-COMPOUND") !=
	    IsdaIndex::USD_Federal_Funds_H_15_OIS_COMPOUND)
		failure_count++;
	if (get_default_converter()->currency_from_string("JPY") != Currency::JPY)
		failure_count++;
	if (get_default_converter()->index_family_from_string("PRIBOR") != PRIBOR)
		failure_count++;
	if (failure_count == 0)
		printf("Conversion Tests OK\n");
	else
		printf("Conversion Tests FAILED\n");
	return failure_count;
}

} // namespace redukti