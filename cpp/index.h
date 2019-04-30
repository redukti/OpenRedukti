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

#ifndef _REDUKTI_INDEX_H
#define _REDUKTI_INDEX_H

#include <enums.pb.h>

#include <calendars.h>
#include <date.h>
#include <dayfractions.h>

#include <inttypes.h>
#include <memory>

namespace redukti
{
class IndexDefinition;

// Unique identifier for an index
typedef uint32_t IndexId;

// Makes a unique identifier from the give ISDA index identifier and
// tenor
inline IndexId make_index_id(IsdaIndex isda_index, Tenor tenor)
{
	assert(IsdaIndex_IsValid(isda_index));
	assert(Tenor_IsValid(tenor));
	return isda_index << 8 | tenor;
}

// Base type for all indices
class Index
{
	public:
	virtual ~Index() {}
	virtual IndexId id() const = 0;
};

// An interest rate index representation. A requirement of
// OpenRedukti is that an each unique IndexId should map to one
// InterestRateIndex instance - as this allows the code to freely
// reference such instances without fear of the reference going away.
// Additionally a requirement is that the instance is immutable.
class InterestRateIndex : public Index
{
	public:
	virtual ~InterestRateIndex() {}
	virtual Currency currency() const = 0;
	virtual IndexFamily family() const = 0;
	virtual Tenor tenor() const = 0;
	virtual IsdaIndex isda_index() const = 0;

	// Given a fixing date, calculate the value date
	// by applying the calendars, day conventions associated
	// with the index
	virtual Date value_date(Date fixing_date) const = 0;
	// Given a value date, calculate the fixing date
	// by applying the calendars, day conventions associated
	// with the index
	virtual Date fixing_date(Date accrual_start_date) const = 0;
	// Given a value date calculate the maturity date
	// Appropriate calendars, day conventions and EOM rules
	// must be applied
	virtual Date maturity_date(Date value_date) const = 0;
	virtual bool is_valid_fixing_date(Date date) const = 0;
	virtual const Calendar *fixing_calendar() const = 0;
	virtual const DayFraction *day_fraction() const = 0;
	virtual BusinessDayConvention day_convention() const = 0;
};

// The IndexService is responsible for returning instances of InterestRateIndex.
// Note that the index service must ensure the following:
// a) There will only ever be one instance of an InterestRateIndex for a given
//    IndexId.
// b) Clients must be free to hold on to references to such instances without
//    fear of them going out of scope. So essentially these instances can only
//    be deleted at system shutdown.
// c) An InterestRateIndex instance must be immutable.
class IndexService
{
	public:
	virtual ~IndexService() {}

	// Adds a definition for use as a template for generating instances of
	// InterestRateIndex
	virtual bool register_index(const IndexDefinition &definition) = 0;

	// Obtains an instance of IntrestRateIndex - must return an existing
	// instance if already defined
	virtual InterestRateIndex *get_index(IsdaIndex isda_index, Tenor tenor) = 0;

	// Obtains an instance of IntrestRateIndex - must return an existing
	// instance if already defined
	virtual InterestRateIndex *get_index(Currency currency, IndexFamily index_family, Tenor tenor) = 0;
};

extern IndexService *get_default_index_service();

extern int test_index();

} // namespace redukti

#endif