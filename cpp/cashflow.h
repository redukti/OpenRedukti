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
#ifndef _REDUKTI_CASHFLOW_H_
#define _REDUKTI_CASHFLOW_H_

#include <enums.pb.h>

#include <allocators.h>
#include <autodiff.h>
#include <curve.h>
#include <date.h>
#include <fixings.h>
#include <hashtable.h>
#include <index.h>
#include <status.h>

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <string>

namespace redukti
{

// We separate out the concept of a Cashflow
// definition (covered here) versus the valuation of
// cashflows covered in cashflow_pricing.h.

// The protobuf definition for a Cashflow Collection
class CFCollection;

// We need a way to refer to logical curve types
// without having to reference real curves - the PricingCurve
// helps us do that. Each PricingCurve instance represents
// a logical identifier for a curve that will be resolved when
// pricing via a CurveProvider implementation.
class PricingCurve
{
      private:
	uint32_t id_;

      public:
	// Default of 0 is okay as it maps to unspecified
	// values component wise
	PricingCurve() : id_(0) {}
	PricingCurve(PricingCurveType type, Currency currency, IndexFamily index_family = INDEX_FAMILY_UNSPECIFIED,
		     Tenor tenor = TENOR_UNSPECIFIED)
	    : id_((uint32_t)type | ((uint32_t)currency << 8) | ((uint32_t)index_family << 16) | ((uint32_t)tenor << 24))
	{
		assert(Currency_IsValid(currency));
		assert(IndexFamily_IsValid(index_family));
		assert(Tenor_IsValid(tenor) && tenor <= TENOR_12M);
		assert(PricingCurveType_IsValid(type));
	}
	explicit PricingCurve(uint32_t id) : id_(id)
	{
		auto ccy = ((id_ >> 8) & 0xFF);
		auto family = ((id_ >> 16) & 0xFF);
		auto t = ((id_ >> 24) & 0xFF);
		auto ct = (id_ & 0xFF);
		if (!Currency_IsValid(ccy) || !PricingCurveType_IsValid(ct) || !IndexFamily_IsValid(family) ||
		    !Tenor_IsValid(t) || t > TENOR_12M) {
			id_ = 0;
		}
	}

	PricingCurve(const PricingCurve &) = default;
	PricingCurve &operator=(const PricingCurve &) = default;

	Currency currency() const
	{
		auto ccy = ((id_ >> 8) & 0xFF);
		return (Currency)ccy;
	}
	IndexFamily index_family() const
	{
		auto family = ((id_ >> 16) & 0xFF);
		return (IndexFamily)family;
	}
	Tenor tenor() const
	{
		auto t = ((id_ >> 24) & 0xFF);
		return (Tenor)t;
	}
	PricingCurveType curve_type() const
	{
		auto ct = (id_ & 0xFF);
		return (PricingCurveType)ct;
	}
	uint32_t id() const { return id_; }
	bool is_valid() const { return id_ != 0; }
	// Ordering is not meaningful - its purpose is to allow
	// insertion into containers
	bool operator<(const PricingCurve &c2) const { return id_ < c2.id_; }
	bool operator==(const PricingCurve &c2) const { return id_ == c2.id_; }
	bool operator!=(const PricingCurve &c2) const { return id_ != c2.id_; }
	// Get a string representation of the PricingCurve
	// Note that this is an expensive operation so use only for
	// debugging
	std::string name() const;
};

// Create a PricingCurve with specified type, and currency, index family
// and tenor taken from the supplied curve Id.
extern PricingCurve make_pricing_curve(PricingCurveType type, CurveId id);

// When generating cashflows we do not know what actual curves will
// be used - and whether the forward and discount curves map to the same
// curve or different curves, or whether different tenor curves map to
// different curves or the same curve. The CurveMapper allows the caller
// to provide a mapping to the desired 'logical' curve. The mapping is
// logical so that given a logical curve id, another function must obtain
// an instance of the real curve.
class CurveMapper
{
      public:
	virtual ~CurveMapper() {}
	virtual PricingCurve map_index_tenor(PricingCurveType curve_type, Currency currency,
					     IndexFamily family = IndexFamily::INDEX_FAMILY_UNSPECIFIED,
					     Tenor tenor = Tenor::TENOR_UNSPECIFIED) const = 0;
};

class ValuationContext
{
      public:
	virtual ~ValuationContext() {}
	virtual Date evaluation_date() const = 0;
	virtual Date payment_cutoff_date() const = 0;
	virtual int derivative_order() const = 0;
	// include today's fixing (e.g. eod)
	// If false then the curve will be used to determine the
	// rate. The bootstrapper requires this to be false so we set
	// the default value to false
	virtual bool include_todays_fixing() const = 0;
	// Retrieve a fixing.
	// If the fixing date is < evaluation date then the absence of a fixing
	// will be an error reported via status. If the fixing date is ==
	// evaluation date then a missing fixing is not treated as error -
	// instead the method will return false; in all cases a true return
	// value indicates that the fixing was found and is set
	virtual bool get_fixing(IndexId fixing_key, Date fixing_date, double &fixing, StatusCode &status) const = 0;
};

class Cashflows;

// Converts the CFCollection to internal cashflow format
extern Cashflows *construct_cashflows(RegionAllocator *A, const CFCollection *cfcollection, const ValuationContext &ctx,
				      const CurveMapper *curve_mapper);

extern int test_pricing();

} // namespace redukti

#endif
