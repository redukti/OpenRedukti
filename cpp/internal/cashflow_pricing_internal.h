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

#ifndef _REDUKTI_CASHFLOW_PRICING_INTERNAL_H
#define _REDUKTI_CASHFLOW_PRICING_INTERNAL_H

#include <autodiff.h>
#include <cashflow.h>
#include <logger.h>

#include <algorithm>
#include <array>
#include <map>
#include <memory>

namespace redukti
{

class SimpleCurveProvider : public CurveProvider
{
	private:
	std::map<PricingCurve, const CurveReference *> mappings_;

	public:
	void add_mapping(PricingCurve curveid, const CurveReference *ref) { mappings_.insert({curveid, ref}); }
	const CurveReference *get_curve(PricingCurve curve) const override final
	{
		auto iter = mappings_.find(curve);
		if (iter == mappings_.end()) {
			warn("No curve found for %s\n", curve.name().c_str());
			return nullptr;
		}
		return iter->second;
	}
};

// Allocate sensitivities to curve pillars and aggregate.
// The input sensitivities are provided in the auto diff
// variable data; the VarResolvers object provides the mapping
// from variable id to a cashflow date. The aggregated sensitivities
// are stored in container.
extern void aggregate_sensitivity(FixedRegionAllocator *allocator, const VariableResolver *discount_resolver,
				  const CurveReference *dc, const VariableResolver *forward1_resolver,
				  const CurveReference *f1, const VariableResolver *forward2_resolver,
				  const CurveReference *f2, Sensitivities *container, redukti_adouble_t &data,
				  StatusCode &status);

} // namespace redukti

#endif