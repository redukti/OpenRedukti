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

#ifndef _REDUKTI_VALUATION_H
#define _REDUKTI_VALUATION_H

#include <valuation.pb.h>

namespace redukti
{

class ValuationService
{
	public:
	virtual ~ValuationService() {}
	virtual CurveInterpolationReply *
	handle_curve_interpolation_request(google::protobuf::Arena *arena,
					   const CurveInterpolationRequest *request) = 0;
	virtual SetCurveMappingsReply *handle_set_curve_mappings_request(google::protobuf::Arena *arena,
									 const SetCurveMappingsRequest *request) = 0;
	virtual SetZeroCurvesReply *handle_set_zero_curves_request(google::protobuf::Arena *arena,
								   const SetZeroCurvesRequest *request) = 0;
	virtual RegisterCurveDefinitionsReply *
	handle_register_curve_definitions_request(google::protobuf::Arena *arena,
						  const RegisterCurveDefinitionsRequest *request) = 0;
	virtual SetFixingsReply *handle_set_fixings_request(google::protobuf::Arena *arena,
							    const SetFixingsRequest *request) = 0;
	virtual ValuationReply *handle_valuation_request(google::protobuf::Arena *arena,
							 const ValuationRequest *request) = 0;
};
} // namespace redukti

#endif