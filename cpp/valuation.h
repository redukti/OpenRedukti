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

#include <index.h>
#include <calendars.h>

#include <valuation.pb.h>

namespace redukti
{

class ValuationService
{
	public:
	virtual ~ValuationService() {}
	virtual CurveInterpolationReply *
	handle_curve_interpolation_request(
					   const CurveInterpolationRequest *request, CurveInterpolationReply *reply) = 0;
	virtual SetCurveMappingsReply *handle_set_curve_mappings_request(
									 const SetCurveMappingsRequest *request, SetCurveMappingsReply *reply) = 0;
	virtual SetZeroCurvesReply *handle_set_zero_curves_request(
								   const SetZeroCurvesRequest *request, SetZeroCurvesReply *reply) = 0;
	virtual RegisterCurveDefinitionsReply *
	handle_register_curve_definitions_request(
						  const RegisterCurveDefinitionsRequest *request, RegisterCurveDefinitionsReply *reply) = 0;
	virtual SetFixingsReply *handle_set_fixings_request(
							    const SetFixingsRequest *request, SetFixingsReply *reply) = 0;
	virtual ValuationReply *handle_valuation_request(
							 const ValuationRequest *request, ValuationReply *reply) = 0;
	virtual ResetValuationServiceReply *
	handle_reset_valuation_service_request(
					       const ResetValuationServiceRequest *request, ResetValuationServiceReply *reply) = 0;
};

std::unique_ptr<ValuationService> get_valuation_service(IndexService *index_service, CalendarService *calendar_service);

} // namespace redukti

#endif