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

#ifndef _REDUKTI_REQUEST_PROCESSOR_H
#define _REDUKTI_REQUEST_PROCESSOR_H

#include <valuation.h>
#include <bootstrap.h>

#include <services.pb.h>

#include <memory>

namespace redukti
{

class RequestProcessor
{
	public:
	virtual ~RequestProcessor() {}
	virtual Response *process(const Request *request, Response *response) = 0;
	virtual void set_shutdown_handler(void *p, void (*funcptr)(void *)) = 0;
};

// Note that the RequestProcessor will take ownership
// of the CurveBuilderService and ValuationService
std::unique_ptr<RequestProcessor> get_request_processor(
	std::unique_ptr<CurveBuilderService> bootstrapper,
	std::unique_ptr<ValuationService> valuation_service);

} // namespace

#endif