#ifndef _REDUKTI_REQUEST_PROCESSOR_H
#define _REDUKTI_REQUEST_PROCESSOR_H

#include <valuation.h>
#include <bootstrap.h>

#include <services.pb.h>

#include <memory>

namespace redukti
{

/* This is the API that the request processor must implement.
All messages allocated by the processor using the supplied Arena will
be automatically deallocated when the work is deallocated. In particular
the response object should be Arena allocated.
*/
class RequestProcessor
{
	public:
	virtual ~RequestProcessor() {}
	virtual Response *process(const Request *request, Response *response) = 0;
};

std::unique_ptr<RequestProcessor> get_request_processor(
	std::unique_ptr<CurveBuilderService> bootstrapper,
	std::unique_ptr<ValuationService> valuation_service);

} // namespace

#endif