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

#include <bootstrap.h>
#include <cashflow.h>
#include <logger.h>
#include <request_processor.h>
#include <valuation.h>

#include <inttypes.h>

#include <chrono>
#include <thread>

namespace redukti
{

// The request processor is just a gateway
// it passes on the requests to respective handlers
// The only requests that it directly handles are:
// a) Hello request
// b) Unknown request
class RequestProcessorImpl : public RequestProcessor
{
	private:
	std::unique_ptr<CurveBuilderService> bootstrapper_;
	std::unique_ptr<ValuationService> valuation_service_;
	void *shutdown_data_;
	void (*shutdown_func_)(void *);

	public:
	RequestProcessorImpl(std::unique_ptr<CurveBuilderService> bootstrapper,
			     std::unique_ptr<ValuationService> valuation_service)
	    : bootstrapper_(std::move(bootstrapper)), valuation_service_(std::move(valuation_service)),
	      shutdown_data_(nullptr), shutdown_func_(nullptr)
	{
	}
	Response *process(const Request *request, Response *response) override;

	void set_shutdown_handler(void *p, void (*funcptr)(void *))
	{
		shutdown_data_ = p;
		shutdown_func_ = funcptr;
	}

	private:
	Response *make_response(const ReplyHeader &reply_header, Response *response);
	Response *handle_hello_request(const Request *request, Response *response);
	Response *handle_shutdown_request(const Request *request, Response *response);
	Response *handle_unknown_request(const Request *request, Response *response);
};

std::unique_ptr<RequestProcessor> get_request_processor(std::unique_ptr<CurveBuilderService> bootstrapper,
							std::unique_ptr<ValuationService> valuation_service)
{
	return std::make_unique<RequestProcessorImpl>(std::move(bootstrapper), std::move(valuation_service));
}

std::unique_ptr<RequestProcessor> get_request_processor(const char *pricing_scipt)
{
	return std::make_unique<RequestProcessorImpl>(
	    get_curve_builder_service(std::string(pricing_scipt)),
	    get_valuation_service(get_default_index_service(), get_calendar_factory()));
}

Response *RequestProcessorImpl::handle_hello_request(const Request *request, Response *response)
{
	HelloReply *helloReply = response->mutable_hello_reply();
	helloReply->set_message(request->hello_request().name());
	ResponseHeader *header = response->mutable_header();
	header->set_response_code(StandardResponseCode::SRC_OK);
	int32_t delayfor = request->hello_request().delay_for();
	if (delayfor > 0) {
		inform("Sleeping for %d milliseconds\n", delayfor);
		std::this_thread::sleep_for(std::chrono::milliseconds(delayfor));
	}
	return response;
}

// This handler simply returns a success response
// The actual shutdown is initiated in the work handler
// TODO this needs to perform some sort of validation of the request
// Possibly password verification
Response *RequestProcessorImpl::handle_shutdown_request(const Request *request, Response *response)
{
	inform("Received shutdown request\n");
	ResponseHeader *header = response->mutable_header();
	if (shutdown_func_) {
		// Use a different thread to invoke shutdown
		header->set_response_code(StandardResponseCode::SRC_OK);
		std::thread shutdown_thread(shutdown_func_, shutdown_data_);
		shutdown_thread.detach();
	} else {
		warn("No shutdown hook registered, do not know how to shutdown\n");
		header->set_response_code(StandardResponseCode::SRC_ERROR);
		header->set_response_message("No shutdown hook registered, do not know how to shutdown");
	}
	return response;
}

Response *RequestProcessorImpl::handle_unknown_request(const Request *request, Response *response)
{
	ResponseHeader *header = response->mutable_header();
	header->set_response_code(StandardResponseCode::SRC_UNKNOWN_REQUEST);
	header->set_response_message("Unknown request type");
	return response;
}

Response *RequestProcessorImpl::make_response(const ReplyHeader &reply_header, Response *response)
{
	ResponseHeader *header = response->mutable_header();
	header->set_response_code(reply_header.response_code());
	header->set_response_sub_code(reply_header.response_sub_code());
	header->set_response_message(reply_header.response_message());
	return response;
}

Response *RequestProcessorImpl::process(const Request *request, Response *response)
{
	get_threadspecific_allocators()->reset();

	auto start = std::chrono::high_resolution_clock::now();
	switch (request->request_case()) {
	case Request::RequestCase::kHelloRequest: {
		response = handle_hello_request(request, response);
		break;
	}
	case Request::RequestCase::kShutdownRequest: {
		response = handle_shutdown_request(request, response);
		break;
	}
	case Request::RequestCase::kBootstrapCurvesRequest: {
		auto reply = bootstrapper_->handle_bootstrap_request(&request->bootstrap_curves_request(),
								     response->mutable_bootstrap_curves_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kCurveInterpolationRequest: {
		auto reply = valuation_service_->handle_curve_interpolation_request(
		    &request->curve_interpolation_request(), response->mutable_curve_interpolation_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kResetValuationServiceRequest: {
		auto reply = valuation_service_->handle_reset_valuation_service_request(
		    &request->reset_valuation_service_request(), response->mutable_reset_valuation_service_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kRegisterCurveDefinitionsRequest: {
		auto reply = valuation_service_->handle_register_curve_definitions_request(
		    &request->register_curve_definitions_request(),
		    response->mutable_register_curve_definitions_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kSetZeroCurvesRequest: {
		auto reply = valuation_service_->handle_set_zero_curves_request(
		    &request->set_zero_curves_request(), response->mutable_set_zero_curves_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kSetCurveMappingsRequest: {
		auto reply = valuation_service_->handle_set_curve_mappings_request(
		    &request->set_curve_mappings_request(), response->mutable_set_curve_mappings_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kSetFixingsRequest: {
		auto reply = valuation_service_->handle_set_fixings_request(&request->set_fixings_request(),
									    response->mutable_set_fixings_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kValuationRequest: {
		auto reply = valuation_service_->handle_valuation_request(&request->valuation_request(),
									  response->mutable_valuation_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kRegisterCalendarRequest: {
		auto reply = get_calendar_factory()->handle_register_calendar_request(
		    &request->register_calendar_request(), response->mutable_register_calendar_reply());
		response = make_response(reply->header(), response);
	}
	case Request::RequestCase ::kRegisterIndexDefinitionRequest: {
		auto reply = get_default_index_service()->handle_register_index_definition_request(
		    &request->register_index_definition_request(), response->mutable_register_index_definition_reply());
		response = make_response(reply->header(), response);
	}
	default: {
		response = handle_unknown_request(request, response);
		break;
	}
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration<int64_t, std::nano>(end - start);
	if (response && response->has_header()) {
		response->mutable_header()->set_elapsed_time(duration.count());
	}
	trace("Request completed in %lld nanosecs", duration.count());
	trace("Response %s\n", response->DebugString().c_str());
	return response;
}

} // namespace redukti