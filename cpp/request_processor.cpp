#include <request_processor.h>
#include <bootstrap.h>
#include <valuation.h>
#include <cashflow.h>
#include <logger.h>

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

	public:
	RequestProcessorImpl(std::unique_ptr<CurveBuilderService> bootstrapper,
			     std::unique_ptr<ValuationService> valuation_service)
	    : bootstrapper_(std::move(bootstrapper)), valuation_service_(std::move(valuation_service))
	{
	}
	Response *process(const Request *request, Response* response) override;

	private:
	Response *make_response(const ReplyHeader &reply_header, Response *response);
	Response *handle_hello_request(const Request *request, Response *response);
	Response *handle_shutdown_request(const Request *request, Response *response);
	Response *handle_unknown_request(const Request *request, Response *response);
};

std::unique_ptr<RequestProcessor> get_request_processor(
	std::unique_ptr<CurveBuilderService> bootstrapper,
							std::unique_ptr<ValuationService> valuation_service)
{
	return std::make_unique<RequestProcessorImpl>(bootstrapper, valuation_service);
}

Response *RequestProcessorImpl::handle_hello_request(const Request *request, Response *response)
{
	HelloReply *helloReply = response->mutable_hello_reply();
	helloReply->set_message(request->hello_request().name());
	ResponseHeader *header = response->mutable_header();
	header->set_response_code(StandardResponseCode::SRC_OK);
	response->set_allocated_header(header);
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
	header->set_response_code(StandardResponseCode::SRC_OK);
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

	Response *response = nullptr;
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
		    &request->register_curve_definitions_request(), response->mutable_register_curve_definitions_reply());
		response = make_response(reply->header(), response);
		break;
	}
	case Request::RequestCase::kSetZeroCurvesRequest: {
		auto reply =
		    valuation_service_->handle_set_zero_curves_request(&request->set_zero_curves_request(), response->set_zero_curves_reply());
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