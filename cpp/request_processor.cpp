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
	Response *make_response(google::protobuf::Arena *arena, const ReplyHeader &reply_header);
	Response *handle_hello_request(const Request *request, Response *response);
	Response *handle_shutdown_request(google::protobuf::Arena *arena, const Request *request);
	Response *handle_unknown_request(google::protobuf::Arena *arena, const Request *request);
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
Response *RequestProcessorImpl::handle_shutdown_request(google::protobuf::Arena *arena, const Request *request)
{
	Response *response = google::protobuf::Arena::CreateMessage<Response>(arena);
	inform("Received shutdown request\n");
	ResponseHeader *header = google::protobuf::Arena::CreateMessage<ResponseHeader>(arena);
	header->set_response_code(StandardResponseCode::SRC_OK);
	response->set_allocated_header(header);
	return response;
}

Response *RequestProcessorImpl::handle_unknown_request(google::protobuf::Arena *arena, const Request *request)
{
	Response *response = google::protobuf::Arena::CreateMessage<Response>(arena);
	ResponseHeader *header = google::protobuf::Arena::CreateMessage<ResponseHeader>(arena);
	header->set_response_code(StandardResponseCode::SRC_UNKNOWN_REQUEST);
	header->set_response_message("MyCCP.Risk: Unknown request type");
	response->set_allocated_header(header);
	return response;
}

Response *RequestProcessorImpl::make_response(google::protobuf::Arena *arena, const ReplyHeader &reply_header)
{
	Response *response = google::protobuf::Arena::CreateMessage<Response>(arena);
	ResponseHeader *header = google::protobuf::Arena::CreateMessage<ResponseHeader>(arena);
	response->set_allocated_header(header);
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
		response = handle_hello_request(arena, request);
		break;
	}
	case Request::RequestCase::kShutdownRequest: {
		response = handle_shutdown_request(arena, request);
		break;
	}
	case Request::RequestCase::kBootstrapCurvesRequest: {
		auto reply = bootstrapper_->handle_bootstrap_request(arena, &request->bootstrap_curves_request());
		response = make_response(arena, reply->header());
		response->set_allocated_bootstrap_curves_reply(reply);
		break;
	}
	case Request::RequestCase::kCurveInterpolationRequest: {
		auto reply = valuation_service_->handle_curve_interpolation_request(
		    arena, &request->curve_interpolation_request());
		response = make_response(arena, reply->header());
		response->set_allocated_curve_interpolation_reply(reply);
		break;
	}
	case Request::RequestCase::kResetValuationServiceRequest: {
		auto reply = valuation_service_->handle_reset_valuation_service_request(
		    arena, &request->reset_valuation_service_request());
		response = make_response(arena, reply->header());
		response->set_allocated_reset_valuation_service_reply(reply);
		break;
	}
	case Request::RequestCase::kRegisterCurveDefinitionsRequest: {
		auto reply = valuation_service_->handle_register_curve_definitions_request(
		    arena, &request->register_curve_definitions_request());
		response = make_response(arena, reply->header());
		response->set_allocated_register_curve_definitions_reply(reply);
		break;
	}
	case Request::RequestCase::kSetZeroCurvesRequest: {
		auto reply =
		    valuation_service_->handle_set_zero_curves_request(arena, &request->set_zero_curves_request());
		response = make_response(arena, reply->header());
		response->set_allocated_set_zero_curves_reply(reply);
		break;
	}
	case Request::RequestCase::kSetCurveMappingsRequest: {
		auto reply = valuation_service_->handle_set_curve_mappings_request(
		    arena, &request->set_curve_mappings_request());
		response = make_response(arena, reply->header());
		response->set_allocated_set_curve_mappings_reply(reply);
		break;
	}
	case Request::RequestCase::kSetFixingsRequest: {
		auto reply = valuation_service_->handle_set_fixings_request(arena, &request->set_fixings_request());
		response = make_response(arena, reply->header());
		response->set_allocated_set_fixings_reply(reply);
		break;
	}
	case Request::RequestCase::kValuationRequest: {
		auto reply = valuation_service_->handle_valuation_request(arena, &request->valuation_request());
		response = make_response(arena, reply->header());
		response->set_allocated_valuation_reply(reply);
		break;
	}
	default: {
		response = handle_unknown_request(arena, request);
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