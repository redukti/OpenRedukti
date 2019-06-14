#include <grpcpp/grpcpp.h>

#include "services.grpc.pb.h"

#include <request_processor.h>
#include <bootstrap.h>
#include <cashflow.h>
#include <valuation.h>
#include <logger.h>

#include <inttypes.h>

#include <chrono>
#include <thread>
#include <iostream>
#include <memory>
#include <string>

using grpc::Server;
using grpc::ServerBuilder;
using grpc::ServerContext;
using grpc::Status;

namespace redukti
{

// Logic and data behind the server's behavior.
class OpenReduktiServiceImpl final : public OpenReduktiServices::Service
{
	private:
	std::unique_ptr<RequestProcessor> request_processor_;

	public:
	OpenReduktiServiceImpl(std::unique_ptr<RequestProcessor> request_processor)
	    : request_processor_(std::move(request_processor))
	{
	}

	Status serve(ServerContext *context, const Request *request, Response *response) override { 
		auto response = request_processor_->process(request, response);		
		return Status::OK; 
	}
};

void RunServer(OpenReduktiServiceImpl service)
{
	std::string server_address("0.0.0.0:50051");

	ServerBuilder builder;
	// Listen on the given address without any authentication mechanism.
	builder.AddListeningPort(server_address, grpc::InsecureServerCredentials());
	// Register "service" as the instance through which we'll communicate with
	// clients. In this case it corresponds to an *synchronous* service.
	builder.RegisterService(&service);
	// Finally assemble the server.
	std::unique_ptr<Server> server(builder.BuildAndStart());
	std::cout << "Server listening on " << server_address << std::endl;

	// Wait for the server to shutdown. Note that some other thread must be
	// responsible for shutting down the server for this call to ever return.
	server->Wait();
}

} // namespace redukti

using namespace redukti;

static void usage(const char *progname)
{
	fprintf(stderr,
		"%s: [-a <address, default 127.0.0.1>] [-p <port, default 9001>] "
		"-s <scriptpath> "
		"[-loglevel <d|t>]\n",
		progname);
}


int
main(int argc, char **argv)
{
	const char *progname = argv[0];
	char server_address[80];
	int port = 9001;
	char script[1024];

	/* Default server address */
	strcpy(server_address, "127.0.0.1");
	script[0] = 0;

	for (int i = 1; i + 1 < argc; i += 2) {
		const char *opt = argv[i];
		const char *arg = argv[i + 1];
		if (strcmp(opt, "-a") == 0) {
			strncpy(server_address, arg, sizeof server_address);
			server_address[sizeof server_address - 1] = 0; /* safety */
		} else if (strcmp(opt, "-p") == 0) {
			port = atoi(arg);
		} else if (strcmp(opt, "-s") == 0) {
			strncpy(script, arg, sizeof script);
			script[sizeof script - 1] = 0; /* safety */
		} else if (strcmp(opt, "-loglevel") == 0) {
			if (arg[0] && (arg[0] == 'd' || arg[0] == 'D')) {
				Redukti_log_mask |= LOG_DEBUG;
				inform("Log level set to DEBUG\n");
			} else if (arg[0] && (arg[0] == 't' || arg[0] == 'T')) {
				Redukti_log_mask |= (LOG_DEBUG | LOG_TRACE);
				inform("Log level set to TRACE\n");
			}
		} else {
			error("Error: unrecognized option '%s'\n", opt);
			usage(progname);
			return 1;
		}
	}

	std::unique_ptr<CurveBuilderService> bootstrapper = get_curve_builder_service(script);
	IndexService *index_service = get_default_index_service();
	CalendarService *calendar_service = get_calendar_factory();
	std::unique_ptr<ValuationService> valuation_service = get_valuation_service(index_service, calendar_service);

	redukti::OpenReduktiServiceImpl service;
	redukti::RunServer(service);

	return 0;
}