#include <grpcpp/grpcpp.h>

#include "services.grpc.pb.h"

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
	Status serve(ServerContext *context, const Request *request, Response *response) override { return Status::OK; }
};

void RunServer()
{
	std::string server_address("0.0.0.0:50051");
	OpenReduktiServiceImpl service;

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

int main(int argc, char **argv)
{
	redukti::RunServer();

	return 0;
}