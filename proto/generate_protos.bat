mkdir ..\generated
mkdir ..\generated\cpp
protoc.exe -I. --cpp_out=../generated/cpp enums.proto
protoc.exe -I. --cpp_out=../generated/cpp index.proto
protoc.exe -I. --cpp_out=../generated/cpp schedule.proto
protoc.exe -I. --cpp_out=../generated/cpp cashflow.proto
protoc.exe -I. --cpp_out=../generated/cpp curve.proto
protoc.exe -I. --cpp_out=../generated/cpp shared.proto
protoc.exe -I. --cpp_out=../generated/cpp bootstrap.proto
protoc.exe -I. --cpp_out=../generated/cpp valuation.proto

