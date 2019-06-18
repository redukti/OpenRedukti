cmake_minimum_required(VERSION 3.3)
project(OpenRedukti)

cmake_policy(SET CMP0074 NEW)

enable_language(CXX)
enable_language(C)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

option(GRPC_SERVER  "Controls whether GRPC server will be built" OFF)
if (WIN32)
	option(USE_OPENBLAS "Controls whether OpenBLAS will be used" ON)
	option(USE_MKL      "Controls whether Intel MKL will be used" OFF)
endif(WIN32)

## On Windows
## If using vcpkg, put protoc on the path as follows
## set PATH=C:\work\vcpkg\installed\x64-windows\tools\protobuf;%PATH%
## If using custom build, pass this to Cmake:
## -DPROTOBUF_SRC_ROOT_FOLDER=c:\Software\protobuf371d

find_package(Protobuf REQUIRED)
#LAPACK and BLAS
if (WIN32)
	if (USE_OPENBLAS)
		find_package(OpenBLAS REQUIRED)
	endif(USE_OPENBLAS)
	if (USE_MKL)
		find_package(MKL REQUIRED)
	endif(USE_MKL)
else()
	find_package(LAPACK REQUIRED)
endif()
if (GRPC_SERVER)
	if (WIN32)
		find_package(gRPC CONFIG REQUIRED)
	else()
		find_package(GRPC REQUIRED)
	endif()
endif()

message( STATUS "Protobuf_LIBRARIES      : " ${Protobuf_LIBRARIES} )
message( STATUS "Protobuf_INCLUDE_DIRS   : " ${Protobuf_INCLUDE_DIRS} )
message( STATUS "LAPACK_LINKER_FLAGS     : " ${LAPACK_LINKER_FLAGS} )
message( STATUS "LAPACK_LIBRARIES        : " ${LAPACK_LIBRARIES} )
message( STATUS "BLAS_LIBRARIES          : " ${BLAS_LIBRARIES} )
message( STATUS "BLAS_LINKER_FLAGS       : " ${BLAS_LINKER_FLAGS} )
message( STATUS "OPENBLAS_LIBRARIES      : " ${OPENBLAS_LIBRARIES} )
message( STATUS "OPENBLAS_INCLUDE_DIRS   : " ${OPENBLAS_INCLUDE_DIRS} )
message( STATUS "MKL_LIBRARIES           : " ${MKL_LIBRARIES} )
message( STATUS "MKL_INCLUDE_DIRS        : " ${MKL_INCLUDE_DIRS} )
message( STATUS "GRPC_LIBRARIES          : " ${GRPC_LIBRARIES} )

# TODO Clang/C2 support  
if (MSVC AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	message( STATUS "MSVC selected" )
	set(CMAKE_C_FLAGS_DEBUG "/O1 /D_DEBUG /MDd /Zi /EHsc")
	set(CMAKE_C_FLAGS_RELEASE "/DNDEBUG /O2 /MD /EHsc")
	set(CMAKE_C_FLAGS_RELWITHDEBINFO "/DNDEBUG /O1 /MD /Zi /EHsc")
	set(CMAKE_C_FLAGS_MINSIZEREL "/DNDEBUG /O2 /MD /EHsc")

	set(CMAKE_CXX_FLAGS "/EHsc")
	set(CMAKE_CXX_FLAGS_DEBUG "/Od /D_DEBUG /MDd /Zi /EHsc")
	set(CMAKE_CXX_FLAGS_RELEASE "/DNDEBUG /O2 /MD /EHsc")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/DNDEBUG /O1 /MD /Zi /EHsc")
	set(CMAKE_CXX_FLAGS_MINSIZEREL "/DNDEBUG /O2 /MD /EHsc")

	# Disable signed/unsigned int conversion warnings.
	add_definitions("/wd4018")
	# Disable warning about using struct/class for the same symobl.
	add_definitions("/wd4099")
	# Disable performance warning about int-to-bool conversion.
	add_definitions("/wd4800")
	# Disable signed/unsigned int conversion warnings.
	add_definitions("/wd4018")
	# Disable warning about the insecurity of using "std::copy".
	add_definitions("/wd4996")
	add_definitions("/wd4291")
	add_definitions("/wd4624")
elseif (UNIX AND NOT APPLE)
	# -fsanitize=bounds -fsanitize=alignment -fsanitize=object-size
	# set(SANITIZER_FLAGS "-fsanitize=address")
	set(CMAKE_C_FLAGS_DEBUG "${SANITIZER_FLAGS} -fno-omit-frame-pointer -std=c99 -O0 -g3 -Wall -Wextra ${CXX_OPTIONS}")
	set(CMAKE_C_FLAGS_RELEASE "-std=c99 -O3 -Wall -Wextra ${CXX_OPTIONS}")
	set(CMAKE_C_FLAGS_RELWITHDEBINFO "${SANITIZER_FLAGS} -std=c99 -O1 -g3 -Wall -Wextra ${CXX_OPTIONS}")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wno-sign-compare -std=c++14 -fno-exceptions ${CXX_OPTIONS}")
	set(CMAKE_CXX_FLAGS_DEBUG "${SANITIZER_FLAGS} -fno-omit-frame-pointer -O0 -g3 -Wall -Wno-sign-compare -std=c++14 -fno-exceptions ${CXX_OPTIONS}")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${SANITIZER_FLAGS} -O1 -g3 -Wall -Wno-sign-compare -std=c++14 -fno-exceptions ${CXX_OPTIONS}")
elseif (APPLE)
	set(CMAKE_C_FLAGS "-std=c99 -O1 -g3 -Wall -Wextra")
	set(CMAKE_C_FLAGS_DEBUG "-std=c99 -O0 -g3 -Wall -Wextra")
	set(CMAKE_C_FLAGS_RELEASE "-std=c99 -O3 -Wall -Wextra")
	set(CMAKE_C_FLAGS_RELWITHDEBINFO "-std=c99 -O1 -g3 -Wall -Wextra")
	set(CMAKE_CXX_FLAGS "-O1 -g3 -Wall -Wno-sign-compare -std=c++14 -fno-exceptions")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wno-sign-compare -std=c++14 -fno-exceptions")
	set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -Wall -Wno-sign-compare -std=c++14 -fno-exceptions")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O1 -g3 -Wall -Wno-sign-compare -std=c++14 -fno-exceptions")
endif ()

include( cmake/ravi.cmake )
include( cmake/cminpack.cmake )
include( cmake/generators.cmake )

if (WIN32)
	# LAPACK defines  
	add_definitions(-DADD_)
	add_definitions(-DHAVE_LAPACK_CONFIG_H)
	add_definitions(-DLAPACK_COMPLEX_STRUCTURE)
	if (USE_OPENBLAS)
		add_definitions(-DUSE_OPENBLAS)
		include_directories("${OPENBLAS_INCLUDE_DIRS}")
		set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};${OPENBLAS_LIBRARIES}")
	elseif (USE_MKL)
		add_definitions(-DUSE_MKL)
		include_directories("${MKL_INCLUDE_DIRS}")
		set(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
	else()
		set(LAPACK_LIBRARIES "liblapack.lib;libblas.lib")
	endif()
elseif (UNIX AND NOT APPLE)
	# FIXME 
	include_directories(/usr/include/openblas)
	add_definitions(-DUSE_OPENBLAS)
	set(EXTRA_LIBRARIES m dl readline)
elseif (APPLE)
	include_directories("/System/Library/Frameworks")
	set(EXTRA_LIBRARIES m readline)
endif()

set(PROTO_FILES 
    proto/redukti/enums.proto
		proto/redukti/common.proto
		proto/redukti/index.proto
		proto/redukti/schedule.proto
		proto/redukti/cashflow.proto
		proto/redukti/curve.proto
		proto/redukti/shared.proto
		proto/redukti/bootstrap.proto
		proto/redukti/valuation.proto
		proto/redukti/infra.proto
		proto/redukti/services.proto
    )

protobuf_generate_cpp_new(PROTO_SRCS PROTO_HDRS ${PROTO_FILES})
if (GRPC_SERVER AND NOT WIN32)
	set(GRPC_PROTO_FILES
		proto/redukti/services.proto)
	protobuf_generate_grpc_cpp_new(services_grpc_srcs services_grpc_hdrs ${GRPC_PROTO_FILES})
	include_directories(${GRPC_INCLUDE_DIR})
	if (MSVC)
		add_definitions(-D_WIN32_WINNT=0x600)
	endif()
endif()

if (GRPC_SERVER AND WIN32)
	set(GRPC_PROTO_FILES
		proto/services.proto)

	get_filename_component(services_proto "proto/services.proto" ABSOLUTE)
	get_filename_component(services_proto_path "${services_proto}" PATH)

	message(STATUS "Using gRPC ${gRPC_VERSION}")
	set(_PROTOBUF_PROTOC $<TARGET_FILE:protobuf::protoc>)
	set(_PROTOBUF_LIBPROTOBUF protobuf::libprotobuf)
	set(_GRPC_GRPCPP_UNSECURE gRPC::grpc++_unsecure)
	set(_GRPC_CPP_PLUGIN_EXECUTABLE $<TARGET_FILE:gRPC::grpc_cpp_plugin>)

	set(GRPC_LIBRARIES gRPC::grpc++_unsecure) 

	set(services_proto_srcs "${CMAKE_CURRENT_BINARY_DIR}/services.pb.cc")
	set(services_proto_hdrs "${CMAKE_CURRENT_BINARY_DIR}/services.pb.h")
	set(services_grpc_srcs "${CMAKE_CURRENT_BINARY_DIR}/services.grpc.pb.cc")
	set(services_grpc_hdrs "${CMAKE_CURRENT_BINARY_DIR}/services.grpc.pb.h")
	#         --cpp_out "${CMAKE_CURRENT_BINARY_DIR}"

	add_custom_command(
      OUTPUT "${services_proto_srcs}" "${services_proto_hdrs}" "${services_grpc_srcs}" "${services_grpc_hdrs}"
      COMMAND ${_PROTOBUF_PROTOC}
      ARGS --grpc_out "${CMAKE_CURRENT_BINARY_DIR}"
        -I "${services_proto_path}"
        --plugin=protoc-gen-grpc="${_GRPC_CPP_PLUGIN_EXECUTABLE}"
        "${services_proto}"
	DEPENDS "${services_proto}")
	# Include generated *.pb.h files
	include_directories(${gRPC_INCLUDE_DIR})
	if (MSVC)
		add_definitions(-D_WIN32_WINNT=0x600)
	endif()
endif()

set(CMAKE_CXX_STANDARD 17)
include_directories("${PROJECT_SOURCE_DIR}/cpp")
include_directories(${Protobuf_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/redukti)

# FIXME
link_directories("${CMAKE_INSTALL_PREFIX}/lib")

set (SRCS 
	cpp/allocators.cpp
	cpp/autodiff.cpp
	cpp/bootstrap.cpp
	cpp/calendars.cpp
	cpp/cashflow.cpp
	cpp/cashflow_pricing.cpp
	cpp/converters.cpp
	cpp/curve.cpp
	cpp/datasource.cpp
	cpp/date.cpp
	cpp/dayfractions.cpp
	cpp/fixings.cpp
	cpp/hashtable.cpp
	cpp/index.cpp
	cpp/interpolators.cpp
	cpp/logger.cpp
	cpp/matrix.c
	cpp/raviapi.cpp
	cpp/request_processor.cpp
	cpp/schedule.cpp
	cpp/sensitivities.cpp
	cpp/statistics.cpp
	cpp/status.cpp
	cpp/valuation.cpp
	)
set (HEADERS 
	cpp/allocators.h
	cpp/autodiff.h
	cpp/bootstrap.h
	cpp/buffer.h
	cpp/calendars.h
	cpp/cashflow.h
	cpp/cashflow_pricing.h
	cpp/converters.h
	cpp/curve.h
	cpp/datasource.h
	cpp/date.h
	cpp/dayfractions.h
	cpp/fixings.h
	cpp/hashtable.h
	cpp/index.h
	cpp/interpolators.h
	cpp/linalg.h
	cpp/logger.h
	cpp/matrix.h
	cpp/request_processor.h
	cpp/schedule.h
	cpp/statistics.h
	cpp/status.h
	cpp/timeseries.h
	cpp/valuation.h
	)
set (INTERNAL_HEADERS
	cpp/internal/raviapi_internal.h 
	cpp/internal/autodiff_helpers.h
	cpp/internal/cashflow_internal.h
	cpp/internal/cashflow_pricing_internal.h
	)
    
add_definitions(-DPRICING_SCRIPT="${PROJECT_SOURCE_DIR}/scripts/pricing.lua")

source_group("OpenRedukti Sources" FILES ${SRCS})
source_group("OpenRedukti Headers" FILES ${HEADERS} ${INTERNAL_HEADERS})
source_group("Protobuf Headers" FILES ${PROTO_HDRS} ${services_proto_hdrs})
if (GRPC_SERVER)
	source_group("GRPC Headers" FILES ${services_grpc_hdrs})
	source_group("GRPC Sources" FILES ${services_grpc_srcs})
endif()

if (MSVC)
	set(LIBRARY_TYPE STATIC)
else()
	set(LIBRARY_TYPE SHARED)
endif()

#Main library
add_library(openredukti ${LIBRARY_TYPE} 
	${HEADERS} ${INTERNAL_HEADERS} ${PROTO_HDRS} ${RAVI_HEADERS}
	${SRCS} ${PROTO_SRCS} ${RAVI_CORE_SRCS} ${RAVI_LIB_SRCS} ${RAVI_JIT_SRCS} 
	${DMR_C_SRCS} ${DMR_C_JIT_SRCS} ${cminpack_srcs} ${cminpack_hdrs})
target_link_libraries(openredukti ${Protobuf_LIBRARIES} ${LUA_LIBRARIES} 
	${SYSTEM_LIBS} ${LAPACK_LIBRARIES} ${EXTRA_LIBRARIES}) 
if(APPLE)
	set_target_properties(openredukti PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

if (GRPC_SERVER)
	set(REDUKTISERV reduktiserv)
	add_executable(${REDUKTISERV} grpc-server/grpcmain.cpp ${services_proto_srcs} ${services_grpc_srcs})
	target_link_libraries(${REDUKTISERV} openredukti ${GRPC_LIBRARIES} ${Protobuf_LIBRARIES})
endif()

add_executable(runtests cpp/runtests.cpp)
target_link_libraries(runtests openredukti) 

add_executable(dylan ${RAVI_CMD_SRCS})
target_link_libraries(dylan openredukti) 

# Ravi VSCode Debug adapter
add_executable(ravidebug ${RAVI_DEBUGGER_HEADERS} ${RAVI_DEBUGGER_SRCS})
target_link_libraries(ravidebug openredukti)

add_test(RunTests runtests)

install(FILES ${HEADERS} ${PROTO_HDRS}
		DESTINATION include/redukti)
install(FILES ${PROTO_FILES} ${GRPC_PROTO_FILES}
	DESTINATION proto/redukti)
install(FILES ${LUA_HEADERS} 
	DESTINATION include/ravi)
install(DIRECTORY scripts/
	DESTINATION scripts)
install(TARGETS openredukti dylan ravidebug ${REDUKTISERV}
	RUNTIME DESTINATION bin
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib)
