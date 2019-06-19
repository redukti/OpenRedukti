# We have to use a custom version the default one in FindProtobuf isn't
# flexible enough to handle relative paths. Unfortunately relative paths are
# needed to generate Python versions with proper package structure
# Our goal here is to ensure that imports of protobuf files include
# the redukti folder; this translates in Python to the redukti package.
# FIXME there is some hardcoding here - this should be parameterised
# Following variables must be set before calling the functions her:
# Protobuf_PROTOC_EXECUTABLE
function(PROTOBUF_GENERATE_CPP_NEW SRCS HDRS)
    if(NOT ARGN)
        message(SEND_ERROR "Error: PROTOBUF_GENERATE_GRPC_CPP() called without any proto files")
        return()
    endif()

    get_filename_component(PROTO_PATH proto ABSOLUTE)

    set(${SRCS})
    set(${HDRS})
    foreach(FIL ${ARGN})
        get_filename_component(ABS_FIL ${FIL} ABSOLUTE)
        get_filename_component(FIL_WE ${FIL} NAME_WE)

        list(APPEND ${SRCS} "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.pb.cc")
        list(APPEND ${HDRS} "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.pb.h")

        add_custom_command(
                OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.pb.cc"
                "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.pb.h"
                COMMAND  ${Protobuf_PROTOC_EXECUTABLE}
                ARGS --cpp_out=${CMAKE_CURRENT_BINARY_DIR}
                -I ${PROTO_PATH} ${ABS_FIL}
                DEPENDS ${ABS_FIL} ${Protobuf_PROTOC_EXECUTABLE}
                COMMENT "Running C++ protocol buffer compiler on ${FIL}"
                VERBATIM)
    endforeach()

    set_source_files_properties(${${SRCS}} ${${HDRS}} PROPERTIES GENERATED TRUE)
    set(${SRCS} ${${SRCS}} PARENT_SCOPE)
    set(${HDRS} ${${HDRS}} PARENT_SCOPE)
endfunction()

# We have to use a custom version the default one in FindProtobuf isn't
# flexible enough to handle relative paths. Unfortunately relative paths are
# needed to generate Python versions with proper package structure
# Our goal here is to ensure that imports of protobuf files include
# the redukti folder; this translates in Python to the redukti package.
# FIXME there is some hardcoding here - this should be parameterised
# Following variables must be set before calling the functions her:
# Protobuf_PROTOC_EXECUTABLE
# GRPC_CPP_PLUGIN
function(PROTOBUF_GENERATE_GRPC_CPP_NEW SRCS HDRS)
    if(NOT ARGN)
        message(SEND_ERROR "Error: PROTOBUF_GENERATE_GRPC_CPP_NEW() called without any proto files")
        return()
    endif()

    get_filename_component(PROTO_PATH proto ABSOLUTE)

    set(${SRCS})
    set(${HDRS})
    foreach(FIL ${ARGN})
        get_filename_component(ABS_FIL ${FIL} ABSOLUTE)
        get_filename_component(FIL_WE ${FIL} NAME_WE)

        list(APPEND ${SRCS} "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.grpc.pb.cc")
        list(APPEND ${HDRS} "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.prpc.pb.h")

        add_custom_command(
                OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.grpc.pb.cc"
                "${CMAKE_CURRENT_BINARY_DIR}/redukti/${FIL_WE}.grpc.pb.h"
                COMMAND  ${Protobuf_PROTOC_EXECUTABLE}
                ARGS --grpc_out=${CMAKE_CURRENT_BINARY_DIR}
                -I ${PROTO_PATH} ${ABS_FIL}
                --plugin=protoc-gen-grpc=${GRPC_CPP_PLUGIN}
                DEPENDS ${ABS_FIL} ${Protobuf_PROTOC_EXECUTABLE}
                COMMENT "Running GRPC C++ compiler on ${FIL}"
                VERBATIM)
    endforeach()

    set_source_files_properties(${${SRCS}} ${${HDRS}} PROPERTIES GENERATED TRUE)
    set(${SRCS} ${${SRCS}} PARENT_SCOPE)
    set(${HDRS} ${${HDRS}} PARENT_SCOPE)
endfunction()