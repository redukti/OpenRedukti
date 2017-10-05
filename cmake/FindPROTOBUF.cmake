if (CMAKE_BUILD_TYPE STREQUAL "Debug")

  find_path(PROTOBUF_INCLUDE_DIR google/protobuf/message.h
    PATHS
    c:/d/protobuf32d/include
    ~/proto3/include
    ~/local/include
  )

  find_library(PROTOBUF_LIBRARY
    NAMES libprotobufd libprotobuf protobufd protobuf
    PATHS
    c:/d/protobuf32d/lib
    ~/proto3/lib
    ~/local/lib
  )

else()

  find_path(PROTOBUF_INCLUDE_DIR google/protobuf/message.h
    PATHS
    c:/d/protobuf32/include
    ~/proto3/include
    ~/local/include
  )

  find_library(PROTOBUF_LIBRARY
    NAMES libprotobuf protobuf
    PATHS
    c:/d/protobuf32/lib
    ~/proto3/lib
    ~/local/lib
  )

endif()

set( PROTOBUF_INCLUDE_DIRS "${PROTOBUF_INCLUDE_DIR}" )
set( PROTOBUF_LIBRARIES "${PROTOBUF_LIBRARY}" )
