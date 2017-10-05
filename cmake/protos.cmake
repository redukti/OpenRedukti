# Protobuf sources 
include_directories("${PROJECT_SOURCE_DIR}/generated/cpp")
file(GLOB GENERATED_HEADERS "${PROJECT_SOURCE_DIR}/generated/cpp/*.h")
file(GLOB GENERATED_SRCS "${PROJECT_SOURCE_DIR}/generated/cpp/*.cc")
if (MSVC)
  source_group("Proto Sources" FILES ${GENERATED_SRCS})
  source_group("Proto Headers" FILES ${GENERATED_HEADERS})
endif(MSVC)