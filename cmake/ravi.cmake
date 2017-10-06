include_directories("${PROJECT_SOURCE_DIR}/ravi/include")
file(GLOB RAVI_HEADERS "${PROJECT_SOURCE_DIR}/ravi/include/*.h")
set(LUA_HEADERS 
  ${PROJECT_SOURCE_DIR}/ravi/include/lua.h 
  ${PROJECT_SOURCE_DIR}/ravi/include/luaconf.h 
  ${PROJECT_SOURCE_DIR}/ravi/include/lualib.h 
  ${PROJECT_SOURCE_DIR}/ravi/include/lauxlib.h)

add_definitions(-DLUA_COMPAT_5_2)
add_definitions(-DLUA_COMPAT_MODULE)
if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    message(STATUS "Enabling Lua extended test harness 'ltests'")
    add_definitions(-DLUA_USER_H="ltests.h")
endif ()
if (APPLE)
  add_definitions(-DLUA_USE_MACOSX)
else()
	if (NOT WIN32)
		add_definitions(-DLUA_USE_LINUX)
	endif()
endif()

# define the lua core source files
set(RAVI_CORE_SRCS ravi/src/lapi.c ravi/src/lcode.c ravi/src/lctype.c ravi/src/ldebug.c ravi/src/ldo.c ravi/src/ldump.c
        ravi/src/lfunc.c ravi/src/lgc.c ravi/src/llex.c ravi/src/lmem.c ravi/src/lobject.c ravi/src/lopcodes.c
        ravi/src/lparser.c ravi/src/lstate.c ravi/src/lstring.c ravi/src/ltable.c ravi/src/ltm.c ravi/src/lundump.c
        ravi/src/lvm.c ravi/src/lzio.c ravi/src/ravijit.cpp ravi/src/ltests.c)
# define the lua lib source files
set(RAVI_LIB_SRCS ravi/src/lauxlib.c ravi/src/lbaselib.c ravi/src/lbitlib.c ravi/src/lcorolib.c ravi/src/ldblib.c ravi/src/liolib.c
        ravi/src/lmathlib.c ravi/src/loslib.c ravi/src/ltablib.c ravi/src/lstrlib.c ravi/src/loadlib.c ravi/src/linit.c ravi/src/lutf8lib.c)
set(RAVI_CMD_SRCS ravi/src/lua.c)
set(RAVI_JIT_SRCS ravi/src/ravi_nojit.c)
set(RAVI_DEBUGGER_SRCS ravi/debugger/src/ravidebug.c ravi/debugger/src/json.c ravi/debugger/src/protocol.c)
set(RAVI_DEBUGGER_HEADERS ravi/debugger/src/json.h ravi/debugger/src/protocol.h)

if (MSVC OR APPLE)
  source_group("Ravi Headers" FILES ${RAVI_HEADERS})
  source_group("Ravi Sources" FILES ${RAVI_CORE_SRCS} ${RAVI_LIB_SRCS} ${RAVI_JIT_SRCS})
  source_group("Ravi Debugger Sources" FILES ${RAVI_DEBUGGER_SRCS})
  source_group("Ravi Debugger Headers" FILES ${RAVI_DEBUGGER_HEADERS})
endif()
