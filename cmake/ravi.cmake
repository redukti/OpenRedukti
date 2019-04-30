include_directories("${PROJECT_SOURCE_DIR}/ravi/include")
include_directories("${PROJECT_SOURCE_DIR}/ravi/dmr_c/src")
file(GLOB RAVI_HEADERS "${PROJECT_SOURCE_DIR}/ravi/include/*.h")
set(LUA_HEADERS 
  ${PROJECT_SOURCE_DIR}/ravi/include/lua.h 
  ${PROJECT_SOURCE_DIR}/ravi/include/luaconf.h 
  ${PROJECT_SOURCE_DIR}/ravi/include/lualib.h 
  ${PROJECT_SOURCE_DIR}/ravi/include/lauxlib.h)

add_definitions(-DLUA_COMPAT_5_2)
add_definitions(-DLUA_COMPAT_MODULE)
add_definitions(-DNO_LUA_DEBUG)

option(COMPUTED_GOTO "Controls whether the interpreter switch will use computed gotos on gcc/clang, default is ON" ON)

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
        ravi/src/lvm.c ravi/src/lzio.c ravi/src/ravijit.cpp ravi/src/ltests.c ravi/src/ravi_profile.c 
        ravi/src/ravi_ast.c ravi/src/ravi_membuf.c ravi/src/ravi_jitshared.c ravi/src/bit.c ravi/src/ravi_alloc.c)
if (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set_source_files_properties(ravi/src/lvm.c PROPERTIES COMPILE_FLAGS -DRAVI_USE_COMPUTED_GOTO)
elseif(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set_source_files_properties(ravi/src/lvm.c PROPERTIES COMPILE_FLAGS "-fno-crossjumping -fno-gcse -DRAVI_USE_COMPUTED_GOTO")
endif()
# define the lua lib source files
set(RAVI_LIB_SRCS ravi/src/lauxlib.c ravi/src/lbaselib.c ravi/src/lbitlib.c ravi/src/lcorolib.c ravi/src/ldblib.c ravi/src/liolib.c
        ravi/src/lmathlib.c ravi/src/loslib.c ravi/src/ltablib.c ravi/src/lstrlib.c ravi/src/loadlib.c ravi/src/linit.c ravi/src/lutf8lib.c)
set(RAVI_CMD_SRCS ravi/src/lua.c)
set(RAVI_JIT_SRCS ravi/src/ravi_nojit.c)
set(RAVI_DEBUGGER_SRCS ravi/debugger/src/ravidebug.c ravi/debugger/src/json.c ravi/debugger/src/protocol.c)
set(RAVI_DEBUGGER_HEADERS ravi/debugger/src/json.h ravi/debugger/src/protocol.h)

if (CMAKE_COMPILER_IS_GNUCC)
    execute_process(COMMAND ${CMAKE_C_COMPILER} --print-file-name=
        OUTPUT_VARIABLE GCC_BASE OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -print-multiarch
        OUTPUT_VARIABLE MULTIARCH_TRIPLET ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

    add_definitions(-DGCC_BASE="${GCC_BASE}")
    add_definitions(-DMULTIARCH_TRIPLET="${MULTIARCH_TRIPLET}")
endif()

message( STATUS "GCC_BASE_DIR      : " ${GCC_BASE})
message( STATUS "MULTIARCH_TRIPLET : " ${MULTIARCH_TRIPLET} )

set(DMR_C_HEADERS 
    ravi/dmr_c/src/allocate.h
    ravi/dmr_c/src/char.h
    ravi/dmr_c/src/expression.h
    ravi/dmr_c/src/flow.h
    ravi/dmr_c/src/ident-list.h
    ravi/dmr_c/src/linearize.h
    ravi/dmr_c/src/lib.h        
    ravi/dmr_c/src/parse.h
    ravi/dmr_c/src/port.h 
    ravi/dmr_c/src/ptrlist.h
    ravi/dmr_c/src/scope.h
    ravi/dmr_c/src/symbol.h
    ravi/dmr_c/src/target.h
    ravi/dmr_c/src/token.h 
    ravi/dmr_c/src/walksymbol.h         
    )

set(DMR_C_SRCS 
    ravi/dmr_c/src/allocate.c
    ravi/dmr_c/src/builtin.c
    ravi/dmr_c/src/char.c
    ravi/dmr_c/src/expression.c
    ravi/dmr_c/src/evaluate.c
    ravi/dmr_c/src/expand.c
    ravi/dmr_c/src/inline.c
    ravi/dmr_c/src/lib.c
    ravi/dmr_c/src/linearize.c
    ravi/dmr_c/src/liveness.c
    ravi/dmr_c/src/parse.c
    ravi/dmr_c/src/target.c
    ravi/dmr_c/src/tokenize.c
    ravi/dmr_c/src/pre-process.c
    ravi/dmr_c/src/ptrlist.c
    ravi/dmr_c/src/scope.c
    ravi/dmr_c/src/show-parse.c
    ravi/dmr_c/src/symbol.c
    ravi/dmr_c/src/walksymbol.c
    ravi/src/ravi_dmrc_parsesymbols.c
    )

set(DMR_C_JIT_HEADERS
    ravi/dmr_c/null-backend/dmr_c.h            
    )
include_directories("${PROJECT_SOURCE_DIR}/ravi/dmr_c/null-backend")    

if (MSVC OR APPLE)
  source_group("Ravi Headers" FILES ${RAVI_HEADERS})
  source_group("Ravi Sources" FILES ${RAVI_CORE_SRCS} ${RAVI_LIB_SRCS} ${RAVI_JIT_SRCS})
  source_group("Ravi Debugger Sources" FILES ${RAVI_DEBUGGER_SRCS})
  source_group("Ravi Debugger Headers" FILES ${RAVI_DEBUGGER_HEADERS})
  source_group("dmrC Headers" FILES ${DMR_C_HEADERS} ${DMR_C_JIT_HEADERS})
  source_group("dmrC Source Files" FILES ${DMR_C_SRCS} ${DMR_C_JIT_SRCS})
endif()
