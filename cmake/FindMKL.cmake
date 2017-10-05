find_path(MKL_INCLUDE_DIR mkl.h
  PATHS
  c:/d/Intel/compilers_and_libraries/windows/mkl/include
)

find_library(MKL_CORE_LIBRARY
  NAMES mkl_core
  PATHS
  C:/d/Intel/compilers_and_libraries/windows/mkl/lib/intel64
)

find_library(MKL_INTERFACE_LIBRARY 
	NAMES mkl_intel_lp64
	PATHS
	C:/d/Intel/compilers_and_libraries/windows/mkl/lib/intel64
)

find_library(MKL_THREADING_LIBRARY  
	NAMES mkl_sequential
	PATHS
	C:/d/Intel/compilers_and_libraries/windows/mkl/lib/intel64
)


set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
set(MKL_LIBRARIES ${MKL_CORE_LIBRARY} ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY})