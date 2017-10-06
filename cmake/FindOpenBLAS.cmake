find_path(OPENBLAS_INCLUDE_DIR openblas_config.h
  PATHS
  c:/OpenRedukti/include/openblas
  ~/OpenRedukti/include/openblas
)

find_library(OPENBLAS_LIBRARY
  NAMES openblas libopenblas
  PATHS
  c:/OpenRedukti/lib
  ~/OpenRedukti/lib
)

find_library(LAPACK_LIBRARY
  NAMES lapack liblapack
  PATHS
  c:/OpenRedukti/lib
  ~/OpenRedukti/lib
)


set( OPENBLAS_INCLUDE_DIRS "${OPENBLAS_INCLUDE_DIR}" )
set( OPENBLAS_LIBRARIES "${OPENBLAS_LIBRARY}" )
set( LAPACK_LIBRARIES "${LAPACK_LIBRARY}" )
