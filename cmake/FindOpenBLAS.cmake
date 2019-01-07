find_path(OPENBLAS_INCLUDE_DIR openblas_config.h
  PATHS
  c:/Software/OpenRedukti/include/openblas
  ~/Software/OpenRedukti/include/openblas
)

find_library(OPENBLAS_LIBRARY
  NAMES openblas libopenblas
  PATHS
  c:/Software/OpenRedukti/lib
  ~/Software/OpenRedukti/lib
)

find_library(LAPACK_LIBRARY
  NAMES lapack liblapack
  PATHS
  c:/Software/OpenRedukti/lib
  ~/Software/OpenRedukti/lib
)


set( OPENBLAS_INCLUDE_DIRS "${OPENBLAS_INCLUDE_DIR}" )
set( OPENBLAS_LIBRARIES "${OPENBLAS_LIBRARY}" )
set( LAPACK_LIBRARIES "${LAPACK_LIBRARY}" )
